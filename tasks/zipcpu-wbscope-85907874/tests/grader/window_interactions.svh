// No single-slot capture can accumulate more than the drain allowance on
// any channel: the next accepted beat can close it and each following clock
// can discharge at most one single-beat transaction. This applies regardless
// of how many samples or transactions the implementation can store.
function automatic integer one_slot_debt;
    integer retirements;
    begin
        retirements = (w_count < aw_count) ? w_count : aw_count;
        if (b_count < retirements) retirements = b_count;
        one_slot_debt = total_samples - retirements;
    end
endfunction

always @(negedge clk)
    if (SCENARIO==3 && active && ideal_stream && capacity==1 && !triggered
            && one_slot_debt() > 4096)
        $fatal(1,"one-slot write debt exceeds possible drain service after closure");

task automatic sustain_single_slot;
    integer prior_debt, debt, windows, flat_windows;
    begin
        // A finite trace cannot certify all future cycles. Sample intervals
        // covering two 16K periods, then extend a growing backlog. Growth is
        // a stimulus choice, never a failure criterion by itself.
        debt=one_slot_debt(); windows=0; flat_windows=0;
        while (flat_windows<2 && windows<32) begin
            prior_debt=debt;
            repeat (32768) @(negedge clk);
            debt=one_slot_debt();
            flat_windows = debt<=prior_debt ? flat_windows+1 : 0;
            windows=windows+1;
        end
    end
endtask

task automatic control_fields_case;
    integer saved;
    reg [31:0] stopped_word;
    begin
        saved=irq_count;
        stage_ring(32'h2500,7);
        write_reg(0,32'h7bf00002,15);
        feed(10); await_completion(saved);
        read_reg(0); stopped_word=last_read;
        // A control read-modify-write can include the readback status.
        write_reg(0,stopped_word,15);
        saved=irq_count;
        feed(7); read_reg(0);
        descriptor(allocation_base(32'h2500,7),0);
        write_reg(0,32'h77f00000,15);
        read_reg(0); read_reg(8); read_reg(12);
        feed(3); await_completion(saved);
        read_reg(4);
    end
endtask

task automatic count_to_stop(input integer n, input integer h);
    integer saved;
    begin
        saved=irq_count;
        descriptor(0,n);
        tv=1;
        write_reg(0,32'h08000000 | h,15);
        // With continuous data and ideal memory this newly opened capture
        // must take every beat, independently of the concurrent status AR.
        ideal_stream=1;
        fork
            begin
                repeat(n+1+h) @(negedge clk);
                tv=0; ideal_stream=0;
            end
            begin
                repeat(n-1) @(negedge clk);
                read_reg(0); read_reg(0);
            end
        join
        await_completion(saved);
        read_reg(4);
    end
endtask

task automatic counter_boundaries_case;
    integer bitpos, n;
    begin
        // Upper capacity bits must also govern a short live prefix. Setting a
        // high bit alongside low bits distinguishes it from a truncated N.
        for (bitpos=12;bitpos<A-SIZE;bitpos=bitpos+1) begin
            n=(1<<bitpos)+3;
            descriptor(0,n);
            write_reg(0,32'h08000000 | (bitpos<20 ? (1<<bitpos) : 20'hfffff),15);
            feed(BURST_LIMIT+4); read_reg(0);
        end
        count_to_stop(34,33);
        count_to_stop(smaller(2049,MAX_SLOTS/2+1),0);
        if (MAX_SLOTS<=131072)
            count_to_stop(MAX_SLOTS,MAX_SLOTS-1);
        else begin
            count_to_stop(131073,0);
            n=smaller(524290,MAX_SLOTS);
            count_to_stop(n,n-1);
        end
    end
endtask

task automatic backpressured_trigger_case;
    integer saved, offered;
    reg pulse_sent;
    begin
        saved=irq_count;
        stage_ring(32'h3200,7); write_reg(0,2,15);
        feed(7); read_reg(0);
        allow_w=0; tv=1; offered=0; pulse_sent=0;
        while (offered<4112 && !pulse_sent) begin
            @(negedge clk);
            offered=offered+1;
            if (tv && !tr) begin
                trig=1;
                @(negedge clk);
                trig=0;
                pulse_sent=1;
            end
        end
        allow_w=1;
        // Leave the queued word valid. Observe the status using the actual
        // handshake history even if the receiver needs time to resume.
        repeat(32) @(negedge clk);
        read_reg(0);
        trig=1;
        repeat(64) @(negedge clk);
        if (active) begin
            // A nonideal open capture has no acquisition-rate deadline.
            // Firmware can still close and drain the entire accepted prefix.
            write_reg(0,0,15);
        end else begin
            await_completion(saved);
            read_reg(4);
        end
    end
endtask

task automatic response_debt_case;
    integer replace_capture, start_count, saved, limit, observations;
    begin
        for (limit=32;limit<=64;limit=limit+32)
        for (replace_capture=0;replace_capture<=1;replace_capture=replace_capture+1) begin
            saved=irq_count;
            stage_ring(32'h3120,1); write_reg(0,0,15);
            stage_ring(32'h6220,1);
            allow_b=0; tv=1; trig=0;
            start_count=source_index;
            // Offer a controlled prefix with no write responses. Derive
            // the fence from the observed issues, including credit throttling.
            repeat(4*limit+64) begin
                trig=!replace_capture && source_index-start_count==limit-1;
                @(negedge clk);
                if (source_index-start_count==limit) tv=0;
            end
            if (replace_capture) begin
                command(0,32'h08000000,15,0,0);
                tv=1; trig=0;
                repeat(23) @(negedge clk);
                allow_b=1; finish_write();
                ideal_stream=1; repeat(2) @(negedge clk);
                tv=0; ideal_stream=0;
            end else begin
                allow_b=1;
                // A credit-limited writer can have accepted fewer than
                // the offered prefix. Resume with a finite transfer trace.
                observations=0;
                while (!natural_close && observations<4*limit+64) begin
                    tv=1; trig=source_index-start_count>=limit-1;
                    @(negedge clk);
                    observations=observations+1;
                    if (natural_close && !source_stalled) tv=0;
                end
                trig=0;
                if (!natural_close) begin
                    // ARM fences precisely what was accepted. Establish a
                    // stopped image under continuous ideal peers afterwards.
                    command(0,32'h08000000,15,0,0);
                    tv=1; finish_write();
                    ideal_stream=1; repeat(2) @(negedge clk);
                    tv=0; ideal_stream=0;
                end
            end
            await_completion(saved);
            read_reg(4);
        end
    end
endtask

task automatic closing_arm_case;
    integer h, saved;
    begin
        for (h=0;h<=2;h=h+2) begin
            stage_ring(32'h3300,7);
            write_reg(0,32'h08000000 | h,15);
            stage_ring(32'h5500,5);
            feed(7+h);
            // Present the closing beat and the two write components together.
            // A later command grant closes whichever prefix then exists.
            saved=irq_count;
            allow_b=0; tv=1;
            command(0,32'h08000001,15,0,0);
            repeat(7) @(negedge clk);
            allow_b=1; finish_write();
            feed(7); await_completion(saved);
            read_reg(4);
        end
    end
endtask

task automatic invalid_drain_case;
    integer saved;
    begin
        saved=irq_count; trig=0;
        stage_ring(32'h23f0,9);
        tv=1; ideal_stream=1;
        write_reg(0,2,15);
        // The staged base changes while the old allocation stays active.
        write_reg(8,3,15);
        repeat(9) @(negedge clk);
        trig=1; @(negedge clk); trig=0;
        @(negedge clk);
        allow_b=0;
        @(negedge clk);
        // BVALID on the closing edge was already driven from the preceding
        // write transfers. Only subsequent responses have been withheld.
        if (!natural_close || source_stalled)
            $fatal(1,"closing sample missed under uninterrupted ideal writes");
        tv=0; ideal_stream=0;
        command(0,32'h0c000002,15,4,0);
        repeat(11) @(negedge clk);
        allow_b=1; finish_write();
        await_completion(saved);
        read_reg(8); read_reg(12);
        repeat(9*LANES+1) read_reg(4);
    end
endtask

task automatic held_write_commands_case;
    integer direction, op, start_count;
    reg [3:0] addr;
    reg [31:0] value;
    begin
        for (direction=0;direction<2;direction=direction+1)
        for (op=0;op<4;op=op+1) begin
            stage_ring(32'h3190,1);
            trig=0; tv=1;
            write_reg(0,0,15); ideal_stream=1;
            // No B can be withheld until a burst has paired. Up to that
            // point all ideal acquisition premises hold.
            start_count=source_index;
            allow_b=0;
            while (paired_bursts==b_count) begin
                @(negedge clk);
                if (source_index-start_count>4096 && paired_bursts==b_count)
                    $fatal(1,"one-slot backlog already exceeds possible drain service");
            end
            ideal_stream=0;
            case (op)
              0: begin addr=8; value=allocation_base(32'h52a0,1); end
              1: begin addr=12; value=257; end
              2: begin addr=0; value=32'h84013579; end
              3: begin addr=0; value=32'h08000001; end
            endcase
            command(addr,value,15,4*direction,4*(1-direction));
            allow_b=1; finish_write();
            read_reg(8); read_reg(12); read_reg(0);
            if (!source_stalled) tv=0;
        end
        // Retire the accepted prefix under the same allocation.
        write_reg(0,0,15);
    end
endtask

task automatic open_gap_case;
    integer resume_stream;
    begin
        stage_ring(32'h7100,7);
        for (resume_stream=0;resume_stream<=1;resume_stream=resume_stream+1) begin
            // The accepted prefix alone cannot shorten an open burst.
            tv=1; ideal_stream=1; write_reg(0,0,15);
            @(negedge clk);
            tv=0; ideal_stream=0;
            repeat(257) @(negedge clk);
            if (resume_stream) begin
                tv=1;
                repeat(11) @(negedge clk);
                if (!source_stalled) tv=0;
                repeat(17) @(negedge clk);
            end
            write_reg(0,0,15);
        end
    end
endtask

// Start the response allowance when the host releases BREADY. The command
// and response drivers run independently when a staging request straddles ARM.
task automatic take_response;
    integer start;
    reg received;
    begin
        start=clocks; received=0; lbr=1;
        while (!received) begin
            @(posedge clk); received=lbv;
            @(negedge clk);
            if (!received && clocks-start>=4096)
                $fatal(1,"missing ordered Lite write response after 4096 ready clocks");
        end
        lbr=0;
    end
endtask

task automatic check_staged_capture;
    integer irqs, lanes_to_read;
    begin
        irqs=irq_count;
        // Stop the source at either side of the priming boundary, then
        // overwrite slots before the trigger. Use every word of the image.
        feed(capacity-sample_count-1); read_reg(0);
        feed(1); read_reg(0);
        feed(capacity+3);
        trig=1; feed(1); trig=0; read_reg(0);
        feed(holdoff); await_completion(irqs);
        lanes_to_read=capacity*LANES+2;
        repeat(lanes_to_read) read_reg(4);
    end
endtask

task automatic overlap_staging_trial(input integer at_handoff, input integer change_capacity);
    integer outgoing_epoch, tick, skew_aw, skew_w;
    reg [3:0] addr, strobes;
    reg [31:0] value;
    begin
        // A one-slot outgoing capture supplies actual unretired writes.
        // Reserve room for both possible capacities even on one-page builds.
        trig=0; stage_ring(32'h91430240,1);
        tv=1; write_reg(0,32'h04000000,15); ideal_stream=1;
        descriptor(allocation_base(32'h12342ff0,19),13);
        allow_b=0;
        // The first suppressed response changes the peer at the next sampling
        // edge. Until then every paired write already has its response offered.
        while (paired_bursts <= b_count + (bvalid ? 1 : 0)) begin
            @(negedge clk);
            if (paired_bursts <= b_count + (bvalid ? 1 : 0)
                    && sample_count - (b_count-first_burst) > 4096)
                $fatal(1,"one-slot backlog exceeds the available drain service");
        end
        ideal_stream=0;
        command(0,3,15,0,0);
        outgoing_epoch=epoch;
        addr=change_capacity ? 12 : 8;
        strobes=change_capacity ? 4'b0001 : 4'b0010;
        value=change_capacity ? 32'h9fde0013 : 32'hface0637;
        skew_aw=change_capacity ? 4 : 0;
        skew_w=change_capacity ? 0 : 4;

        fork
            begin
                // The first B handshake starts acquisition even when the
                // staging response is held. Sample the very next ideal edge.
                while (epoch==outgoing_epoch) @(negedge clk);
                @(negedge clk);
                tv=0; ideal_stream=0;
            end
            begin
                ideal_stream=1;
                if (!at_handoff) begin
                    // Keep real B debt until both staging components transfer.
                    // An AR offered alongside the later component can see either
                    // side of its update, according to the AR handshake.
                    fork
                        begin command(addr,value,strobes,skew_aw,skew_w); allow_b=1; end
                        begin repeat(4) @(negedge clk); read_reg(addr); end
                    join
                    read_reg(addr); read_reg(0);
                    repeat(5) @(negedge clk);
                    take_response();
                end else begin
                    allow_b=1;
                    // Align the later component's offer with BREADY. When the
                    // response is already valid, both may transfer together.
                    // Slower arbitration can complete either operation first.
                    for (tick=0; tick<64 && !lbv; tick=tick+1) @(negedge clk);
                    fork
                        command(addr,value,strobes,skew_aw,skew_w);
                        begin repeat(4) @(negedge clk); take_response(); end
                        begin repeat(4) @(negedge clk); read_reg(addr); end
                    join
                end
                repeat(5) @(negedge clk);
                take_response();
                read_reg(addr);
            end
        join
        check_staged_capture();
        // Firmware does not restage either field. This ARM uses the value
        // merged by the preceding overlapping write.
        write_reg(0,2,15);
        check_staged_capture();
    end
endtask

task automatic overlap_staging_case;
    integer alignment, field;
    begin
        for (alignment=0; alignment<2; alignment=alignment+1)
            for (field=0; field<2; field=field+1)
                overlap_staging_trial(alignment,field);
    end
endtask

task automatic stage_field(input integer addr, input [31:0] location,
                          input integer aw_skew, input integer w_skew);
    begin
        if (addr == 8)
            command(8,{16'hb4d2,location[15:0]},4'b0011,aw_skew,w_skew);
        else
            command(12,32'hfc720013,4'b0001,aw_skew,w_skew);
    end
endtask

task automatic dual_stage_trial(input integer first_offset, input integer at_handoff,
                                input integer invalid_pair);
    integer second_offset, outgoing_epoch, tick;
    reg second_accepted;
    reg [31:0] arm_base, future_base;
    begin
        // The intermediate allocations each fit either new field by itself.
        // In the rejection trial only the combination exceeds the address span.
        arm_base = A==12 ? BYTE_SPAN-24*W : BYTE_SPAN-4096-3*W;
        if (invalid_pair) future_base = BYTE_SPAN-16*W;
        else future_base = A==12 ? BYTE_SPAN/2-W : BYTE_SPAN-4096+3*W;
        second_offset = 20-first_offset;
        second_accepted=0;

        trig=0; stage_ring(32'h91430240,1);
        tv=1; write_reg(0,32'h04000000,15); ideal_stream=1;
        descriptor(arm_base,13);
        allow_b=0;
        while (paired_bursts <= b_count + (bvalid ? 1 : 0)) begin
            @(negedge clk);
            if (paired_bursts <= b_count + (bvalid ? 1 : 0)
                    && sample_count-(b_count-first_burst)>4096)
                $fatal(1,"one-slot write backlog cannot drain within the published allowance");
        end
        ideal_stream=0;
        command(0,3,15,0,0);
        outgoing_epoch=epoch;
        stage_field(first_offset,future_base,0,3);
        if (at_handoff) begin
            allow_b=1;
            for (tick=0; tick<64 && !lbv; tick=tick+1) @(negedge clk);
        end

        fork
            begin
                stage_field(second_offset,future_base,4,0);
                second_accepted=1;
            end
            begin
                // AR is offered with the later staging component. Both orders
                // of grant and a coincident pre-update snapshot are legal.
                repeat(4) @(negedge clk);
                read_reg(second_offset);
                while (!second_accepted) @(negedge clk);
                read_reg(8); read_reg(12); read_reg(0);
            end
            begin
                if (!at_handoff) begin
                    while (!second_accepted) @(negedge clk);
                    // Both commands must execute independently of the old B.
                    // Reads may stall: release DDR and Lite B independently.
                    repeat(18) @(negedge clk);
                    allow_b=1;
                    repeat(8) @(negedge clk);
                end else begin
                    // Present BREADY with the later AW offer. A design can
                    // grant the staging component on this edge or a later one.
                    repeat(4) @(negedge clk);
                end
                take_response();
                while (!second_accepted) @(negedge clk);
                repeat(5) @(negedge clk); take_response();
                repeat(3) @(negedge clk); take_response();
            end
            begin
                ideal_stream=1;
                while (epoch==outgoing_epoch) @(negedge clk);
                @(negedge clk);
                tv=0; ideal_stream=0;
            end
        join

        check_staged_capture();
        // The older ARM used its execution-time snapshot. The next ARM sees
        // both merged staging words. An invalid pair preserves the stopped
        // image and the advanced word cursor through rejection.
        write_reg(0,2,15);
        read_reg(8); read_reg(12); read_reg(0);
        if (invalid_pair) begin
            repeat(capacity*LANES+1) read_reg(4);
            check_quiet();
        end else check_staged_capture();
    end
endtask

task automatic dual_stage_case;
    begin
        dual_stage_trial(8,0,0);
        dual_stage_trial(12,0,1);
        dual_stage_trial(12,1,0);
    end
endtask

// Initial data requests may reserve a DDR transaction before granting Lite
// AR. Make the initial stall finite if necessary to obtain that first grant.
// Once granted, hold the chosen memory channel through the second Lite AR.
task automatic begin_held_data_read(input integer hold_address);
    integer elapsed;
    reg granted;
    begin
        elapsed=0; granted=0;
        lrr=1; lra=4; lrv_addr=1;
        allow_ar=!hold_address; allow_r=hold_address;
        while (!granted) begin
            if (elapsed==13) begin allow_ar=1; allow_r=1; end
            @(posedge clk); granted=lrr_addr;
            @(negedge clk);
            elapsed=elapsed+1;
            if (!granted && elapsed>=13+4096)
                $fatal(1,"first data AR did not transfer after releasing DDR");
        end
        lrr=0; lrv_addr=0; lra=idle_addr(4);
        allow_ar=!hold_address; allow_r=hold_address;
    end
endtask

// This deadline belongs to the second *offer*. Neither the earlier read
// response nor the availability of the DDR read channels restarts it.
task automatic offer_second_read(input [3:0] addr);
    integer elapsed;
    reg granted;
    begin
        elapsed=0; granted=0; lra=addr; lrv_addr=1;
        while (!granted && elapsed<4096) begin
            @(posedge clk); granted=lrr_addr;
            @(negedge clk); elapsed=elapsed+1;
        end
        if (!granted)
            $fatal(1,"second Lite AR at offset %0d was not accepted within 4096 clocks of its offer",addr);
        lrv_addr=0; lra=idle_addr(addr);
    end
endtask

task automatic receive_two_reads;
    integer elapsed, received;
    begin
        elapsed=0; received=0; lrr=1;
        while (received<2 && elapsed<4096) begin
            @(posedge clk);
            if (lrv) received=received+1;
            @(negedge clk); elapsed=elapsed+1;
        end
        if (received!=2) $fatal(1,"two accepted Lite reads did not return with ready peers");
        lrr=0;
    end
endtask

task automatic paired_read_handoff(input [3:0] second_offset, input integer hold_address);
    integer saved, j;
    reg [31:0] origin;
    begin
        saved=irq_count;
        origin=allocation_base(32'hca650140,7);
        descriptor(origin,7); write_reg(0,32'h08000003,15);
        feed(11); await_completion(saved);
        // Cursor now crosses both a 32-bit lane boundary and a beat boundary
        // on the two data ARs. The replacement writes over the same bytes.
        for (j=0;j<LANES-1;j=j+1) read_reg(4);
        descriptor(origin,5);
        begin_held_data_read(hold_address);
        offer_second_read(second_offset);
        command(0,32'h08000001,15,0,2);
        // Neither queued Lite response can transfer. Release DDR only after
        // both Lite addresses and the replacing command have been accepted.
        // The common fence uses the count of ALL owed fetches at this point.
        allow_ar=1; allow_r=1; ideal_stream=1; tv=1;
        finish_write();
        // Five samples prime, the next triggers, and one more closes. The
        // first of these sampling edges follows the actual ARM B handshake.
        repeat(7) @(negedge clk);
        tv=0; ideal_stream=0;
        receive_two_reads();
        await_completion(saved+1); read_reg(4); check_quiet();
    end
endtask

task automatic two_outstanding_reads_case;
    integer saved;
    begin
        // The second control AR must snapshot the outgoing STOPPED state.
        paired_read_handoff(0,1);
        // Two old words can require two DDR responses even with Lite RREADY
        // held low; the controller may dispatch those fetches in either order.
        paired_read_handoff(4,0);
        paired_read_handoff(4,1);

        saved=irq_count;
        stage_ring(32'h34680100,7);
        write_reg(0,32'h08000002,15);
        feed(8);
        // The first control AR sees the triggered capture. Its Lite R is
        // held across natural completion. Wait through the published drain
        // allowance before using the stopped data port for the second read.
        lrr=1; read_request(0); lrr=0;
        feed(2);
        repeat(4096) @(negedge clk);
        allow_ar=0;
        offer_second_read(4);
        allow_ar=1;
        receive_two_reads();
        await_completion(saved); check_quiet();
    end
endtask
