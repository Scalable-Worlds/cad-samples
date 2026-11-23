// The host drives changes on falling edges. Its two request channels have
// independent schedules and VALIDs stay asserted until accepted.
function automatic [D-1:0] sample_word(input integer index);
    integer lane, group_index, phase, step, bit_to_flip;
    reg [31:0] v;
    begin
        group_index = index / 16;
        phase = index % 16;
        sample_word = 0;
        if (phase >= 2) begin
            for (lane=0; lane<LANES; lane=lane+1) begin
                // Independent lane sequences, with adjacent identical full
                // words and a pair differing in a single rotating data bit.
                step = phase<=4 ? 2 : phase;
                if (phase>=8) step = (phase-8) >> lane;
                v = (group_index*16 + step) + (lane+1)*32'h9e3779b9;
                v = (v ^ (v >> 16)) * 32'h85ebca6b;
                v = (v ^ (v >> 13)) * 32'hc2b2ae35;
                sample_word[32*lane +: 32] = v ^ (v >> 16);
            end
            if (phase>=12 && LANES>1)
                sample_word[31:0] = 0;
            if (phase==4) begin
                bit_to_flip = (group_index*37+32) % D;
                sample_word[bit_to_flip] = !sample_word[bit_to_flip];
            end
        end
    end
endfunction
integer source_index = 0;
assign td = tv ? sample_word(source_index) : {LANES{pattern ^ 32'h5d19ab47}};
always @(posedge clk)
    if (resetn && tv && tr) source_index <= source_index + 1;

// ACTIVE requests keep their address, data and strobes until the sampled edge.
// Inactive buses carry unrelated deterministic values, including on skewed writes.
function automatic [3:0] idle_addr(input [3:0] addr);
    idle_addr = addr ^ (pattern[3:0] | 4'h4);
endfunction

integer write_start, write_requests_done, read_start, read_address_done;
reg wait_lite_valid = 0;
reg memory_available = 1;
always @(posedge clk)
    memory_available <= allow_aw && allow_w && allow_b && allow_ar && allow_r && !jitter;

// A receiver that first observes VALID introduces a real stall. Its first
// transfer releases that stall and starts the ready-peer allowance again.
reg held_peer_aw=0, held_peer_w=0, held_peer_ar=0, released_peer=0;
always @(posedge clk) begin
    released_peer <= dest_mode!=0 && ((held_peer_aw && awvalid && awready)
                  || (held_peer_w && wvalid && wready) || (held_peer_ar && arvalid && arready));
    held_peer_aw <= awvalid && !awready;
    held_peer_w <= wvalid && !wready;
    held_peer_ar <= arvalid && !arready;
end

task automatic command(input [3:0] addr, input [31:0] value,
                       input [3:0] mask, input integer aw_delay, input integer w_delay);
    integer elapsed, skew;
    reg got_a, got_w;
    begin
        got_a=0; got_w=0; elapsed=0;
        skew = aw_delay>w_delay ? aw_delay : w_delay;
        write_start = clocks+skew;
        while (!got_a || !got_w) begin
            lav = !got_a && elapsed>=aw_delay;
            lwv = !got_w && elapsed>=w_delay;
            la = lav ? addr : idle_addr(addr);
            lwd = lwv ? value : value ^ (pattern|32'h80000001);
            lws = lwv ? mask : mask ^ (pattern[3:0]|4'h5);
            @(posedge clk);
            if (lav && lar) got_a=1;
            if (lwv && lwr) got_w=1;
            elapsed=elapsed+1;
            @(negedge clk);
            if (!memory_available && elapsed>=skew) write_start=clocks;
            if (elapsed>=4096+skew && (!got_a || !got_w))
                $fatal(1,"Lite AW/W did not transfer independently of DDR at offset %0d",addr);
        end
        lav=0; lwv=0; la=idle_addr(addr);
        lwd=value ^ (pattern|32'h80000001); lws=mask ^ (pattern[3:0]|4'h5);
        write_requests_done=clocks;
    end
endtask

task automatic finish_write;
    reg done, offered;
    begin
        if (clocks>write_requests_done) write_start=clocks; // Deliberately held BREADY.
        done=0; offered=0;
        while (!done) begin
            lbr=!wait_lite_valid || offered;
            @(posedge clk); done=lbr && lbv===1'b1;
            if (!lbr && lbv===1'b1) offered=1;
            @(negedge clk);
            if (!memory_available || released_peer || (!lbr && offered)) write_start=clocks;
            if ((!done && clocks-write_start>=4096) || (done && clocks-write_start>4096))
                $fatal(1,"Lite write exceeded 4096 ready clocks from request to response");
        end
        lbr=0;
    end
endtask

task automatic write_reg(input [3:0] addr,input [31:0] value,input [3:0] mask);
    begin
        lbr=!wait_lite_valid;
        command(addr,value,mask,0,0);
        finish_write();
    end
endtask

task automatic read_request(input [3:0] addr);
    reg done;
    begin
        lra=addr; lrv_addr=1; done=0; read_start=clocks;
        while (!done) begin
            @(posedge clk); done=lrr_addr===1'b1;
            @(negedge clk);
            if (!memory_available || released_peer) read_start=clocks;
            if (!done && clocks-read_start>=4096) $fatal(1,"Lite AR did not transfer at offset %0d",addr);
        end
        lrv_addr=0; lra=idle_addr(addr);
        read_address_done=clocks;
    end
endtask

reg [31:0] last_read;
task automatic finish_read;
    reg done, offered;
    begin
        if (clocks>read_address_done) read_start=clocks; // Deliberately held RREADY.
        done=0; offered=0;
        while (!done) begin
            lrr=!wait_lite_valid || offered;
            @(posedge clk); done=lrr && lrv===1'b1;
            if (!lrr && lrv===1'b1) offered=1;
            if (done) last_read=lrd;
            @(negedge clk);
            if (!memory_available || released_peer || (!lrr && offered)) read_start=clocks;
            if ((!done && clocks-read_start>=4096) || (done && clocks-read_start>4096))
                $fatal(1,"Lite read exceeded 4096 ready clocks from request to response");
        end
        lrr=0;
    end
endtask

task automatic read_reg(input [3:0] addr);
    begin lrr=!wait_lite_valid; read_request(addr); finish_read(); end
endtask
task automatic descriptor(input [31:0] location,input [31:0] slots);
    begin write_reg(8,location,15); write_reg(12,slots,15); end
endtask

task automatic stage_ring(input [31:0] preferred, input integer slots);
    begin descriptor(allocation_base(preferred, slots), slots); end
endtask

// Interrupted acquisitions are advanced by transfers; only the uninterrupted
// ideal traces place a clock-by-clock requirement on the receiver.
task automatic feed(input integer count);
    integer target;
    begin
        target=source_index+count; tv=1;
        while (source_index<target) @(negedge clk);
        tv=0;
    end
endtask

task automatic await_completion(input integer old_irqs);
    begin
        // The monitor times the interrupt and STOPPED from closure or stall
        // release, and timestamps each status snapshot at its own AR handshake.
        while (irq_count==old_irqs) @(negedge clk);
        if (irq_count!=old_irqs+1) $fatal(1,"Expected one completion interrupt after drain");
        read_reg(0);
        while (last_read[30:28]!=3) read_reg(0);
        check_quiet();
    end
endtask

task automatic staging_case;
    integer m;
    reg [31:0] pattern;
    begin
        read_reg(8); read_reg(12); read_reg(0);
        for (m=0; m<16; m=m+1) begin
            pattern=32'hcf2159b7 ^ (m*32'h72cf953b);
            command(8,pattern,m[3:0],m%4,(15-m)%4); finish_write();
            command(12,~pattern,(15-m),((15-m)%5),m%5); finish_write();
            read_reg(8); read_reg(12);
        end
        // Observe command execution before its B response is taken.
        command(8,32'h543210ee,15,5,0);
        fork
            read_reg(8);
            begin repeat(17) @(negedge clk); finish_write(); end
        join
        command(12,32'h457120fd,3,0,4);
        fork
            read_reg(12);
            begin repeat(17) @(negedge clk); finish_write(); end
        join
        write_reg(0,32'h8c0abcde,15);
        read_reg(0);
    end
endtask

task automatic try_descriptor(input [31:0] location,input [31:0] slots,input [19:0] h);
    begin descriptor(location,slots); write_reg(0,{12'b0,h},15); read_reg(8); read_reg(12); read_reg(0); end
endtask
task automatic descriptor_case;
    reg [31:0] edge_base;
    reg [19:0] max_h;
    integer bitpos;
    begin
        max_h=smaller(MAX_SLOTS-1,20'hfffff);
        edge_base=BYTE_SPAN-5*W;
        try_descriptor(0,0,0);
        for (bitpos=0;bitpos<SIZE;bitpos=bitpos+1)
            try_descriptor(W+(1<<bitpos),7,0);
        try_descriptor(edge_base,5,4);
        try_descriptor(edge_base,5,5);
        try_descriptor(edge_base,5,6);
        try_descriptor(edge_base,6,0);
        try_descriptor(0,MAX_SLOTS,max_h);
        try_descriptor(0,MAX_SLOTS+1,0);
        try_descriptor(0,32'hffffffff,0);
        try_descriptor(BYTE_SPAN-W,1,0);
        try_descriptor(BYTE_SPAN-W,32'h80000001,0);
        for (bitpos=A;bitpos<32;bitpos=bitpos+1)
            try_descriptor(32'd1<<bitpos,2,1);
        for (bitpos=A-$clog2(W);bitpos<32;bitpos=bitpos+1)
            try_descriptor(0,(32'd1<<bitpos)|1,0);
        try_descriptor(0,7,20'h80000);
        if (MAX_SLOTS >= 32'h00080001) begin
            try_descriptor(0,32'h00080000,20'h80000);
            try_descriptor(0,32'h00080001,20'h80000);
        end
        if (MAX_SLOTS >= 32'h00100000)
            try_descriptor(0,32'h00100000,20'hfffff);
    end
endtask

task automatic throughput_case(input integer slots,input [31:0] address,input integer h);
    integer saved, stop_at;
    begin
        saved=irq_count; descriptor(address,slots);
        // The full-rate write guarantee has no read-channel prerequisite.
        ideal_stream=1; tv=1; trig=0;
        write_reg(0,h,15);
        allow_ar=0; allow_r=0;
        if (slots==1) sustain_single_slot();
        else repeat (2*slots+29) @(negedge clk);
        trig=1;
        stop_at=source_index + 1 + h;
        repeat(1+h) @(negedge clk);
        if (source_index!=stop_at) $fatal(1,"The ideal stream lost a sampling edge");
        tv=0; trig=0; allow_ar=1; allow_r=1;
        await_completion(saved);
        read_reg(4);
    end
endtask

task automatic flags_case;
    integer saved;
    begin
        saved=irq_count; stage_ring(32'h2400,7);
        ideal_stream=1; tv=1; trig=1;
        write_reg(0,32'h04000003,15);
        repeat(11) @(negedge clk);
        read_reg(0);
        // External trigger is disabled. The flags write raises MANUAL while
        // leaving DISABLE set. Its own transfer still sees the old bits.
        tv=1; command(0,32'h8c0fffff,15,0,4);
        finish_write();
        await_completion(saved);
        saved=irq_count;
        // A still-asserted source resumes on the first edge of the next
        // acquisition. Disable the external trigger on its qualifying edge.
        trig=0; write_reg(0,2,15);
        repeat(10) @(negedge clk);
        la=0; lwd=32'h840aaaaa; lws=15; lav=1; lwv=1;
        // While the Lite write is pending, the old controls permit trigger.
        trig=1; command(0,32'h840aaaaa,15,0,0);
        finish_write(); await_completion(saved);
        read_reg(0);
        // Keep the next offered word stable across STOPPED and the new ARM.
        trig=0; ideal_stream=0;
        saved=irq_count;
        write_reg(0,2,15);
        feed(1);
        write_reg(0,32'h8c0fffff,15); read_reg(0);
        write_reg(0,32'h88012345,15); read_reg(0);
        write_reg(0,32'h80012345,15); read_reg(0);
        write_reg(0,32'h840fffff,15); read_reg(0);
        write_reg(0,32'h8c012345,15); read_reg(0);
        write_reg(0,32'h84012345,15); read_reg(0);
        trig=1; feed(9); read_reg(0);
        write_reg(0,32'h800abcde,15); read_reg(0);
        feed(3); trig=0; await_completion(saved);
    end
endtask

task automatic early_arm_case;
    integer saved;
    begin
        saved=irq_count; stage_ring(32'h12340fe0,smaller(409,MAX_SLOTS/2));
        ideal_stream=1; tv=1;
        write_reg(0,17,15);
        repeat(9) @(negedge clk);
        stage_ring(32'hc6791fc0,17);
        tv=1;
        command(0,32'h08000004,15,5,0);
        fork
            begin read_reg(0); end
            begin repeat(7) @(negedge clk); finish_write(); end
        join
        // Incoming MANUAL must prime from a fresh first index.
        await_completion(saved);
        read_reg(4);
    end
endtask

task automatic pending_writes_case;
    integer i;
    begin
        for (i=0;i<3;i=i+1) begin
            stage_ring(32'h3000+256*i,13);
            write_reg(0,9,15); trig=0;
            // Part of the final ring pass is accepted post-trigger.
            feed(13); trig=1; feed(2); trig=0;
            stage_ring(32'h3070+256*i,7);
            if (i==0) allow_aw=0;
            if (i==1) allow_w=0;
            if (i==2) allow_b=0;
            command(0,2,15,i,3-i);
            fork
                begin read_reg(0); end
                begin
                    repeat(23) @(negedge clk);
                    allow_aw=1; allow_w=1; allow_b=1;
                    finish_write();
                end
            join
            read_reg(0);
            if (irq_count!=0) $fatal(1,"Replacement of triggered capture raised an interrupt");
        end
    end
endtask

task automatic suppressed_irq_case;
    integer saved;
    begin
        saved=irq_count;
        stage_ring(32'h3ff0,9);
        write_reg(0,32'h08000002,15);
        stage_ring(32'h4010,5);
        feed(11);
        // Hold B before offering the closing sample. A receiver with no
        // spare credit may backpressure it, in which case ARM closes the
        // accepted prefix and the stalled word stays valid through handoff.
        allow_b=0; tv=1;
        repeat(17) begin
            @(negedge clk);
            if (sample_count==12 && !source_stalled) tv=0;
        end
        command(0,32'h08000001,15,0,0);
        allow_b=1;
        fork
            begin read_reg(0); end
            begin repeat(31) @(negedge clk); finish_write(); end
        join
        if (sample_count<7) feed(7-sample_count);
        await_completion(saved);
    end
endtask

task automatic readout_case;
    integer saved,j,ring;
    reg [31:0] address;
    begin
        for (ring=0;ring<2;ring=ring+1) begin
            saved=irq_count;
            address = ring==0 ? BYTE_SPAN-11*W
                              : allocation_base(A>16 ? 32'h12340000 : 32'h2930,11);
            descriptor(address,11);
            // Build a wrapped, nonzero-origin cursor from an independent
            // uninterrupted acquisition. Readout runs after both completions.
            trig=0; tv=1; ideal_stream=1;
            write_reg(0,10,15);
            repeat(29) @(negedge clk);
            trig=1; @(negedge clk); trig=0;
            repeat(10) @(negedge clk);
            tv=0; ideal_stream=0;
            await_completion(saved);
            for (j=0;j<2*11*LANES+3;j=j+1) read_reg(4);
            command(0,32'h8400cdef,15,0,5); finish_write();
            read_reg(0); read_reg(4);
            command(4,0,15,5,0); finish_write();
            for (j=0;j<11*LANES+1;j=j+1) read_reg(4);
            command(4,ring ? 32'hffffffff : 32'h51271234,15,0,5); finish_write();
            repeat(LANES+1) read_reg(4);
            read_reg(0);
            if (irq_count!=saved+1) $fatal(1,"stopped cursor or control writes raised an interrupt");
        end
        write_reg(0,0,15);
    end
endtask

// A stopped-data request can reserve a DDR fetch before the Lite address
// transfers. Observe both interfaces while offering ARM under a memory stall.
task automatic read_fence_case(input integer mode, input integer simultaneous);
    integer saved, j, old_ar, k;
    reg got_read, got_aw, got_w, arm_done;
    begin : trial
        saved=irq_count;
        stage_ring(32'h6000,7); trig=0;
        write_reg(0,32'h08000003,15);
        feed(11); await_completion(saved);
        for (j=0;j<LANES-1;j=j+1) read_reg(4);
        stage_ring(32'h6000+W,5);
        old_ar=ar_count; got_read=0; got_aw=0; got_w=0;
        allow_ar=mode!=0; allow_r=mode!=1;
        if (!simultaneous) begin
            lrv_addr=1; lra=4;
            // Wait for a genuine blocked DDR address or an accepted address
            // awaiting R. A delayed Lite address grant is allowed either way.
            k=0;
            while (!(mode==0 ? arvalid===1'b1 : ar_count>old_ar)) begin
                @(posedge clk);
                if (lrv_addr && lrr_addr) got_read=1;
                @(negedge clk);
                if (got_read) begin lrv_addr=0; lra=idle_addr(4); end
                k=k+1;
                if (k>4096) $fatal(1,"stopped read did not offer the required DDR fetch");
            end
        end
        // In the simultaneous case ARVALID is first offered alongside WVALID,
        // independently of whether the address has already been accepted.
        k=0;
        while (!got_aw || !got_w) begin
            lav=!got_aw; lwv=!got_w && k>=2;
            la=lav ? 0 : idle_addr(0);
            lwd=lwv ? 32'h08000001 : ~32'h08000001; lws=lwv ? 15 : 0;
            lrv_addr=!got_read && (!simultaneous || k>=2);
            lra=lrv_addr ? 4 : idle_addr(4);
            @(posedge clk);
            if (lrv_addr && lrr_addr) got_read=1;
            if (lav && lar) got_aw=1;
            if (lwv && lwr) got_w=1;
            @(negedge clk);
            k=k+1;
            if ((!got_aw || !got_w) && k>=4098)
                $fatal(1,"Lite ARM command waited for the old DDR read to retire");
        end
        lav=0; lwv=0; la=idle_addr(0); lwd=32'hf7fffffe; lws=0;
        // A read arbitrated behind ARM has no stopped snapshot. If a fetch
        // was already in flight, it must still retire before the write reply.
        if (!got_read) begin
            // The old data request lost arbitration. Its eventual payload is
            // unspecified, but ARM must still retire and enable its new ring.
            late_reserved_data=1;
            ideal_stream=1; tv=1;
            k=0; arm_done=0;
            while (!arm_done) begin
                if (k>=3) begin allow_ar=1; allow_r=1; lbr=1; end
                @(posedge clk);
                if (lrv_addr && lrr_addr) got_read=1;
                arm_done=lbr && lbv;
                @(negedge clk);
                if (got_read) begin lrv_addr=0; lra=idle_addr(4); end
                k=k+1;
                if (!arm_done && k>=4096+3) $fatal(1,"ARM did not finish after DDR read release");
            end
            lbr=0;
            while (irq_count==saved+1) begin
                @(posedge clk);
                if (lrv_addr && lrr_addr) got_read=1;
                @(negedge clk);
                if (got_read) begin lrv_addr=0; lra=idle_addr(4); end
                if (natural_close && !source_stalled) tv=0;
            end
            // Keep the request offered. The granting edge determines its
            // status. RREADY now releases any reply already being held.
            lrr=1; read_start=clocks; read_address_done=clocks;
            if (!got_read) read_request(4);
            finish_read();
            late_reserved_data=0;
            read_reg(8); read_reg(0);
            while (last_read[30:28]!=3) read_reg(0);
            // A read granted between owners has an unspecified cursor effect.
            // Rewind before checking the next stopped read.
            write_reg(4,0,15); read_reg(4);
            write_reg(0,0,15);
            disable trial;
        end
        lrv_addr=0; lra=idle_addr(4);
        repeat(3) @(negedge clk);
        allow_ar=1; allow_r=1;
        write_start=clocks; write_requests_done=clocks;
        finish_write();
        // Old Lite R stays stalled as the new capture overwrites shared DDR.
        feed(7);
        repeat(33) @(negedge clk);
        read_start=clocks; read_address_done=clocks;
        finish_read(); await_completion(saved+1);
        read_reg(4);
    end
endtask

task automatic snapshot_case;
    integer saved;
    begin
        // Each target read stays outstanding across an update. Concurrent
        // request channels can also transfer on a common clock edge.
        stage_ring(32'h1000,19);
        read_request(8);
        fork
            begin command(8,allocation_base(32'h20e0,23),15,4,0); finish_write(); end
            begin repeat(13) @(negedge clk); finish_read(); end
        join
        fork
            begin read_request(12); repeat(7) @(negedge clk); finish_read(); end
            begin command(12,23,15,0,0); finish_write(); end
        join
        write_reg(0,4,15);
        fork
            begin read_request(0); repeat(5) @(negedge clk); finish_read(); end
            begin command(0,32'h840bbbba,15,0,0); finish_write(); end
        join
        feed(4);
        // Preserve a control AR snapshot across ARM and fresh samples.
        read_request(0); command(0,32'h08000002,15,0,3);
        finish_write(); feed(3); finish_read();
        saved=irq_count; feed(23+1+2-3); await_completion(saved);
        fork
            begin read_request(0); repeat(7) @(negedge clk); finish_read(); end
            begin command(0,1,15,0,0); finish_write(); end
        join
    end
endtask

task automatic invalid_live_case;
    integer saved;
    begin
        saved=irq_count;
        stage_ring(32'h17e0,13); write_reg(0,4,15);
        feed(12);
        // A staged wrap in the descriptor arithmetic and H >= N both fail
        // without resetting the accepted prefix or changing active controls.
        descriptor(BYTE_SPAN-W,3);
        tv=1;
        command(0,32'h0c000000,15,0,5); finish_write();
        read_reg(0); feed(3); trig=1; feed(1); trig=0;
        stage_ring(32'h2100,1); write_reg(0,32'h0c000001,15);
        read_reg(8); read_reg(12); read_reg(0); feed(4); await_completion(saved);
        repeat(3*LANES+1) read_reg(4);
        descriptor(3,7);
        write_reg(0,32'h0c000006,15);
        read_reg(8); read_reg(12); read_reg(0);
        repeat(13*LANES) read_reg(4);
        if (irq_count!=saved+1) $fatal(1,"Invalid ARM changed completion interrupt history");
    end
endtask

task automatic valid_dependency_case;
    integer mode, sink, saved;
    begin
        wait_lite_valid=1;
        for (mode=1;mode<=3;mode=mode+1) begin
            dest_mode=0;
            stage_ring(32'h7100,7); write_reg(0,0,15);
            feed(1);
            // Closing forces the final one-beat burst. A legal receiver may
            // await address or data VALID before accepting its first transfer.
            dest_mode=mode;
            write_reg(0,32'h08000000,15);
            dest_mode=0;
            saved=irq_count;
            feed(8); await_completion(saved);
            dest_mode=mode;
            read_reg(4);
        end
        // One sink accepts W through WLAST before granting its AW. The
        // other waits for the matching AW before accepting any of its W.
        for (sink=4;sink<=5;sink=sink+1) begin
            dest_mode=0; stage_ring(32'h7180,7); write_reg(0,0,15);
            feed(4);
            dest_mode=sink; write_reg(0,0,15);
            for (mode=0; mode<3; mode=mode+1) begin
                tv=1;
                repeat(35) @(negedge clk);
                if (!source_stalled) tv=0;
                repeat(21) @(negedge clk);
            end
            write_reg(0,0,15);
        end
    end
endtask

task automatic live_count_case;
    integer saved, slots;
    begin
        saved=irq_count;
        // The full address span must also work as a live capture, including
        // its unprimed status, first bursts and short outgoing prefix.
        descriptor(0,MAX_SLOTS);
        write_reg(0,32'h08000000,15);
        feed(smaller(BURST_LIMIT+7,MAX_SLOTS-1)); read_reg(0);
        slots = smaller(1025,MAX_SLOTS/2+1);
        descriptor(0,slots); write_reg(0,32'h08000010,15);
        feed(slots-1); read_reg(0);
        feed(1); read_reg(0);
        feed(1); read_reg(0);
        feed(15); read_reg(0);
        feed(1); await_completion(saved);
        repeat(2*LANES+1) read_reg(4);
    end
endtask

task automatic edge_status_case;
    integer saved;
    begin
        saved=irq_count; stage_ring(32'h1280,7);
        write_reg(0,2,15);
        // Trigger is ignored both on the first beat and on the priming beat.
        trig=1; feed(1); trig=0; feed(5); read_reg(0);
        trig=1;
        fork feed(1); read_reg(0); join
        trig=0; read_reg(0);
        // An external pulse with no accepted sample has no lasting effect.
        trig=1; @(negedge clk); trig=0; feed(1); read_reg(0);
        trig=1;
        fork feed(1); read_reg(0); join
        trig=0; read_reg(0); feed(1);
        fork feed(1); read_reg(0); join
        await_completion(saved);
        stage_ring(32'h20b0,9);
        command(0,32'h04000005,15,0,0);
        fork
            begin read_reg(8); read_reg(12); read_reg(0); end
            begin repeat(23) @(negedge clk); finish_write(); end
        join
        command(0,32'h08000007,15,0,0);
        // An AR may coincide with B handoff or be arbitrated to a later edge.
        fork read_reg(0); finish_write(); join
    end
endtask

initial begin
    resetn=0; trig=0; tv=0; lav=0; la=0; lwv=0; lwd=0; lws=0;
    lbr=0; lrv_addr=0; lra=0; lrr=0;
    repeat(5) @(negedge clk);
    resetn=1;
    @(negedge clk);
    begin : directed_traffic
        case (SCENARIO)
            1: staging_case();
            2: descriptor_case();
            3: throughput_case(1,BYTE_SPAN-W,0);
            // One-page builds wrap the entire page; larger builds cross
            // an interior page boundary three beats after the opening AW.
            4: if (A==12) throughput_case(MAX_SLOTS,0,5);
               else throughput_case(smaller(2*BURST_LIMIT+9,(BYTE_SPAN-4096)/W+3),
                                    4096-3*W,5);
            5: flags_case();
            6: early_arm_case();
            7: pending_writes_case();
            8: suppressed_irq_case();
            9: readout_case();
            10: read_fence_case(0,0);
            11: snapshot_case();
            12: invalid_live_case();
            13: mixed_stalls_case();
            14: valid_dependency_case();
            15: live_count_case();
            16: edge_status_case();
            17: read_fence_case(0,1);
            18: read_fence_case(1,1);
            19: read_fence_case(1,0);
            20: control_fields_case();
            21: counter_boundaries_case();
            22: backpressured_trigger_case();
            23: response_debt_case();
            24: closing_arm_case();
            25: invalid_drain_case();
            26: held_write_commands_case();
            27: open_gap_case();
            28: overlap_staging_case();
            29: dual_stage_case();
            30: two_outstanding_reads_case();
            default: $fatal(1,"Unknown stimulus");
        endcase
    end
    $display("BENCH_OK"); $finish;
end
