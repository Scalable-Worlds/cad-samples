// Sustained backpressure combines long channel-specific holds with per-cycle
// jitter. The source changes TVALID only when idle or following a transfer.
task automatic mixed_stalls_case;
    integer slots, trigger_at, stop_at, previous, tick, saved, h, word_index;
    begin
        slots = BURST_LIMIT/4 + 19;
        h = 5;
        trigger_at = ((1 << G)/slots + 8)*slots + 5;
        stop_at = trigger_at + 1 + h;
        saved = irq_count;
        stage_ring(32'h5000 - 3*W, slots);
        write_reg(0,h,15);
        previous = source_index;
        pattern = 32'h768af31d ^ (A<<23) ^ (D<<11) ^ G;
        jitter = 1;
        for (tick=0; tick < 9*stop_at && source_index<stop_at; tick=tick+1) begin
            allow_aw = (tick%113 >= 39);
            allow_w = (tick%79 >= 17);
            allow_b = (tick%193 >= 71);
            trig = source_index >= trigger_at;
            if (!tv || source_index != previous)
                tv = (tick%53 >= 9) && (pattern[4] || pattern[27]);
            previous = source_index;
            @(negedge clk);
        end
        // Release the peers at the end of the finite traffic sequence.
        // If the capture remains open, ARM closes the prefix that actually
        // transferred. A stalled stream word remains offered through handoff.
        jitter=0; allow_aw=1; allow_w=1; allow_b=1;
        trig=0;
        if (!source_stalled) tv=0;
        if (natural_close) begin
            await_completion(saved);
            for (word_index=0; word_index<slots*LANES+2; word_index=word_index+1)
                read_reg(4);
        end
        // A final ARM allows a directly forwarded Lite response to outlive
        // its host handshake while still requiring the DDR transaction to retire.
        write_reg(0,0,15);
    end
endtask
