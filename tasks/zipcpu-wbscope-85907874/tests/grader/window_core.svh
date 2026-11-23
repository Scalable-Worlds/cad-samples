localparam W = D/8;
localparam SIZE = $clog2(W);
localparam LANES = D/32;
localparam BURST_LIMIT = (1 << (G-1)) < 256 ? (1 << (G-1)) : 256;

localparam [63:0] BYTE_SPAN = 64'd1 << A;
localparam integer MAX_SLOTS = BYTE_SPAN / W;

function automatic integer smaller(input integer a, input integer b);
    smaller = a < b ? a : b;
endfunction

// Size every live allocation in samples, then place its full byte interval.
// Keep the preferred address bits that exist in this build. At the address
// ceiling the final byte of the allocation can coincide with the final byte
// of DDR. Invalid descriptors are driven separately, as raw register words.
function automatic [31:0] allocation_base(input [31:0] preferred, input integer slots);
    reg [63:0] address, last_start;
    begin
        address = preferred[A-1:0];
        last_start = BYTE_SPAN - slots * W;
        allocation_base = address <= last_start ? address : last_start;
    end
endfunction

reg clk = 0, resetn = 0;
always #5 clk = !clk;
reg trig = 0, tv = 0;
wire [D-1:0] td;
wire tr, intr;

reg lav = 0, lwv = 0, lbr = 0, lrv_addr = 0, lrr = 0;
reg [3:0] la = 0, lra = 0, lws = 0;
reg [31:0] lwd = 0;
reg [1:0] host_awprot = 2, host_arprot = 3;
always @(posedge clk) begin
    if (lav && lar) host_awprot <= host_awprot+1'b1;
    if (lrv_addr && lrr_addr) host_arprot <= host_arprot+1'b1;
end
wire lar, lwr, lbv, lrr_addr, lrv;
wire [1:0] lbresp, lrresp;
wire [31:0] lrd;

reg allow_aw = 1, allow_w = 1, allow_b = 1, allow_ar = 1, allow_r = 1;
reg jitter = 0, ideal_stream = 0;
reg [31:0] pattern = 32'h3a1c57b9;

wire awvalid, awid, awlock, wvalid, wlast, bready;
wire arvalid, arid, arlock, rready;
wire [A-1:0] awaddr, araddr;
wire [7:0] awlen, arlen;
wire [2:0] awsize, arsize, awprot, arprot;
wire [1:0] awburst, arburst;
wire [3:0] awcache, arcache, awqos, arqos;
wire [W-1:0] wstrb;
wire [D-1:0] wdata;
reg bvalid = 0, bid = 0, rvalid = 0, rid = 0;
reg [D-1:0] rdata = 0;
wire [1:0] bresp = bvalid ? 2'b00 : pattern[4:3];
wire [1:0] rresp = rvalid ? 2'b00 : pattern[19:18];
wire response_last = rvalid || pattern[27];

// These destinations can wait for a previously sampled VALID. In each
// cross-channel mode, record the first offers so that an address or data
// transfer that has already finished still qualifies its partner.
integer dest_mode = 0, last_dest_mode = 0;
reg saw_aw = 0, saw_w = 0, saw_ar = 0;
reg offered_aw = 0, offered_w = 0;
always @(posedge clk) begin
    saw_aw <= awvalid;
    saw_w <= wvalid;
    saw_ar <= arvalid;
    offered_aw <= awvalid || (dest_mode == last_dest_mode && offered_aw);
    offered_w <= wvalid || (dest_mode == last_dest_mode && offered_w);
    last_dest_mode <= dest_mode;
end
// The receiver may accept every data beat before granting that burst's AW.
integer granted_addresses = 0, received_wlasts = 0;
always @(posedge clk) begin
    if (awvalid && awready) granted_addresses <= granted_addresses + 1;
    if (wvalid && wready && wlast) received_wlasts <= received_wlasts + 1;
end
wire awready = allow_aw && (!jitter || pattern[2] || pattern[7])
             && (!dest_mode || (saw_aw && last_dest_mode == dest_mode
                               && (dest_mode != 2 || offered_w)
                               && (dest_mode != 4 || received_wlasts > granted_addresses)));
wire wready = allow_w && (!jitter || pattern[5])
            && (!dest_mode || (saw_w && last_dest_mode == dest_mode
                              && (dest_mode != 3 || offered_aw)
                              && (dest_mode != 5 || granted_addresses > received_wlasts)));
wire arready = allow_ar && (!jitter || pattern[11]) && (!dest_mode || saw_ar);

memscope
`ifndef FIXED_WINDOW_NETLIST
         #(.OPT_WINDOWED(1), .OPT_TREADY_WHILE_IDLE(0),
           .C_AXI_ADDR_WIDTH(A), .C_AXI_DATA_WIDTH(D), .LGFIFO(G))
`endif
         device (
    .S_AXI_ACLK(clk), .S_AXI_ARESETN(resetn),
    .i_trigger(trig), .S_AXIS_TVALID(tv), .S_AXIS_TREADY(tr), .S_AXIS_TDATA(td),
    .S_AXIL_AWVALID(lav), .S_AXIL_AWREADY(lar), .S_AXIL_AWADDR(la), .S_AXIL_AWPROT({1'b0,host_awprot}),
    .S_AXIL_WVALID(lwv), .S_AXIL_WREADY(lwr), .S_AXIL_WDATA(lwd), .S_AXIL_WSTRB(lws),
    .S_AXIL_BVALID(lbv), .S_AXIL_BREADY(lbr), .S_AXIL_BRESP(lbresp),
    .S_AXIL_ARVALID(lrv_addr), .S_AXIL_ARREADY(lrr_addr), .S_AXIL_ARADDR(lra), .S_AXIL_ARPROT({1'b0,host_arprot}),
    .S_AXIL_RVALID(lrv), .S_AXIL_RREADY(lrr), .S_AXIL_RDATA(lrd), .S_AXIL_RRESP(lrresp),
    .M_AXI_AWVALID(awvalid), .M_AXI_AWREADY(awready), .M_AXI_AWID(awid),
    .M_AXI_AWADDR(awaddr), .M_AXI_AWLEN(awlen), .M_AXI_AWSIZE(awsize), .M_AXI_AWBURST(awburst),
    .M_AXI_AWLOCK(awlock), .M_AXI_AWCACHE(awcache), .M_AXI_AWPROT(awprot), .M_AXI_AWQOS(awqos),
    .M_AXI_WVALID(wvalid), .M_AXI_WREADY(wready), .M_AXI_WDATA(wdata),
    .M_AXI_WSTRB(wstrb), .M_AXI_WLAST(wlast),
    .M_AXI_BVALID(bvalid), .M_AXI_BREADY(bready), .M_AXI_BID(bid), .M_AXI_BRESP(bresp),
    .M_AXI_ARVALID(arvalid), .M_AXI_ARREADY(arready), .M_AXI_ARID(arid),
    .M_AXI_ARADDR(araddr), .M_AXI_ARLEN(arlen), .M_AXI_ARSIZE(arsize), .M_AXI_ARBURST(arburst),
    .M_AXI_ARLOCK(arlock), .M_AXI_ARCACHE(arcache), .M_AXI_ARPROT(arprot), .M_AXI_ARQOS(arqos),
    .M_AXI_RVALID(rvalid), .M_AXI_RREADY(rready), .M_AXI_RID(rid), .M_AXI_RDATA(rdata),
    .M_AXI_RLAST(response_last), .M_AXI_RRESP(rresp), .o_int(intr)
);

// Separate journals for the accepted stream, the two independent AXI write
// channels, and the physical DDR locations committed by pairing those channels.
reg [D-1:0] sample_log[$], write_log[$];
reg [A-1:0] written_addresses[$], burst_addresses[$];
integer burst_lengths[$];
reg write_lasts[$], write_ids[$];
reg [D-1:0] fetched_words[$];
reg read_ids[$];
reg [A-1:0] expected_fetches[$], observed_fetches[$];
integer matched_fetches = 0;

integer sample_count = 0, total_samples = 0;
integer aw_count = 0, w_count = 0, b_count = 0, ar_count = 0, r_count = 0, irq_count = 0;
integer paired_beats = 0, paired_bursts = 0, beat_in_burst = 0, addressed_beats = 0;
integer first_sample = 0, first_burst = 0;
reg [31:0] staged_b = 0, staged_n = 0, base = 0, capacity = 0;
reg [19:0] holdoff = 0;
reg manual = 0, disabled = 0;
reg [31:0] next_base, next_capacity, next_control;
reg active = 0, triggered = 0, pending_arm = 0;
reg natural_close = 0, sent_interrupt = 0, published_stop = 0;
integer tail_left = 0;
reg [63:0] word_cursor = 0;

reg have_address = 0, have_data = 0, owe_r = 0;
reg [3:0] command_address, command_mask;
reg [31:0] command_data;
// In command order, each entry records whether its B handshake installs
// a descriptor and which BRESP the master must see.
reg [2:0] replies[$], current_reply;
reg [31:0] expected_r;
reg pending_control_read;
reg [31:0] incoming_control_r;
reg stop_choice, stop_was_required;
reg late_reserved_data = 0, read_outside_stopped = 0, optional_fetch_claimed = 0;
integer epoch = 0, read_epoch, read_fetch_index = -1;
// Store the complete AR-time expectation when another address transfers
// before the first response. Queuing only the returned word would lose
// the epoch, fetch obligation, and permissible pending status choices.
typedef struct packed {
    reg [31:0] word_value, other_control;
    integer fetch_index, owner;
    reg pending_control, may_stop, must_stop, unspecified;
} lite_read_t;
lite_read_t oldest_read;
reg have_older_read = 0;
function automatic lite_read_t latest_read;
    latest_read = {expected_r, incoming_control_r, read_fetch_index, read_epoch,
                   pending_control_read, stop_choice, stop_was_required, read_outside_stopped};
endfunction

integer settled_clocks = 0;
reg arm_reply_checked = 0;
integer clocks = 0;

function automatic writes_retired;
    writes_retired = paired_beats == total_samples && w_count == paired_beats
                  && paired_bursts == aw_count && b_count == aw_count;
endfunction

function automatic [31:0] merge_strobes(input [31:0] oldword, input [31:0] newword, input [3:0] enables);
    integer lane;
    begin
        merge_strobes = oldword;
        for (lane=0; lane<4; lane=lane+1)
            if (enables[lane]) merge_strobes[8*lane +: 8] = newword[8*lane +: 8];
    end
endfunction

function automatic descriptor_ok(input [31:0] start_addr, input [31:0] n, input [19:0] h);
    reg [63:0] finish_addr;
    begin
        finish_addr = {32'b0,start_addr} + ({32'b0,n} << SIZE);
        descriptor_ok = (start_addr % W) == 0 && n != 0 && finish_addr <= BYTE_SPAN && h < n;
    end
endfunction

function automatic integer prescribed_beats(input integer offset, input reg is_closed);
    reg [63:0] slot, addr, ring_remaining;
    integer page_remaining;
    begin
        slot = offset % capacity;
        addr = {32'b0,base} + W*slot;
        ring_remaining = capacity-slot;
        page_remaining = (4096-addr[11:0])/W;
        prescribed_beats = BURST_LIMIT;
        if (ring_remaining < prescribed_beats) prescribed_beats = ring_remaining;
        if (page_remaining < prescribed_beats) prescribed_beats = page_remaining;
        if (is_closed && sample_count-offset < prescribed_beats)
            prescribed_beats = sample_count-offset;
    end
endfunction

// Check the entire closed prefix. The burst chosen at each byte address is
// calculated from its physical distance to the two boundaries and the part of
// the stream actually accepted before closure.
task automatic audit_capture;
    integer sample_i, burst_i, beats, p, b;
    reg [63:0] slot, address;
    begin
        if (!writes_retired())
            $fatal(1,"allocation released before write retirement: stream=%0d W=%0d paired=%0d AW=%0d B=%0d",
                     total_samples,w_count,paired_beats,aw_count,b_count);
        sample_i = 0;
        burst_i = first_burst;
        while (sample_i < sample_count) begin
            slot = sample_i % capacity;
            address = {32'b0,base} + W*slot;
            beats = prescribed_beats(sample_i,1);
            if (burst_i >= aw_count || burst_addresses[burst_i] !== address[A-1:0]
                     || burst_lengths[burst_i] != beats)
                $fatal(1,"burst at sample %0d: expected %0d beats at %h, observed %0d at %h",
                         sample_i, beats,address,burst_lengths[burst_i],burst_addresses[burst_i]);
            for (b=0; b<beats; b=b+1) begin
                p = first_sample+sample_i+b;
                if (write_log[p] !== sample_log[p] || written_addresses[p] !== address[A-1:0]+W*b)
                    $fatal(1,"DDR sample %0d at %h differs: expected %h, observed %h at %h",
                           sample_i+b,address+W*b,sample_log[p],write_log[p],written_addresses[p]);
            end
            burst_i = burst_i+1;
            sample_i = sample_i+beats;
        end
        if (burst_i != aw_count) $fatal(1,"unexpected additional write transaction");
    end
endtask

// A stopped read can reserve its beat while the host still holds ARVALID.
// The Lite AR handshake fixes the cursor and lane. Compare the two address
// journals as soon as both entries exist, in either order.
// The address, count, and supply time of each fetched word matter. With
// two stopped reads, DDR issue order can differ from Lite AR order.
integer fetch_binding[$];
reg observed_claimed[$];
task automatic match_fetches;
    integer i, j;
    begin
        if (SCENARIO==30) begin
            while (fetch_binding.size() < expected_fetches.size()) fetch_binding.push_back(-1);
            while (observed_claimed.size() < observed_fetches.size()) observed_claimed.push_back(0);
            for (i=0; i<fetch_binding.size(); i=i+1) begin
                for (j=0; fetch_binding[i]<0 && j<observed_claimed.size(); j=j+1) begin
                    if (!observed_claimed[j] && expected_fetches[i]===observed_fetches[j]) begin
                        fetch_binding[i]=j;
                        observed_claimed[j]=1;
                        matched_fetches=matched_fetches+1;
                    end
                end
            end
        end else while (matched_fetches < expected_fetches.size()
                && matched_fetches < observed_fetches.size()) begin
            if (observed_fetches[matched_fetches] !== expected_fetches[matched_fetches])
                $fatal(1,"DDR read %0d: expected address %h, observed %h",
                       matched_fetches,expected_fetches[matched_fetches],observed_fetches[matched_fetches]);
            matched_fetches = matched_fetches+1;
        end
    end
endtask

// A fetch must belong to a request already accepted or currently offered.
// The latter can reserve its cursor address before the Lite AR handshake.
task automatic check_paired_fetch_offer;
    reg permitted;
    integer i;
    reg [63:0] slot, held_addr;
    begin
        permitted=0;
        for (i=0; i<fetch_binding.size(); i=i+1)
            if (fetch_binding[i]<0 && expected_fetches[i]===araddr) permitted=1;
        if (!permitted && lrv_addr && !lrr_addr && lra==4 && capacity!=0
                && observed_fetches.size()==matched_fetches) begin
            slot=(sample_count-capacity + word_cursor/LANES) % capacity;
            held_addr={32'b0,base}+slot*W;
            permitted=araddr===held_addr[A-1:0];
        end
        if (!permitted) $fatal(1,"DDR fetched address %h without a matching Lite data request",araddr);
    end
endtask

// An accepted read outside STOPPED has unspecified data. If it fetched a
// word, record that observation without advancing either stopped cursor.
task automatic account_optional_fetch;
    begin
        if (!optional_fetch_claimed && observed_fetches.size() > expected_fetches.size()) begin
            expected_fetches.push_back(observed_fetches[expected_fetches.size()]);
            match_fetches();
            optional_fetch_claimed = 1;
        end
    end
endtask

task automatic check_quiet;
    begin
        audit_capture();
        if (matched_fetches != expected_fetches.size() || matched_fetches != ar_count || r_count != ar_count)
            $fatal(1,"data reads did not retire one DDR fetch per Lite request");
        if (natural_close && !sent_interrupt)
            $fatal(1,"completed capture did not pulse the interrupt");
    end
endtask

task automatic observe_lite_ar;
    integer index_in_log, lane;
    reg [63:0] sample_number, addr;
    begin
        if (owe_r) begin
            if (SCENARIO!=30 || have_older_read) $fatal(1,"host exceeded its outstanding read limit");
            oldest_read=latest_read();
            have_older_read=1;
        end
        owe_r = 1;
        read_epoch = epoch;
        read_fetch_index = -1;
        read_outside_stopped = 0;
        optional_fetch_claimed = 0;
        expected_r = 0;
        pending_control_read = 0;
        stop_choice = 0;
        stop_was_required = 0;
        case (lra)
        8: expected_r = staged_b;
        12: expected_r = staged_n;
        0: begin
            expected_r[27] = manual;
            expected_r[26] = disabled;
            expected_r[19:0] = holdoff;
            if (pending_arm) begin
                expected_r[31] = 1;
                // The pending status is specified, as is the first new
                // sampling edge. Firmware may see either complete descriptor's
                // control fields here. Snapshot both words on this AR edge.
                pending_control_read = 1;
                incoming_control_r = {4'b1000, next_control[27:26], 6'b0, next_control[19:0]};
            end else if (triggered) begin
                expected_r[30:28] = 2;
                if (natural_close && writes_retired()) begin
                    stop_choice = 1;
                    stop_was_required = published_stop || settled_clocks > 4096;
                end
            end else if (capacity != 0 && sample_count >= capacity) begin
                expected_r[30:28] = 1;
            end
        end
        // If ARM won arbitration first, only a physically outstanding
        // reservation has a specified result: it must retire before BVALID.
        // The host keeps its original data request stable through any grant.
        4: if (late_reserved_data && !(natural_close && writes_retired() &&
                   (published_stop || sent_interrupt || settled_clocks >= 4096))) begin
            read_outside_stopped = 1;
            account_optional_fetch();
        end
        else begin
            if (!natural_close || !writes_retired() || sample_count < capacity || !capacity)
                $fatal(1,"data port requested before a stopped image existed");
            sample_number = sample_count-capacity + word_cursor/LANES;
            index_in_log = first_sample + sample_number;
            lane = word_cursor % LANES;
            expected_r = sample_log[index_in_log] >> (32*lane);
            addr = {32'b0,base} + W*(sample_number % capacity);
            read_fetch_index = expected_fetches.size();
            expected_fetches.push_back(addr[A-1:0]);
            match_fetches();
            word_cursor = (word_cursor+1) % ({32'b0,capacity}*LANES);
        end
        default: $fatal(1,"unknown Lite read offset");
        endcase
    end
endtask

task automatic execute_command;
    reg [2:0] reply;
    begin
        reply = 0;
        case (command_address)
        8: staged_b = merge_strobes(staged_b,command_data,command_mask);
        12: staged_n = merge_strobes(staged_n,command_data,command_mask);
        4: word_cursor = 0;
        0: begin
            if (command_data[31]) begin
                manual = command_data[27];
                disabled = command_data[26];
            end else if (!descriptor_ok(staged_b,staged_n,command_data[19:0])) begin
                reply[1:0] = 2'b10;
            end else begin
                reply[2] = 1;
                next_base = staged_b;
                next_capacity = staged_n;
                next_control = command_data;
                pending_arm = 1;
                arm_reply_checked = 0;
                active = 0;
                natural_close = 0;
            end
        end
        default: $fatal(1,"unknown Lite write offset");
        endcase
        replies.push_back(reply);
        have_address = 0;
        have_data = 0;
    end
endtask

task automatic accept_sample;
    begin
        sample_log.push_back(td);
        if (triggered) begin
            tail_left = tail_left-1;
        end else if (sample_count >= capacity && (manual || (trig && !disabled))) begin
            triggered = 1;
            tail_left = holdoff;
        end
        if (triggered && tail_left == 0) begin
            active = 0;
            natural_close = 1;
        end
        sample_count = sample_count+1;
        total_samples = total_samples+1;
    end
endtask

task automatic accept_descriptor;
    begin
        base = next_base;
        capacity = next_capacity;
        holdoff = next_control[19:0];
        manual = next_control[27];
        disabled = next_control[26];
        pending_arm = 0;
        active = 1;
        triggered = 0;
        natural_close = 0;
        sent_interrupt = 0;
        published_stop = 0;
        tail_left = 0;
        word_cursor = 0;
        sample_count = 0;
        first_sample = total_samples;
        first_burst = aw_count;
        epoch = epoch + 1;
    end
endtask

reg was_aw_stalled = 0, was_w_stalled = 0, was_ar_stalled = 0;
reg source_stalled = 0;
reg [D-1:0] saved_td;
reg was_lr_stalled = 0, was_lb_stalled = 0;
wire [A+25:0] aw_payload = {awid,awaddr,awlen,awsize,awburst,awlock,awcache,awprot,awqos};
wire [A+25:0] ar_payload = {arid,araddr,arlen,arsize,arburst,arlock,arcache,arprot,arqos};
wire [D+W:0] w_payload = {wdata,wstrb,wlast};
reg [A+25:0] saved_aw,saved_ar;
reg [D+W:0] saved_w;
reg [33:0] saved_lr;
reg [1:0] saved_lb;
task automatic check_stability;
    begin
        if (source_stalled && (tv !== 1'b1 || td !== saved_td))
            $fatal(1,"stream source withdrew or changed a stalled beat");
        source_stalled = tv && !tr;
        saved_td = td;
        if ((awvalid && (^aw_payload) === 1'bx) || (arvalid && (^ar_payload) === 1'bx))
            $fatal(1,"undefined address or control on a valid DDR address channel");
        if (was_aw_stalled && (awvalid !== 1'b1 || aw_payload !== saved_aw)) $fatal(1,"stalled DDR AW changed");
        if (was_w_stalled && (wvalid !== 1'b1 || w_payload !== saved_w)) $fatal(1,"stalled DDR W changed");
        if (was_ar_stalled && (arvalid !== 1'b1 || ar_payload !== saved_ar)) $fatal(1,"stalled DDR AR changed");
        if (was_lr_stalled && (lrv !== 1'b1 || {lrresp,lrd} !== saved_lr)) $fatal(1,"stalled Lite R changed");
        if (was_lb_stalled && (lbv !== 1'b1 || lbresp !== saved_lb)) $fatal(1,"stalled Lite B changed");
        was_aw_stalled = awvalid && !awready;
        was_w_stalled = wvalid && !wready;
        was_ar_stalled = arvalid && !arready;
        was_lr_stalled = lrv && !lrr;
        was_lb_stalled = lbv && !lbr;
        saved_aw = aw_payload;
        saved_w = w_payload;
        saved_ar = ar_payload;
        saved_lr = {lrresp,lrd};
        saved_lb = lbresp;
    end
endtask

// AXI4 reserves allocation bits on a non-modifiable transaction (AxCACHE[1]=0).
function automatic legal_attributes(input id, input [3:0] cache);
    legal_attributes = id === 1'b0 &&
                        (^cache) !== 1'bx && (cache[1] || cache[3:2] == 0);
endfunction

task automatic memory_edge;
    reg [A-1:0] phys;
    reg [D-1:0] value;
    integer prior_write, offset, beats;
    begin
        if (bvalid && bready) b_count = b_count+1;
        if (rvalid && rready) r_count = r_count+1;
        if (awvalid && awready) begin
            if (awburst !== 2'b01 || awsize !== SIZE[2:0] || (^awlen) === 1'bx)
                $fatal(1,"DDR write does not describe a full-width INCR burst");
            if (awlock !== 1'b0 || !legal_attributes(awid,awcache))
                $fatal(1,"invalid DDR write attributes: AWID=%b AWLOCK=%b AWCACHE=%b",awid,awlock,awcache);
            if (!capacity) $fatal(1,"DDR write without a configured allocation");
            offset = addressed_beats-first_sample;
            beats = prescribed_beats(offset,!active);
            if (1+int'(unsigned'(awlen)) != beats)
                $fatal(1,"DDR burst committed at sample %0d while %s: expected %0d beats, observed %0d",
                       offset,active ? "open" : "closed",beats,1+int'(unsigned'(awlen)));
            burst_addresses.push_back(awaddr);
            burst_lengths.push_back(1+int'(unsigned'(awlen)));
            addressed_beats = addressed_beats + 1 + unsigned'(awlen);
            write_ids.push_back(awid);
            aw_count = aw_count+1;
        end
        if (wvalid && wready) begin
            if (wstrb !== {W{1'b1}}) $fatal(1,"DDR write masked a byte lane");
            write_log.push_back(wdata);
            write_lasts.push_back(wlast);
            w_count = w_count+1;
        end
        if (!active && (w_count > total_samples || addressed_beats > total_samples))
            $fatal(1,"DDR write exceeds the entire accepted prefix of a closed capture");
        while (paired_beats < w_count && paired_bursts < aw_count) begin
            phys = burst_addresses[paired_bursts]+W*beat_in_burst;
            written_addresses.push_back(phys);
            beat_in_burst = beat_in_burst+1;
            if (write_lasts[paired_beats] !== (beat_in_burst == burst_lengths[paired_bursts]))
                $fatal(1,"WLAST does not coincide with the transferred AWLEN");
            paired_beats = paired_beats+1;
            if (beat_in_burst == burst_lengths[paired_bursts]) begin
                paired_bursts = paired_bursts+1;
                beat_in_burst = 0;
            end
        end
        if (!bvalid || bready) begin
            bvalid <= allow_b && (!jitter || pattern[15]) && b_count < paired_bursts;
            if (allow_b && (!jitter || pattern[15]) && b_count < paired_bursts)
                bid <= write_ids[b_count];
            else bid <= pattern[8];
        end
        if (arvalid && arready) begin
            if (arburst !== 2'b00 && arburst !== 2'b01)
                $fatal(1,"single-beat DDR read used an invalid burst type");
            if (!legal_attributes(arid,arcache))
                $fatal(1,"invalid DDR read attributes: ARID=%b ARLOCK=%b ARCACHE=%b",arid,arlock,arcache);
            if (arsize !== SIZE[2:0] || arlen !== 0)
                $fatal(1,"DDR read at %h is not one full beat: ARLEN=%0d ARSIZE=%0d",araddr,arlen,arsize);
            if (ar_count >= expected_fetches.size() &&
                    !(ar_count == expected_fetches.size() && (late_reserved_data || (lrv_addr && lra == 4))))
                $fatal(1,"DDR read without a corresponding held or accepted Lite data request");
            if (SCENARIO==30) check_paired_fetch_offer();
            observed_fetches.push_back(araddr);
            if (late_reserved_data && read_outside_stopped) account_optional_fetch();
            match_fetches();
            prior_write = written_addresses.size()-1;
            while (prior_write >= 0 && written_addresses[prior_write] !== araddr)
                prior_write = prior_write-1;
            // Unwritten memory also returns a word. The address and payload
            // journals independently enforce every specified stopped read.
            value = prior_write < 0 ? {LANES{pattern}} : write_log[prior_write];
            fetched_words.push_back(value);
            read_ids.push_back(arid);
            ar_count = ar_count+1;
        end
        if (!rvalid || rready) begin
            rvalid <= allow_r && (!jitter || pattern[21]) && r_count < ar_count;
            if (allow_r && (!jitter || pattern[21]) && r_count < ar_count) begin
                rdata <= fetched_words[r_count];
                rid <= read_ids[r_count];
            end else begin
                rdata <= {LANES{pattern ^ 32'h2e1d539b}};
                rid <= pattern[24];
            end
        end
    end
endtask

reg handoff_now, continuous_source = 0;
reg sampled_reset = 0;
always @(posedge clk) begin
    if (!resetn) begin
        if (sampled_reset && {awvalid,wvalid,arvalid,lbv,lrv} !== 5'b0)
            $fatal(1,"AXI source VALIDs must be low during reset");
        sampled_reset = 1;
    end else begin
        clocks = clocks+1;
        pattern <= {pattern[30:0],pattern[31]^pattern[21]^pattern[1]^pattern[0]};
        check_stability();
        if ((^{awvalid,wvalid,arvalid,lbv,lrv,intr}) === 1'bx)
            $fatal(1,"undefined AXI VALID or interrupt following reset");
        if ((^{lar,lwr,lrr_addr,bready,rready,tr}) === 1'bx)
            $fatal(1,"undefined interface READY following reset");
        handoff_now = 0;
        // These outputs describe state established strictly before this edge.
        if (!active && tr !== 1'b0) $fatal(1,"TREADY asserted outside acquisition");
        if (active && !tv) continuous_source = 0;
        if (ideal_stream && continuous_source && active && tv && tr !== 1'b1) $fatal(1,"gap during ideal continuous acquisition");
        if (intr) begin
            if (!natural_close || pending_arm || sent_interrupt || !writes_retired())
                $fatal(1,"early, suppressed, repeated or unsolicited completion interrupt");
            audit_capture();
            sent_interrupt = 1;
            irq_count = irq_count+1;
        end
        if (natural_close && !jitter && allow_aw && allow_w && allow_b && allow_ar && allow_r)
            settled_clocks = settled_clocks + 1;
        else settled_clocks = 0;
        if (settled_clocks > 4096 && !sent_interrupt)
            $fatal(1,"completed writes did not produce a completion interrupt with ready peers");
        if ((pending_arm && lbv || natural_close && (sent_interrupt || published_stop))
                && (awvalid && awready || wvalid && wready))
            $fatal(1,"DDR write transferred after declaring the closed allocation retired");
        if (lbv) begin
            if (!replies.size()) $fatal(1,"unsolicited Lite write response");
            current_reply = replies[0];
            if (lbresp !== current_reply[1:0])
                $fatal(1,"Lite BRESP %b, expected %b",lbresp,current_reply[1:0]);
            handoff_now = current_reply[2] && lbr;
            if (current_reply[2] && (!arm_reply_checked || lbr)) begin
                audit_capture();
                if (r_count != ar_count || r_count < expected_fetches.size() || arvalid)
                    $fatal(1,"ARM acknowledged before old DDR reads retired");
                arm_reply_checked = 1;
            end
            if (lbr) current_reply = replies.pop_front();
        end
        if (lrv) begin : lite_read_response
            lite_read_t head;
            integer fetch;
            head = have_older_read ? oldest_read : latest_read();
            if (!owe_r || (!head.unspecified && lrresp !== 2'b00))
                $fatal(1,"unsolicited or unsuccessful Lite R");
            // A pass-through return can be visible on both buses in this cycle.
            fetch = head.fetch_index;
            if (SCENARIO==30 && fetch>=0) begin
                fetch=fetch_binding[fetch];
                if (fetch<0) $fatal(1,"Lite data returned without its matching DDR address");
            end
            if (fetch >= 0 && r_count <= fetch && !(rvalid && r_count == fetch))
                $fatal(1,"Lite returned data before DDR supplied the requested beat");
            if (head.may_stop && lrd[30:28] == 3) begin
                if (head.owner == epoch) published_stop = 1;
                head.word_value[30:28] = 3;
            end
            if (head.must_stop && lrd[30:28] !== 3) $fatal(1,"STOPPED regressed before ARM");
            if (!head.unspecified && lrd !== head.word_value) begin
                if (!head.pending_control)
                    $fatal(1,"Lite read snapshot differs: expected %h, got %h",head.word_value,lrd);
                else if (lrd !== head.other_control)
                    $fatal(1,"pending control read: expected %h or %h, got %h",
                           head.word_value,head.other_control,lrd);
            end
            if (lrr) begin
                if (have_older_read) have_older_read=0;
                else owe_r=0;
            end
        end
        if (lrv_addr && lrr_addr) observe_lite_ar();
        if (tv && tr) accept_sample();
        if (lav && lar) begin command_address=la; have_address=1; end
        if (lwv && lwr) begin command_data=lwd; command_mask=lws; have_data=1; end
        if (have_address && have_data) execute_command();
        memory_edge();
        if (handoff_now) begin
            // A pre-handoff DDR reservation for an as-yet-unaccepted Lite
            // request has retired. A subsequent stopped grant uses the new ring.
            if (late_reserved_data && lrv_addr && !owe_r) account_optional_fetch();
            accept_descriptor();
            continuous_source = 1;
        end
    end
end
