#!/bin/bash
set -euo pipefail
cd /app
rm -rf -- out/rtl out/sim
mkdir -p out/rtl out/sim
cp inputs/{axisrle,memscopc,sfifo,skidbuffer}.v out/rtl/

/usr/bin/python3 - <<'PY_INTEGRATE'
from pathlib import Path
import re

src = Path("/app/inputs/memscope.v").read_text()
header, body = src.split("\n\t// Local parameters\n", 1)
body = "\n\t// Local parameters\n" + body
# The first nine parameters stay in their original positional order.
header = header.replace(
    "DEF_HOLDOFF = 0,",
    "DEF_HOLDOFF = 0,\n"
    "\t\t// Enable host-configured DDR capture rings and descriptor handoff.\n"
    "\t\tparameter [0:0]\tOPT_WINDOWED = 1'b0,")
ports = re.findall(r"^\s*(?:input|output)\s+(?:wire|reg)\s+"
                   r"(?:\[[^\n]+?\]\s*)?(\w+),?\s*$", header, re.M)
assert len(ports) == 63
assert ports[0] == "S_AXI_ACLK" and ports[-1] == "o_int"
connections = ",\n".join("\t\t\t." + p + "(" +
                         ("irq" if p == "o_int" else p) + ")" for p in ports)

selection = """
\tgenerate if (OPT_WINDOWED) begin : g_windowed
\t\twire irq;
\t\tmemscope_window #(
\t\t\t.A(C_AXI_ADDR_WIDTH), .D(C_AXI_DATA_WIDTH), .G(LGFIFO),
\t\t\t.ID_WIDTH(C_AXI_ID_WIDTH), .AXI_ID(AXI_ID)
\t\t) windowed (
""" + connections + """
\t\t);
\t\talways @* o_int = irq;
\tend else begin : g_legacy
"""
# Inner conditional generates remain within the single outer generate region.
body = body.replace("generate if (", "if (").replace("end endgenerate", "end")
body = body.replace("endmodule", "\tend endgenerate\nendmodule")
Path("/app/out/rtl/memscope.v").write_text(header + selection + body)
PY_INTEGRATE

cat > out/rtl/memscope_window.v <<'WINDOW_CONTROLLER'
// Windowed mode of the memory-backed WBScope.
`default_nettype none
module memscope_window #(
    parameter A = 32,
    parameter D = 32,
    parameter G = 9,
    parameter ID_WIDTH = 1,
    parameter [ID_WIDTH-1:0] AXI_ID = 0
) (
    input  wire                 S_AXI_ACLK,
    input  wire                 S_AXI_ARESETN,

    input  wire                 i_trigger,
    input  wire                 S_AXIS_TVALID,
    output wire                 S_AXIS_TREADY,
    input  wire [D-1:0]         S_AXIS_TDATA,

    input  wire                 S_AXIL_AWVALID,
    output wire                 S_AXIL_AWREADY,
    input  wire [3:0]           S_AXIL_AWADDR,
    input  wire [2:0]           S_AXIL_AWPROT,
    input  wire                 S_AXIL_WVALID,
    output wire                 S_AXIL_WREADY,
    input  wire [31:0]          S_AXIL_WDATA,
    input  wire [3:0]           S_AXIL_WSTRB,
    output wire                 S_AXIL_BVALID,
    input  wire                 S_AXIL_BREADY,
    output wire [1:0]           S_AXIL_BRESP,
    input  wire                 S_AXIL_ARVALID,
    output wire                 S_AXIL_ARREADY,
    input  wire [3:0]           S_AXIL_ARADDR,
    input  wire [2:0]           S_AXIL_ARPROT,
    output wire                 S_AXIL_RVALID,
    input  wire                 S_AXIL_RREADY,
    output wire [31:0]          S_AXIL_RDATA,
    output wire [1:0]           S_AXIL_RRESP,

    output wire                 M_AXI_AWVALID,
    input  wire                 M_AXI_AWREADY,
    output wire [ID_WIDTH-1:0]  M_AXI_AWID,
    output wire [A-1:0]         M_AXI_AWADDR,
    output wire [7:0]           M_AXI_AWLEN,
    output wire [2:0]           M_AXI_AWSIZE,
    output wire [1:0]           M_AXI_AWBURST,
    output wire                 M_AXI_AWLOCK,
    output wire [3:0]           M_AXI_AWCACHE,
    output wire [2:0]           M_AXI_AWPROT,
    output wire [3:0]           M_AXI_AWQOS,
    output wire                 M_AXI_WVALID,
    input  wire                 M_AXI_WREADY,
    output wire [D-1:0]         M_AXI_WDATA,
    output wire [D/8-1:0]       M_AXI_WSTRB,
    output wire                 M_AXI_WLAST,
    input  wire                 M_AXI_BVALID,
    output wire                 M_AXI_BREADY,
    input  wire [ID_WIDTH-1:0]  M_AXI_BID,
    input  wire [1:0]           M_AXI_BRESP,

    output wire                 M_AXI_ARVALID,
    input  wire                 M_AXI_ARREADY,
    output wire [ID_WIDTH-1:0]  M_AXI_ARID,
    output wire [A-1:0]         M_AXI_ARADDR,
    output wire [7:0]           M_AXI_ARLEN,
    output wire [2:0]           M_AXI_ARSIZE,
    output wire [1:0]           M_AXI_ARBURST,
    output wire                 M_AXI_ARLOCK,
    output wire [3:0]           M_AXI_ARCACHE,
    output wire [2:0]           M_AXI_ARPROT,
    output wire [3:0]           M_AXI_ARQOS,
    input  wire                 M_AXI_RVALID,
    output wire                 M_AXI_RREADY,
    input  wire [ID_WIDTH-1:0]  M_AXI_RID,
    input  wire [D-1:0]         M_AXI_RDATA,
    input  wire                 M_AXI_RLAST,
    input  wire [1:0]           M_AXI_RRESP,
    output reg                  o_int
);

    localparam integer BYTE_SHIFT = $clog2(D/8);
    localparam LANES = D/32;
    localparam LANE_BITS = LANES > 1 ? $clog2(LANES) : 1;
    localparam integer WIDTH_BYTES = D/8;
    localparam [A-1:0] BEAT_BYTES = WIDTH_BYTES[A-1:0];
    localparam integer LAST_LANE = LANES-1;

    wire clk = S_AXI_ACLK;
    wire reset = !S_AXI_ARESETN;

    // Assemble one Lite command from its two independent channel handshakes.
    // The half arriving last bypasses its holding register on that very edge.
    reg address_saved, data_saved;
    reg [1:0] saved_address;
    reg [31:0] saved_data;
    reg [3:0] saved_strb;
    reg arm_pending, arm_drained;

    // Admit four ordinary commands. Once ARM executes, its two staging
    // writes have reserved space even if earlier Lite replies remain held.
    // At most six of eight slots are occupied, so the pointer difference
    // gives the occupancy across wrap.
    localparam [1:0] REPLY_OKAY = 0, REPLY_ARM = 1, REPLY_SLVERR = 2;
    reg [1:0] replies [0:7];
    reg [2:0] reply_wr, reply_rd;
    wire reply_empty = reply_wr == reply_rd;
    wire [2:0] reply_count = reply_wr - reply_rd;
    wire [1:0] oldest_reply = replies[reply_rd];
    wire allow_command = !reset && reply_count < (arm_pending ? 3'd6 : 3'd4);
    assign S_AXIL_AWREADY = allow_command && !address_saved;
    assign S_AXIL_WREADY  = allow_command && !data_saved;
    wire take_aw = S_AXIL_AWVALID && S_AXIL_AWREADY;
    wire take_w  = S_AXIL_WVALID && S_AXIL_WREADY;
    wire execute = (address_saved || take_aw) && (data_saved || take_w);
    wire [1:0] word_address = address_saved ? saved_address : S_AXIL_AWADDR[3:2];
    wire [31:0] command_data = data_saved ? saved_data : S_AXIL_WDATA;
    wire [3:0] command_strb = data_saved ? saved_strb : S_AXIL_WSTRB;

    always @(posedge clk) begin
        if (reset || execute) begin
            address_saved <= 0;
            data_saved <= 0;
        end else begin
            if (take_aw) begin
                saved_address <= S_AXIL_AWADDR[3:2];
                address_saved <= 1;
            end
            if (take_w) begin
                saved_data <= S_AXIL_WDATA;
                saved_strb <= S_AXIL_WSTRB;
                data_saved <= 1;
            end
        end
    end

    reg [31:0] staged_base, staged_capacity;
    integer byte_lane;
    always @(posedge clk) begin
        if (reset) begin
            staged_base <= 0;
            staged_capacity <= 0;
        end else if (execute) begin
            for (byte_lane = 0; byte_lane < 4; byte_lane = byte_lane+1) begin
                if (command_strb[byte_lane] && word_address == 2)
                    staged_base[8*byte_lane +: 8] <= command_data[8*byte_lane +: 8];
                if (command_strb[byte_lane] && word_address == 3)
                    staged_capacity[8*byte_lane +: 8] <= command_data[8*byte_lane +: 8];
            end
        end
    end

    // Keep the carry above a shifted 32-bit count plus a 32-bit byte address.
    // Equality to the exclusive limit is valid, even for the topmost beat.
    localparam MATH_WIDTH = 33 + BYTE_SHIFT;
    localparam [MATH_WIDTH-1:0] ADDRESS_LIMIT = {{(MATH_WIDTH-1){1'b0}}, 1'b1} << A;
    wire [MATH_WIDTH-1:0] end_exclusive =
        {{(MATH_WIDTH-32){1'b0}}, staged_base}
        + {1'b0, staged_capacity, {BYTE_SHIFT{1'b0}}};
    wire valid_descriptor = staged_base[BYTE_SHIFT-1:0] == 0
                         && staged_capacity != 0 && end_exclusive <= ADDRESS_LIMIT
                         && {12'b0, command_data[19:0]} < staged_capacity;
    wire arm_command = execute && word_address == 0 && !command_data[31];
    wire close_for_arm = arm_command && valid_descriptor;
    wire flags_command = execute && word_address == 0 && command_data[31];
    wire rewind = execute && word_address == 1;

    reg [A-1:0] pending_base, pending_last_address;
    reg [31:0] pending_capacity;
    reg [19:0] pending_holdoff;
    reg [1:0] pending_flags;
    wire writer_idle, ddr_read_busy;
    assign S_AXIL_BVALID = !reply_empty && (oldest_reply != REPLY_ARM || arm_drained);
    assign S_AXIL_BRESP  = {(!reply_empty && oldest_reply == REPLY_SLVERR), 1'b0};
    wire take_b = S_AXIL_BVALID && S_AXIL_BREADY;
    wire start = take_b && oldest_reply == REPLY_ARM;
    always @(posedge clk) begin
        if (reset) begin
            arm_pending <= 0;
            arm_drained <= 0;
            pending_base <= 0;
            pending_capacity <= 0;
            pending_last_address <= 0;
            pending_holdoff <= 0;
            pending_flags <= 0;
        end else begin
            if (close_for_arm) begin
                pending_base <= staged_base[A-1:0];
                pending_capacity <= staged_capacity;
                pending_last_address <= end_exclusive[A-1:0] - BEAT_BYTES;
                pending_holdoff <= command_data[19:0];
                pending_flags <= command_data[27:26];
                arm_pending <= 1;
                arm_drained <= 0;
            end else if (start) begin
                arm_pending <= 0;
                arm_drained <= 0;
            end else if (arm_pending && writer_idle && !ddr_read_busy)
                arm_drained <= 1;
        end
    end

    // Tag the ARM when the command executes. It remains in issue order
    // while DDR drains. Concurrent enqueue and dequeue advance both pointers.
    always @(posedge clk) begin
        if (reset) begin
            reply_wr <= 0;
            reply_rd <= 0;
        end else begin
            if (execute) begin
                replies[reply_wr] <= close_for_arm ? REPLY_ARM
                                   : arm_command ? REPLY_SLVERR : REPLY_OKAY;
                reply_wr <= reply_wr + 1'b1;
            end
            if (take_b)
                reply_rd <= reply_rd + 1'b1;
        end
    end

    localparam [2:0] IDLE = 0, FILL = 1, PRIMED = 2, POST_TRIGGER = 3,
                     DRAIN = 4, STOPPED = 5;
    reg [2:0] state;
    reg manual_trigger, disable_trigger;
    reg [19:0] holdoff, samples_after_trigger;
    reg [31:0] samples_to_prime;
    reg [A-1:0] active_base, last_address, next_sample;
    wire acquiring = state == FILL || state == PRIMED || state == POST_TRIGGER;
    wire space;
    assign S_AXIS_TREADY = !reset && acquiring && space;
    wire sample = S_AXIS_TVALID && S_AXIS_TREADY;
    wire trigger_condition = manual_trigger || (i_trigger && !disable_trigger);
    wire finish = state == DRAIN && writer_idle && !close_for_arm;

    always @(posedge clk) begin
        if (reset) begin
            state <= IDLE;
            manual_trigger <= 0;
            disable_trigger <= 0;
            holdoff <= 0;
            samples_after_trigger <= 0;
            samples_to_prime <= 0;
            active_base <= 0;
            last_address <= 0;
            next_sample <= 0;
            o_int <= 0;
        end else begin
            o_int <= finish;
            // next_sample is the circular insertion pointer. After closure it
            // points at the oldest retained sample of a primed ring.
            if (sample)
                next_sample <= (next_sample == last_address)
                             ? active_base : next_sample + BEAT_BYTES;
            if (flags_command)
                {manual_trigger, disable_trigger} <= command_data[27:26];

            if (start) begin
                state <= close_for_arm ? IDLE : FILL;
                active_base <= pending_base;
                last_address <= pending_last_address;
                next_sample <= pending_base;
                samples_to_prime <= pending_capacity;
                holdoff <= pending_holdoff;
                {manual_trigger, disable_trigger} <= pending_flags;
            end else if (close_for_arm) begin
                state <= IDLE;
            end else if (finish) begin
                state <= STOPPED;
            end else if (sample) begin
                case (state)
                FILL: begin
                    samples_to_prime <= samples_to_prime - 1'b1;
                    if (samples_to_prime == 1)
                        state <= PRIMED;
                end
                PRIMED: if (trigger_condition) begin
                    samples_after_trigger <= holdoff;
                    state <= (holdoff == 0) ? DRAIN : POST_TRIGGER;
                end
                POST_TRIGGER: begin
                    samples_after_trigger <= samples_after_trigger - 1'b1;
                    if (samples_after_trigger == 1)
                        state <= DRAIN;
                end
                default: state <= IDLE;
                endcase
            end
        end
    end

    memscope_writer #(.A(A), .D(D), .G(G)) writer (
        .clk(clk), .reset(reset), .start(start), .closed(!acquiring),
        .base(pending_base), .capacity(pending_capacity),
        .push(sample), .data(S_AXIS_TDATA), .space(space), .idle(writer_idle),
        .awvalid(M_AXI_AWVALID), .awready(M_AXI_AWREADY),
        .awaddr(M_AXI_AWADDR), .awlen(M_AXI_AWLEN),
        .wvalid(M_AXI_WVALID), .wready(M_AXI_WREADY),
        .wdata(M_AXI_WDATA), .wlast(M_AXI_WLAST),
        .bvalid(M_AXI_BVALID), .bready(M_AXI_BREADY)
    );
    assign M_AXI_AWID = AXI_ID;
    assign M_AXI_AWSIZE = BYTE_SHIFT[2:0];
    assign M_AXI_AWBURST = 2'b01;
    assign M_AXI_AWLOCK = 0;
    assign M_AXI_AWCACHE = 0;
    assign M_AXI_AWPROT = 0;
    assign M_AXI_AWQOS = 0;
    assign M_AXI_WSTRB = {(D/8){1'b1}};

    // Reserve a response slot on each Lite AR. Register reads complete there
    // immediately. A stopped-data read keeps its own address and lane, then
    // keeps its returned word until it reaches the head and the host takes it.
    wire [2:0] status = (state == PRIMED) ? 3'd1
                      : (state == POST_TRIGGER || state == DRAIN) ? 3'd2
                      : (state == STOPPED) ? 3'd3 : 3'd0;
    wire [31:0] control_word = {arm_pending, status, manual_trigger,
                                disable_trigger, 6'd0, holdoff};
    reg [31:0] register_value;
    always @* begin
        case (S_AXIL_ARADDR[3:2])
        0: register_value = control_word;
        2: register_value = staged_base;
        3: register_value = staged_capacity;
        default: register_value = 0;
        endcase
    end

    reg read_head, read_tail;
    wire [1:0] vacant, needs_address, expects_data, complete;
    wire [31:0] response_word [0:1];
    wire [A-1:0] fetch_address [0:1];

    assign S_AXIL_ARREADY = !reset && vacant[read_tail];
    assign S_AXIL_RVALID = complete[read_head];
    assign S_AXIL_RDATA = response_word[read_head];
    assign S_AXIL_RRESP = 0;
    wire take_read = S_AXIL_ARVALID && S_AXIL_ARREADY;
    wire return_read = S_AXIL_RVALID && S_AXIL_RREADY;
    wire data_port = S_AXIL_ARADDR[3:2] == 1 && state == STOPPED;

    reg [A-1:0] read_cursor;
    reg [LANE_BITS-1:0] lane;
    always @(posedge clk) begin
        if (reset) begin
            read_head <= 0;
            read_tail <= 0;
            lane <= 0;
            read_cursor <= 0;
        end else begin
            if (return_read)
                read_head <= !read_head;
            if (take_read)
                read_tail <= !read_tail;
            if (finish || rewind) begin
                lane <= 0;
                read_cursor <= next_sample;
            end else if (take_read && data_port) begin
                if (lane == LAST_LANE[LANE_BITS-1:0]) begin
                    lane <= 0;
                    read_cursor <= (read_cursor == last_address)
                                 ? active_base : read_cursor + BEAT_BYTES;
                end else
                    lane <= lane + 1'b1;
            end
        end
    end

    // With two slots, the head is older than its complement. Apply that
    // ordering independently to DDR AR, DDR R and Lite R. Thus the second
    // DDR read can retire even while the host stalls the first Lite reply.
    wire ar_slot = needs_address[read_head] ? read_head : !read_head;
    wire r_slot = expects_data[read_head] ? read_head : !read_head;
    assign M_AXI_ARVALID = |needs_address;
    assign M_AXI_ARADDR = fetch_address[ar_slot];
    assign M_AXI_RREADY = |expects_data;
    wire send_address = M_AXI_ARVALID && M_AXI_ARREADY;
    wire receive_data = M_AXI_RVALID && M_AXI_RREADY;

    // A fetch is owed from Lite AR to DDR R, including when Lite AR and
    // ARM execute together. Completed Lite responses do not hold the fence.
    assign ddr_read_busy = (|needs_address) || (|expects_data);

    localparam [1:0] READ_FREE = 0, READ_ADDRESS = 1,
                     READ_MEMORY = 2, READ_COMPLETE = 3;
    genvar slot;
    generate for (slot = 0; slot < 2; slot = slot + 1) begin : read_slots
        reg [1:0] phase;
        reg [31:0] result;
        reg [A-1:0] address;
        reg [LANE_BITS-1:0] saved_lane;
        assign vacant[slot] = phase == READ_FREE;
        assign needs_address[slot] = phase == READ_ADDRESS;
        assign expects_data[slot] = phase == READ_MEMORY;
        assign complete[slot] = phase == READ_COMPLETE;
        assign response_word[slot] = result;
        assign fetch_address[slot] = address;

        always @(posedge clk) begin
            if (reset) begin
                phase <= READ_FREE;
                result <= 0;
                address <= 0;
                saved_lane <= 0;
            end else if (take_read && read_tail == slot) begin
                if (data_port) begin
                    phase <= READ_ADDRESS;
                    address <= read_cursor;
                    saved_lane <= lane;
                end else begin
                    phase <= READ_COMPLETE;
                    result <= register_value;
                end
            end else if (send_address && ar_slot == slot) begin
                phase <= READ_MEMORY;
            end else if (receive_data && r_slot == slot) begin
                phase <= READ_COMPLETE;
                result <= M_AXI_RDATA[32*saved_lane +: 32];
            end else if (return_read && read_head == slot) begin
                phase <= READ_FREE;
            end
        end
    end endgenerate

    assign M_AXI_ARID = AXI_ID;
    assign M_AXI_ARLEN = 0;
    assign M_AXI_ARSIZE = BYTE_SHIFT[2:0];
    assign M_AXI_ARBURST = 2'b01;
    assign M_AXI_ARLOCK = 0;
    assign M_AXI_ARCACHE = 4'b0011;
    assign M_AXI_ARPROT = 0;
    assign M_AXI_ARQOS = 0;

endmodule
WINDOW_CONTROLLER

cat > out/rtl/memscope_writer.v <<'BURST_WRITER'
// Store-and-forward AXI burst writer for the windowed raw capture.
`default_nettype none
module memscope_writer #(
    parameter A = 32,
    parameter D = 32,
    parameter G = 9
) (
    input  wire         clk,
    input  wire         reset,
    input  wire         start,
    input  wire [A-1:0] base,
    input  wire [31:0]  capacity,
    input  wire         closed,
    input  wire         push,
    input  wire [D-1:0] data,
    output wire         space,
    output wire         idle,

    output wire         awvalid,
    input  wire         awready,
    output wire [A-1:0] awaddr,
    output wire [7:0]   awlen,
    output wire         wvalid,
    input  wire         wready,
    output wire [D-1:0] wdata,
    output wire         wlast,
    input  wire         bvalid,
    output wire         bready
);
    localparam BYTE_SHIFT = $clog2(D/8);
    localparam [31:0] MAX_BEATS = 1 << ((G > 8) ? 8 : G-1);

    wire [G:0] fill;
    wire full, empty;
    wire pop = wvalid && wready;
    sfifo #(.BW(D), .LGFLEN(G), .OPT_ASYNC_READ(1'b0)) samples (
        .i_clk(clk), .i_reset(reset || start),
        .i_wr(push), .i_data(data), .o_full(full), .o_fill(fill),
        .i_rd(pop), .o_data(wdata), .o_empty(empty)
    );
    assign space = !full;

    // Keep the ring distance in beats so that an exclusive end address
    // exactly at the top of the address space never needs to be represented.
    reg [A-1:0] ring_base, cursor;
    reg [31:0] ring_size, ring_left;
    wire [31:0] page_beats = (32'd4096 - {20'd0, cursor[11:0]}) >> BYTE_SHIFT;
    wire [31:0] page_limit = (page_beats < MAX_BEATS) ? page_beats : MAX_BEATS;
    wire [31:0] goal = (ring_left < page_limit) ? ring_left : page_limit;
    wire [31:0] buffered_beats = {{(31-G){1'b0}}, fill};
    wire enough_data = buffered_beats >= goal;
    wire [8:0] next_length = enough_data ? goal[8:0] : buffered_beats[8:0];

    // Reserve one B credit as soon as a burst is offered. Two credits let
    // even one-beat bursts transfer every clock with next-cycle responses,
    // while a prolonged B stall places a finite bound on drain time.
    reg [1:0] reserved;
    reg active, address_pending;
    reg [8:0] burst_length, data_pending;
    wire offer = !active && !empty && (closed || enough_data) && reserved < 2;
    wire [8:0] length = active ? burst_length : next_length;
    wire [8:0] remaining = active ? data_pending : next_length;

    wire [31:0] length_beats = {23'd0, length};
    // Include the exclusive end of the address space. A 256-beat burst
    // can advance by a whole 4096-byte page even when A is only 12.
    wire [A:0] cursor_after_burst = {1'b0, cursor}
                                  + {length, {BYTE_SHIFT{1'b0}}};

    // AW and W may make progress in either order. The first cycle of a
    // burst is exposed directly from registered FIFO and cursor state.
    // This allows the next burst to start on the edge following WLAST.
    assign awvalid = active ? address_pending : offer;
    assign awaddr = cursor;
    assign awlen = length[7:0] - 8'd1;
    assign wvalid = active ? (data_pending != 0) : offer;
    assign wlast = wvalid && remaining == 1;
    assign bready = 1'b1;
    wire finish = (active || offer) && (!awvalid || awready)
                                    && (!wvalid || (wready && wlast));

    always @(posedge clk) begin
        if (reset || start) begin
            active <= 1'b0;
            address_pending <= 1'b0;
            burst_length <= 0;
            data_pending <= 0;
            reserved <= 0;
        end else begin
            if (active || offer) begin
                active <= !finish;
                address_pending <= awvalid && !awready;
                burst_length <= length;
                data_pending <= remaining - {8'd0, pop};
            end
            case ({offer, bvalid && bready})
            2'b10: reserved <= reserved + 1'b1;
            2'b01: reserved <= reserved - 1'b1;
            default: reserved <= reserved;
            endcase
        end
    end

    always @(posedge clk) begin
        if (reset) begin
            cursor <= 0;
            ring_base <= 0;
            ring_size <= 0;
            ring_left <= 0;
        end else if (start) begin
            cursor <= base;
            ring_base <= base;
            ring_size <= capacity;
            ring_left <= capacity;
        end else if (finish) begin
            if (ring_left == length_beats) begin
                cursor <= ring_base;
                ring_left <= ring_size;
            end else begin
                cursor <= cursor_after_burst[A-1:0];
                ring_left <= ring_left - length_beats;
            end
        end
    end

    assign idle = empty && !active && reserved == 0;
endmodule
BURST_WRITER

cat > out/sim/tb.sv <<'PIN_LEVEL_BENCH'
// Pin-level regression for windowed raw memscope. Define A_WIDTH, D_WIDTH
// and G_WIDTH for the build under test. A common synchronous reset precedes
// traffic. The scoreboard derives its reference state from interface pins.
`timescale 1ns/1ps
`default_nettype none
`ifndef A_WIDTH
`define A_WIDTH 16
`endif
`ifndef D_WIDTH
`define D_WIDTH 32
`endif
`ifndef G_WIDTH
`define G_WIDTH 5
`endif
module tb;

localparam A = `A_WIDTH, D = `D_WIDTH, G = `G_WIDTH;
localparam BYTES_PER_BEAT = D/8, AXI_SIZE = $clog2(BYTES_PER_BEAT), LANES = D/32;
localparam MAX_BURST = G > 8 ? 256 : 1<<(G-1);
// Only the low word-address bits index the bench RAM. The scoreboard checks
// every full-width bus address independently, including the highest legal beat.
localparam MAX_TRANSFERS = 40000, RAM_WORDS = 65536, ACCESS_BOUND = 4096;
localparam [3:0] CONTROL_REG = 0, DATA_REG = 4, BASE_REG = 8, CAPACITY_REG = 12;
localparam [63:0] LIMIT = 64'h1<<A;
localparam PAGE_BEATS = 4096/BYTES_PER_BEAT;
localparam [31:0] FULL_CAPACITY = LIMIT/BYTES_PER_BEAT;
localparam [19:0] FULL_HOLDOFF = FULL_CAPACITY > 20'hfffff ?
                                 20'hfffff : FULL_CAPACITY-1;
// A=12 has a single physical page. At larger widths the first allocation
// crosses a page after three samples, including at the two-beat burst limit.
localparam [31:0] SPLIT_BASE = A == 12 ? 3*BYTES_PER_BEAT : 4096-3*BYTES_PER_BEAT;
localparam SPLIT_CAPACITY = A == 12 ? PAGE_BEATS-7 :
                            MAX_BURST == PAGE_BEATS ? PAGE_BEATS+1 : MAX_BURST+9;
localparam [31:0] ACTIVE_BASE = A >= 16 ? 32'h2040 : 32'h40;
localparam [31:0] REPLACE_BASE = A == 32 ? 32'h95c34000 :
                                 A >= 16 ? 32'h4000 : 5*BYTES_PER_BEAT;
localparam REPLACE_CAPACITY = PAGE_BEATS == MAX_BURST ? PAGE_BEATS-9 : MAX_BURST+5;

reg clk = 0, rst_n = 0;
always #5 clk = !clk;
reg stream_valid = 0, trigger_input = 0;
reg [D-1:0] stream_data = 0;
wire stream_ready;
reg l_awvalid = 0, l_wvalid = 0, l_bready = 0, l_arvalid = 0, l_rready = 0;
reg [3:0] l_awaddr = 0, l_araddr = 0, l_wstrb = 0;
reg [31:0] l_wdata = 0;
wire l_awready, l_wready, l_bvalid, l_arready, l_rvalid;
wire [1:0] l_bresp, l_rresp;
wire [31:0] l_rdata;
reg m_awready = 0, m_wready = 0, m_bvalid = 0, m_arready = 0, m_rvalid = 0;
reg [D-1:0] m_rdata = 0;
wire m_awvalid, m_wvalid, m_wlast, m_bready, m_arvalid, m_rready, irq;
wire [A-1:0] m_awaddr, m_araddr;
wire [7:0] m_awlen, m_arlen;
wire [2:0] m_awsize, m_arsize;
wire [1:0] m_awburst, m_arburst;
wire [BYTES_PER_BEAT-1:0] m_wstrb;
wire [D-1:0] m_wdata;
wire m_awid, m_arid, m_awlock, m_arlock;

memscope
`ifndef SYNTHESIZED_DUT
 #(.C_AXI_ADDR_WIDTH(A), .C_AXI_DATA_WIDTH(D), .LGFIFO(G),
   .OPT_TREADY_WHILE_IDLE(0), .OPT_WINDOWED(1))
`endif
 dut (
  .S_AXI_ACLK(clk), .S_AXI_ARESETN(rst_n), .S_AXIS_TVALID(stream_valid), .S_AXIS_TREADY(stream_ready), .S_AXIS_TDATA(stream_data), .i_trigger(trigger_input), .o_int(irq),
  .S_AXIL_AWVALID(l_awvalid), .S_AXIL_AWREADY(l_awready), .S_AXIL_AWADDR(l_awaddr), .S_AXIL_AWPROT(3'b0),
  .S_AXIL_WVALID(l_wvalid), .S_AXIL_WREADY(l_wready), .S_AXIL_WDATA(l_wdata), .S_AXIL_WSTRB(l_wstrb),
  .S_AXIL_BVALID(l_bvalid), .S_AXIL_BREADY(l_bready), .S_AXIL_BRESP(l_bresp),
  .S_AXIL_ARVALID(l_arvalid), .S_AXIL_ARREADY(l_arready), .S_AXIL_ARADDR(l_araddr), .S_AXIL_ARPROT(3'b0),
  .S_AXIL_RVALID(l_rvalid), .S_AXIL_RREADY(l_rready), .S_AXIL_RDATA(l_rdata), .S_AXIL_RRESP(l_rresp),
  .M_AXI_AWVALID(m_awvalid), .M_AXI_AWREADY(m_awready), .M_AXI_AWADDR(m_awaddr), .M_AXI_AWLEN(m_awlen), .M_AXI_AWSIZE(m_awsize),
  .M_AXI_AWBURST(m_awburst), .M_AXI_AWID(m_awid), .M_AXI_AWLOCK(m_awlock),
  .M_AXI_WVALID(m_wvalid), .M_AXI_WREADY(m_wready), .M_AXI_WLAST(m_wlast), .M_AXI_WDATA(m_wdata), .M_AXI_WSTRB(m_wstrb),
  .M_AXI_BVALID(m_bvalid), .M_AXI_BREADY(m_bready), .M_AXI_BID(1'b0), .M_AXI_BRESP(2'b0),
  .M_AXI_ARVALID(m_arvalid), .M_AXI_ARREADY(m_arready), .M_AXI_ARADDR(m_araddr), .M_AXI_ARLEN(m_arlen), .M_AXI_ARSIZE(m_arsize),
  .M_AXI_ARBURST(m_arburst), .M_AXI_ARID(m_arid), .M_AXI_ARLOCK(m_arlock),
  .M_AXI_RVALID(m_rvalid), .M_AXI_RREADY(m_rready), .M_AXI_RDATA(m_rdata), .M_AXI_RLAST(1'b1), .M_AXI_RID(1'b0), .M_AXI_RRESP(2'b0)
);

integer cycle = 0, total_samples = 0, phase = 0, good_arms = 0, bad_arms = 0, irqs = 0, max_drain = 0, arm_cycle = 0;
integer expected_bursts = 0, aw_count = 0, w_count = 0, committed_samples = 0, complete_bursts = 0, beat_in_burst = 0, checked_bursts = 0, b_count = 0;
integer expected_reads = 0, ar_count = 0, r_count = 0;
reg [A-1:0] sample_address[0:MAX_TRANSFERS-1], expected_burst_address[0:MAX_TRANSFERS-1], aw_address[0:MAX_TRANSFERS-1], expected_read_address[0:MAX_TRANSFERS-1];
reg [D-1:0] sample_value[0:MAX_TRANSFERS-1], w_payload[0:MAX_TRANSFERS-1], ram[0:RAM_WORDS-1], read_payload[0:MAX_TRANSFERS-1];
reg wlast_value[0:MAX_TRANSFERS-1];
integer expected_length[0:MAX_TRANSFERS-1], aw_length[0:MAX_TRANSFERS-1];
reg [31:0] staged_base = 0, staged_samples = 0, base = 0, capacity = 0, holdoff = 0, pending_base = 0, pending_capacity = 0, pending_holdoff = 0;
reg [1:0] flags = 0, pending_flags = 0;
reg arm_pending = 0, acquiring = 0, capture_closed = 0, triggered = 0, stopped = 0;
integer capture_samples = 0, first_sample = 0, post_remaining = 0, read_word = 0, assembled_beats = 0, target_length = 0;
reg [A-1:0] commit_address;
integer sample_index, byte_lane;
reg have_lite_aw = 0, have_lite_w = 0;
reg [1:0] lite_addr, response_queue[0:2];
reg arm_response[0:2], previous_edge_pending;
integer issued_writes = 0, accepted_responses = 0;
integer overlapped_writes = 0, handoff_collisions = 0;
integer staging_deadline = -1;
reg time_staging = 0, saw_staging_aw = 0, saw_staging_w = 0;
reg [31:0] lite_data;
reg [31:0] expected_rdata[0:1], pending_status_snapshot[0:1];
reg allow_pending_status[0:1];
integer issued_reads = 0, returned_reads = 0, paired_reads = 0, longest_second_ar = 0;
integer read_tail, read_head;
integer return_limit = MAX_TRANSFERS;
integer max_outstanding_writes = 0, longest_staging_acceptance = 0;
reg [3:0] lite_strb;
reg [63:0] extent;
reg ideal = 1, stall_memory_aw = 0, stall_memory_w = 0, stall_memory_b = 0, stall_memory_ar = 0, stall_memory_r = 0;
reg source_enabled = 0, source_gaps = 0;
integer sample_id = 1;
reg [31:0] rng = 32'hfeff0123;
initial if ($value$plusargs("SEED=%h", rng) && rng == 0)
  $fatal(1, "The stall-generator seed must be nonzero");
reg b_transferred = 0, r_transferred = 0;
reg aw_was_stalled = 0, w_was_stalled = 0, ar_was_stalled = 0, lite_r_was_stalled = 0, lite_b_was_stalled = 0;
reg [A+14:0] prior_aw, prior_ar;
reg [D+BYTES_PER_BEAT:0] prior_w;
reg [33:0] prior_lite_r;
reg [1:0] prior_lite_b;
reg stream_was_stalled = 0, continuous_ideal = 0;
reg [D-1:0] prior_stream;
integer checked_rate = 0, single_slot_rate = 0;
integer checked_reads = 0;
reg [A-1:0] observed_read_address[0:MAX_TRANSFERS-1];

task die(input string msg);
begin
  $display("FAIL A%0d D%0d G%0d phase=%0d cycle=%0d sample=%0d capture_samples=%0d aw=%0d/%0d w=%0d resolved=%0d B=%0d AR=%0d R=%0d busy=%0d: %s",
           A, D, G, phase, cycle, total_samples, capture_samples, aw_count, expected_bursts, w_count, committed_samples, b_count, ar_count, r_count, arm_pending, msg);
  $fatal;
end
endtask

function automatic integer burst_goal(input [31:0] addr, input integer slot);
integer goal, pageleft;
begin
  goal = MAX_BURST;
  if (capacity-slot<goal) goal = capacity-slot;
  pageleft = (4096-addr[11:0])/BYTES_PER_BEAT;
  if (pageleft<goal) goal = pageleft;
  burst_goal = goal;
end
endfunction

task finish_expected_burst;
begin
  if (assembled_beats != 0) begin
    expected_burst_address[expected_bursts] = sample_address[total_samples-assembled_beats];
    expected_length[expected_bursts] = assembled_beats;
    expected_bursts = expected_bursts+1;
    assembled_beats = 0;
  end
end
endtask

function automatic [31:0] status;
reg [2:0] state_bits;
begin
  if (arm_pending) state_bits = 0;
  else if (stopped) state_bits = 3;
  else if (triggered) state_bits = 2;
  else if (capacity != 0 && capture_samples >= capacity) state_bits = 1;
  else state_bits = 0;
  status = {arm_pending, state_bits, flags, 6'b0, holdoff[19:0]};
end
endfunction

function automatic [D-1:0] payload(input integer id);
integer lane;
reg [31:0] value;
begin
  for (lane = 0; lane<LANES; lane = lane+1) begin
    value = id * (32'h9e3779b9 + 2*lane) + 32'h85ebca6b*lane;
    value = (value ^ (value >> 16)) * 32'h7feb352d;
    payload[lane*32+:32] = value ^ (value >> 15);
  end
  case (id % 19)
  0, 1, 2: payload = 0;
  3: begin
    payload = 0;
    payload[(id/19)%D] = 1;
  end
  endcase
end
endfunction

always @(posedge clk) begin
  if (!rst_n) begin
    stream_valid <= 0;
    stream_data <= 0;
    sample_id <= 1;
  end else if (!stream_valid || stream_ready) begin
    stream_valid <= source_enabled && (!source_gaps || rng[2:0] != 0);
    if (stream_valid && stream_ready) sample_id <= sample_id+1;
    stream_data <= payload(sample_id + (stream_valid && stream_ready));
  end
end

// This peer permits independent AW and W progress. B is offered on the next
// cycle after the matching address and all W transfers in the ideal mode.
always @(negedge clk) begin
  rng = {rng[30:0], rng[31]^rng[21]^rng[1]^rng[0]};
  m_awready = rst_n && !stall_memory_aw && (ideal || rng[3:2] != 0);
  m_wready = rst_n && !stall_memory_w && (ideal || rng[8:6] >= 3);
  m_arready = rst_n && !stall_memory_ar && (ideal || rng[15:14] != 0);
  if (!rst_n) begin
    m_bvalid = 0;
    m_rvalid = 0;
  end else begin
    if (!m_bvalid || b_transferred)
      m_bvalid = !stall_memory_b && (b_count<complete_bursts) && (ideal || rng[19:18] != 0);
    if (!m_rvalid || r_transferred) begin
      m_rvalid = !stall_memory_r && r_count<return_limit && (r_count<ar_count) && (ideal || rng[25:24] != 0);
      if (m_rvalid) m_rdata = read_payload[r_count];
      else m_rdata = payload(cycle);
    end
  end
end

// Sample pre-edge state. Monitor CPU AR first, then the outgoing stream
// transfer, then control writes. DDR handshakes retire last so the tests catch
// publication of a completion or an ARM response one cycle too soon.
always @(posedge clk) if (rst_n) begin
  cycle = cycle+1;
  if (time_staging) begin
    if (l_awvalid) saw_staging_aw = 1;
    if (l_wvalid) saw_staging_w = 1;
    if (saw_staging_aw && saw_staging_w && staging_deadline < 0)
      staging_deadline = cycle + ACCESS_BOUND - 1;
    if (staging_deadline >= 0 && cycle > staging_deadline)
      die("staging write did not transfer while DDR was stalled");
  end else begin
    saw_staging_aw = 0;
    saw_staging_w = 0;
    staging_deadline = -1;
  end
  if (cycle > 100000) die("global timeout");
  previous_edge_pending = arm_pending;
  b_transferred = m_bvalid && m_bready;
  r_transferred = m_rvalid && m_rready;
  if (stream_was_stalled && (!stream_valid || stream_data !== prior_stream))
    die("bench source changed stalled beat");
  stream_was_stalled = stream_valid && !stream_ready;
  prior_stream = stream_data;
  if (aw_was_stalled && (!m_awvalid || {m_awaddr, m_awlen, m_awsize, m_awburst, m_awid, m_awlock} !== prior_aw))
    die("AW changed under stall");
  if (w_was_stalled && (!m_wvalid || {m_wdata, m_wstrb, m_wlast} !== prior_w))
    die("W changed under stall");
  if (ar_was_stalled && (!m_arvalid || {m_araddr, m_arlen, m_arsize, m_arburst, m_arid, m_arlock} !== prior_ar))
    die("AR changed under stall");
  if (lite_r_was_stalled && (!l_rvalid || {l_rresp, l_rdata} !== prior_lite_r))
    die("Lite R changed under stall");
  if (lite_b_was_stalled && (!l_bvalid || l_bresp !== prior_lite_b))
    die("Lite B changed under stall");
  aw_was_stalled = m_awvalid && !m_awready;
  prior_aw = {m_awaddr, m_awlen, m_awsize, m_awburst, m_awid, m_awlock};
  w_was_stalled = m_wvalid && !m_wready;
  prior_w = {m_wdata, m_wstrb, m_wlast};
  ar_was_stalled = m_arvalid && !m_arready;
  prior_ar = {m_araddr, m_arlen, m_arsize, m_arburst, m_arid, m_arlock};
  lite_r_was_stalled = l_rvalid && !l_rready;
  prior_lite_r = {l_rresp, l_rdata};
  lite_b_was_stalled = l_bvalid && !l_bready;
  prior_lite_b = l_bresp;

  if (irq) begin
    if (arm_pending || !capture_closed || !triggered || stopped)
      die("spurious completion interrupt");
    if (b_count != expected_bursts || committed_samples != total_samples)
      die("interrupt before write retirement");
    stopped = 1;
    irqs = irqs+1;
  end
  if (stream_ready && !acquiring) die("stream READY outside acquisition");
  if (!ideal || stall_memory_aw || stall_memory_w || stall_memory_b || !stream_valid)
    continuous_ideal = 0;
  if (continuous_ideal && acquiring) begin
    if (!stream_ready) die("ideal stream throughput bubble");
    checked_rate = checked_rate+1;
    if (capacity == 1) single_slot_rate = single_slot_rate+1;
  end

  // Compare responses against the outstanding command/read from earlier edges.
  if (l_bvalid) begin
    if (accepted_responses == issued_writes ||
        l_bresp !== response_queue[accepted_responses%3]) die("unexpected Lite B response");
    if (arm_response[accepted_responses%3] && (b_count != expected_bursts || aw_count != expected_bursts
        || committed_samples != total_samples || ar_count != expected_reads || r_count != expected_reads))
      die("ARM BVALID ahead of DDR retirement");
  end
  read_head = returned_reads%2;
  if (l_rvalid) begin
    if (issued_reads == returned_reads || l_rresp !== 0 ||
        (l_rdata !== expected_rdata[read_head] &&
         !(allow_pending_status[read_head] && l_rdata === pending_status_snapshot[read_head]))) begin
      $display("Lite R #%0d expected %x got %x", returned_reads, expected_rdata[read_head], l_rdata);
      die("incorrect Lite R");
    end
  end
  if (l_rvalid && l_rready) returned_reads = returned_reads+1;

  // AR snapshots the pre-edge model, even on command and B handshake edges.
  if (l_arvalid && l_arready) begin
    if (issued_reads-returned_reads >= 2) die("bench exceeded read concurrency");
    if (issued_reads != returned_reads) paired_reads = paired_reads+1;
    read_tail = issued_reads%2;
    issued_reads = issued_reads+1;
    allow_pending_status[read_tail] = arm_pending && l_araddr[3:2] == 0;
    pending_status_snapshot[read_tail] = {1'b1, 3'b0, pending_flags, 6'b0, pending_holdoff[19:0]};
    case (l_araddr[3:2])
    0: expected_rdata[read_tail] = status();
    1: begin
      if (!stopped || arm_pending) die("bench read while not stopped");
      sample_index = first_sample+capture_samples-capacity+read_word/LANES;
      expected_read_address[expected_reads] = sample_address[sample_index];
      expected_reads = expected_reads+1;
      expected_rdata[read_tail] = sample_value[sample_index] >> (32*(read_word%LANES));
      read_word = (read_word+1)%(capacity*LANES);
    end
    2: expected_rdata[read_tail] = staged_base;
    3: expected_rdata[read_tail] = staged_samples;
    endcase
  end

  if (stream_valid && stream_ready) begin
    if (total_samples >= MAX_TRANSFERS) die("sample log capacity exceeded");
    sample_index = capture_samples%capacity;
    sample_address[total_samples] = base+sample_index*BYTES_PER_BEAT;
    sample_value[total_samples] = stream_data;
    if (assembled_beats == 0) target_length = burst_goal(sample_address[total_samples], sample_index);
    assembled_beats = assembled_beats+1;
    total_samples = total_samples+1;
    if (assembled_beats == target_length) finish_expected_burst();
    if (triggered) begin
      post_remaining = post_remaining-1;
      if (post_remaining == 0) begin
        capture_closed = 1;
        acquiring = 0;
        finish_expected_burst();
      end
    end else if (capture_samples >= capacity && (flags[1] || (trigger_input && !flags[0]))) begin
      triggered = 1;
      post_remaining = holdoff;
      if (holdoff == 0) begin
        capture_closed = 1;
        acquiring = 0;
        finish_expected_burst();
      end
    end
    capture_samples = capture_samples+1;
  end

  if (l_bvalid && l_bready) begin
    accepted_responses = accepted_responses+1;
    if (arm_response[(accepted_responses-1)%3]) begin
      if (cycle-arm_cycle>max_drain) max_drain = cycle-arm_cycle;
      arm_pending = 0;
      acquiring = 1;
      continuous_ideal = ideal && !stall_memory_aw && !stall_memory_w && !stall_memory_b;
      capture_closed = 0;
      triggered = 0;
      stopped = 0;
      capacity = pending_capacity;
      base = pending_base;
      holdoff = pending_holdoff;
      flags = pending_flags;
      capture_samples = 0;
      read_word = 0;
      first_sample = total_samples;
    end
  end

  if (l_awvalid && l_awready) begin
    have_lite_aw = 1;
    lite_addr = l_awaddr[3:2];
  end
  if (l_wvalid && l_wready) begin
    have_lite_w = 1;
    lite_data = l_wdata;
    lite_strb = l_wstrb;
  end
  if (have_lite_aw && have_lite_w) begin
    if (issued_writes - accepted_responses >= 3) die("bench exceeded write concurrency");
    if (previous_edge_pending) begin
      if (lite_addr < 2) die("bench overlapped a command that is not a staging write");
      overlapped_writes = overlapped_writes+1;
      if (l_bvalid && l_bready) handoff_collisions = handoff_collisions+1;
    end
    have_lite_aw = 0;
    have_lite_w = 0;
    response_queue[issued_writes%3] = 0;
    arm_response[issued_writes%3] = 0;
    case (lite_addr)
    0: if (lite_data[31]) flags = lite_data[27:26];
       else begin
         extent = {32'b0, staged_base}+({32'b0, staged_samples}<<AXI_SIZE);
         if (staged_samples == 0 || |(staged_base&(BYTES_PER_BEAT-1)) || extent>LIMIT || lite_data[19:0] >= staged_samples) begin
           response_queue[issued_writes%3] = 2;
           bad_arms = bad_arms+1;
         end else begin
           arm_response[issued_writes%3] = 1;
           good_arms = good_arms+1;
           arm_cycle = cycle;
           finish_expected_burst();
           arm_pending = 1;
           acquiring = 0;
           pending_base = staged_base;
           pending_capacity = staged_samples;
           pending_holdoff = lite_data[19:0];
           pending_flags = lite_data[27:26];
         end
       end
    1: read_word = 0;
    2: for (byte_lane = 0; byte_lane<4; byte_lane = byte_lane+1)
         if (lite_strb[byte_lane])
           staged_base[byte_lane*8+:8] = lite_data[byte_lane*8+:8];
    3: for (byte_lane = 0; byte_lane<4; byte_lane = byte_lane+1)
         if (lite_strb[byte_lane])
           staged_samples[byte_lane*8+:8] = lite_data[byte_lane*8+:8];
    endcase
    issued_writes = issued_writes+1;
    if (issued_writes-accepted_responses > max_outstanding_writes)
      max_outstanding_writes = issued_writes-accepted_responses;
    if (time_staging && cycle-staging_deadline+ACCESS_BOUND > longest_staging_acceptance)
      longest_staging_acceptance = cycle-staging_deadline+ACCESS_BOUND;
  end

  // Join address and data by transaction order. W may finish ahead of AW.
  if (m_awvalid && m_awready) begin
    if ({m_awid, m_awlock, m_awburst, m_awsize} !== {1'b0, 1'b0, 2'b01, 3'(AXI_SIZE)})
      die("incorrect AW attributes");
    aw_address[aw_count] = m_awaddr;
    aw_length[aw_count] = m_awlen+1;
    aw_count = aw_count+1;
  end
  if (m_wvalid && m_wready) begin
    if (m_wstrb !== {BYTES_PER_BEAT{1'b1}}) die("partial DDR write");
    if (w_count >= total_samples || m_wdata !== sample_value[w_count]) begin
      $display("W index=%0d expected=%x got=%x", w_count, sample_value[w_count], m_wdata);
      die("incorrect W payload");
    end
    w_payload[w_count] = m_wdata;
    wlast_value[w_count] = m_wlast;
    w_count = w_count+1;
  end
  while (checked_bursts<aw_count && checked_bursts<expected_bursts) begin
    if (aw_address[checked_bursts] !== expected_burst_address[checked_bursts] || aw_length[checked_bursts] !== expected_length[checked_bursts]) begin
      $display("AW index=%0d expected %x +%0d got %x +%0d", checked_bursts, expected_burst_address[checked_bursts], expected_length[checked_bursts], aw_address[checked_bursts], aw_length[checked_bursts]);
      die("incorrect burst geometry");
    end
    checked_bursts = checked_bursts+1;
  end
  while (complete_bursts<aw_count && committed_samples<w_count) begin
    commit_address = aw_address[complete_bursts]+BYTES_PER_BEAT*beat_in_burst;
    if (commit_address !== sample_address[committed_samples])
      die("DDR write beat belongs to wrong allocation or slot");
    if (wlast_value[committed_samples] !== (beat_in_burst+1 == aw_length[complete_bursts]))
      die("incorrect WLAST");
    ram[(commit_address>>AXI_SIZE)&(RAM_WORDS-1)] = w_payload[committed_samples];
    committed_samples = committed_samples+1;
    beat_in_burst = beat_in_burst+1;
    if (beat_in_burst == aw_length[complete_bursts]) begin
      complete_bursts = complete_bursts+1;
      beat_in_burst = 0;
    end
  end
  if (b_transferred) b_count = b_count+1;
  if (m_arvalid && m_arready) begin
    if ({m_arid, m_arburst, m_arsize, m_arlen} !== {1'b0, 2'b01, 3'(AXI_SIZE), 8'b0})
      die("incorrect DDR read attributes");
    observed_read_address[ar_count] = m_araddr;
    read_payload[ar_count] = ram[(m_araddr>>AXI_SIZE)&(RAM_WORDS-1)];
    ar_count = ar_count+1;
  end
  // With two Lite requests outstanding the younger fetch can issue first.
  // Compare an exchanged pair once both journals contain its two entries.
  begin : match_fetches
    while (checked_reads < ar_count && checked_reads < expected_reads) begin
      if (observed_read_address[checked_reads] === expected_read_address[checked_reads]) begin
        checked_reads = checked_reads+1;
      end else begin
        if (checked_reads+1 >= ar_count || checked_reads+1 >= expected_reads)
          disable match_fetches;
        if (observed_read_address[checked_reads] !== expected_read_address[checked_reads+1] ||
            observed_read_address[checked_reads+1] !== expected_read_address[checked_reads])
          die("incorrect DDR read address");
        checked_reads = checked_reads+2;
      end
    end
  end
  if (r_transferred) r_count = r_count+1;
end

task automatic delay_edges(input integer ticks);
  repeat (ticks) @(negedge clk);
endtask

// Start each offer on a falling edge. Keep its payload until the observed
// rising-edge transfer, then poison the inactive channel.
task automatic offer_aw(input [3:0] addr);
begin
  l_awaddr = addr;
  l_awvalid = 1;
  do @(posedge clk); while (!l_awready);
  @(negedge clk);
  l_awvalid = 0;
  l_awaddr = ~addr;
end
endtask

task automatic offer_w(input [31:0] val, input [3:0] strb);
begin
  l_wdata = val;
  l_wstrb = strb;
  l_wvalid = 1;
  do @(posedge clk); while (!l_wready);
  @(negedge clk);
  l_wvalid = 0;
  l_wstrb = ~strb;
  l_wdata = ~val;
end
endtask

task automatic write_request(input [3:0] addr, input [31:0] val, input [3:0] strb,
                             input integer aw_delay = 0, input integer w_delay = 0);
begin
  @(negedge clk);
  fork
    begin delay_edges(aw_delay); offer_aw(addr); end
    begin delay_edges(w_delay);  offer_w(val,strb); end
  join
end
endtask

task automatic take_b(input integer delay_cycles = 0);
begin
  delay_edges(delay_cycles);
  l_bready = 1;
  do @(posedge clk); while (!l_bvalid);
  @(negedge clk);
  l_bready = 0;
end
endtask

task automatic write_cmd(input [3:0] addr, input [31:0] val, input [3:0] strb,
                         input integer aw_delay = 0, input integer w_delay = 0, input integer b_delay = 0);
begin
  write_request(addr,val,strb,aw_delay,w_delay);
  take_b(b_delay);
end
endtask

task automatic read_addr(input [3:0] addr);
integer offer_clocks;
reg second;
begin
  @(negedge clk);
  l_araddr = addr;
  l_arvalid = 1;
  offer_clocks = 0;
  second = issued_reads != returned_reads;
  do begin
    @(posedge clk);
    offer_clocks = offer_clocks+1;
    if (second && !l_arready && offer_clocks == ACCESS_BOUND)
      die("second Lite AR missed the offer deadline");
  end while (!l_arready);
  if (second && offer_clocks > longest_second_ar) longest_second_ar = offer_clocks;
  @(negedge clk);
  l_arvalid = 0;
  l_araddr = ~addr;
end
endtask

task automatic read_rsp(input integer delay_cycles = 0);
begin
  delay_edges(delay_cycles);
  l_rready = 1;
  do @(posedge clk); while (!l_rvalid);
  @(negedge clk);
  l_rready = 0;
end
endtask

task automatic read_reg(input [3:0] addr, input integer delay_cycles = 0);
begin
  read_addr(addr);
  read_rsp(delay_cycles);
end
endtask

task automatic stage(input [31:0] byte_base, input [31:0] sample_capacity);
begin
  write_cmd(BASE_REG, byte_base, 15, 3, 0);
  write_cmd(CAPACITY_REG, sample_capacity, 15, 0, 4);
end
endtask

// Wait on actual samples, never elapsed clocks, when stream or DDR can stall.
task automatic wait_count(input integer num);
begin
  while (capture_samples<num) @(negedge clk);
end
endtask

task automatic wait_stop;
integer start_cycle;
begin
  start_cycle = cycle;
  while (!stopped) begin
    @(negedge clk);
    if (cycle-start_cycle>ACCESS_BOUND) die("natural drain timeout");
  end
  read_reg(CONTROL_REG);
end
endtask

task automatic check_ring;
integer z, idx;
begin
  for (z = 0; z<capacity; z = z+1) begin
    idx = first_sample+capture_samples-capacity+z;
    if (ram[(sample_address[idx]>>AXI_SIZE)&(RAM_WORDS-1)] !== sample_value[idx])
      die("stopped DDR ring contents");
  end
end
endtask

// Withhold the one-slot owner's B responses and reassign it. The incoming
// allocation crosses an interior page where the address space permits it.
// Restage one or both fields using byte enables while keeping ARM at the
// head of the Lite response stream. The scoreboard keys on actual handshakes.
task automatic test_stage_during_arm(input integer trial);
  reg [3:0] first_offset, last_offset, last_strobes, first_strobes;
  reg [31:0] first_value, last_value, base_write, size_write, snapshot_base;
  integer read_index, start_irq, start_bad_arms, expected_replies;
  reg two_writes, reject_next;
begin
  two_writes = trial >= 2;
  reject_next = two_writes && trial[0];
  start_irq = irqs;
  start_bad_arms = bad_arms;
  trigger_input = 0;
  source_enabled = 1;
  stage(32'h200,1);
  write_cmd(CONTROL_REG,0,15);
  snapshot_base = A == 12 ? LIMIT-32*BYTES_PER_BEAT : LIMIT-4096-3*BYTES_PER_BEAT;
  stage(snapshot_base,9);
  stall_memory_b = 1;
  while (complete_bursts == b_count || m_bvalid) @(negedge clk);
  write_request(CONTROL_REG,32'h08000002,15,1,3);
  if (!arm_pending) die("missing ARM ticket");

  // Only the low two base bytes and low capacity byte are enabled. The
  // valid merged pair ends at the address-space limit. In the invalid
  // trials each field alone is legal, but the pair overflows by one beat.
  base_write = ((LIMIT-12*BYTES_PER_BEAT) & 32'hffff) | 32'hba5e0000;
  size_write = 32'hfeed1200 | (reject_next ? 13 : 12);
  last_offset = trial[0] ? CAPACITY_REG : BASE_REG;
  first_offset = trial[0] ? BASE_REG : CAPACITY_REG;
  last_strobes = trial[0] ? 4'b0001 : 4'b0011;
  first_strobes = trial[0] ? 4'b0011 : 4'b0001;
  last_value = trial[0] ? size_write : base_write;
  first_value = trial[0] ? base_write : size_write;

  if (two_writes) begin
    time_staging = 1;
    write_request(first_offset,first_value,first_strobes,trial[0] ? 0 : 4,trial[0] ? 4 : 0);
    time_staging = 0;
    read_reg(first_offset);
  end
  if (trial < 4) begin
    // Offer AR alongside the later write component. AR observes the old
    // or updated staged value according to the actual relative handshakes.
    time_staging = 1;
    fork
      begin
        write_request(last_offset,last_value,last_strobes,trial[0] ? 5 : 0,trial[0] ? 0 : 5);
        time_staging = 0;
      end
      begin delay_edges(5); read_reg(last_offset); end
    join
    read_reg(BASE_REG);
    read_reg(CAPACITY_REG);
    read_reg(CONTROL_REG);
    if (l_bvalid || !arm_pending) die("staging reply bypassed unpaid ARM");
    stall_memory_b = 0;
    wait (l_bvalid);
    read_reg(CONTROL_REG);
    take_b(3);
  end else begin
    // For the second write, offer one component before the DDR release.
    // Put the partner on the ARM B edge or later. Arbitration may separate
    // events and may postpone the first handshake until both are offered.
    @(negedge clk);
    fork
      begin
        if (trial[0]) offer_w(last_value,last_strobes);
        else offer_aw(last_offset);
      end
      begin
        delay_edges(6);
        stall_memory_b = 0;
        wait (l_bvalid);
        @(negedge clk);
        if (trial < 6) begin
          fork
            begin
              if (trial[0]) offer_aw(last_offset);
              else offer_w(last_value,last_strobes);
            end
            take_b();
          join
        end else begin
          take_b();
          delay_edges(4);
          if (trial[0]) offer_aw(last_offset);
          else offer_w(last_value,last_strobes);
        end
      end
    join
  end

  // Keep the younger B responses stalled while sampling under the ARM-time
  // descriptor. The next ARM must see both merged staging words.
  read_reg(BASE_REG);
  read_reg(CAPACITY_REG);
  expected_replies = two_writes ? 2 : 1;
  repeat (expected_replies) take_b(4);
  wait_stop();
  check_ring();
  for (read_index=0; read_index<capacity*LANES+1; read_index=read_index+1)
    read_reg(DATA_REG);
  write_cmd(CONTROL_REG,32'h08000002,15);
  if (reject_next) begin
    if (bad_arms != start_bad_arms+1 || irqs != start_irq+1)
      die("invalid merged descriptor altered completion history");
    read_reg(BASE_REG);
    read_reg(CAPACITY_REG);
    read_reg(CONTROL_REG);
    // Cursor remains advanced by one CPU word, including the lane on D>32.
    for (read_index=0; read_index<LANES+2; read_index=read_index+1) read_reg(DATA_REG);
    check_ring();
  end else begin
    wait_stop();
    check_ring();
    for (read_index=0; read_index<capacity*LANES+2; read_index=read_index+1) read_reg(DATA_REG);
    if (bad_arms != start_bad_arms || irqs != start_irq+2) die("overlap completion count");
  end
end
endtask


// Each accepted stopped-data AR reserves a distinct CPU word. Hold the host's
// RREADY low past ownership transfer, then read back both old results after
// new samples have overwritten the same physical ring.
task automatic test_two_reads_at_handoff(input integer trial);
  integer start_r, words, z, release_cycle;
  reg second_is_control;
begin
  second_is_control = trial >= 2;
  source_enabled = 1;
  trigger_input = 0;
  stage(LIMIT-5*BYTES_PER_BEAT,5);
  write_cmd(CONTROL_REG,32'h08000001,15);
  wait_stop();
  check_ring();

  // The last physical slot supplies the final lane of the first fetch;
  // its successor wraps to the base. This holds for all three beat widths.
  words = LANES*((capacity-1-(capture_samples%capacity)+capacity)%capacity)+LANES-1;
  for (z=0; z<words; z=z+1) read_reg(DATA_REG);
  start_r = r_count;
  stall_memory_ar = !trial[0];
  stall_memory_r = 1;
  read_addr(DATA_REG);
  if (second_is_control) begin
    // The status AR may arbitrate before, on, or after the ARM command edge.
    fork
      read_addr(CONTROL_REG);
      write_request(CONTROL_REG,0,15);
    join
  end else begin
    read_addr(DATA_REG);
    write_request(CONTROL_REG,0,15);
  end

  if (!arm_pending || issued_reads-returned_reads != 2)
    die("setup of the two-read fence");
  delay_edges(11);
  stall_memory_ar = 0;
  if (!second_is_control) begin
    // Let exactly one DDR R retire. A younger fetch may retire first, but
    // neither staggered ordering permits an ARM acknowledgment yet.
    return_limit = start_r+1;
    stall_memory_r = 0;
    while (r_count == start_r) @(negedge clk);
    delay_edges(7);
    if (l_bvalid || !arm_pending) die("ARM ignored one owed DDR response");
  end
  return_limit = MAX_TRANSFERS;
  stall_memory_r = 0;
  release_cycle = cycle;
  take_b();
  if (cycle-release_cycle > ACCESS_BOUND) die("ARM waited too long with host R withheld");
  if (r_count != start_r+(second_is_control ? 1 : 2)) die("unpaid or duplicate fetch");
  wait_count(3*capacity+2);
  read_rsp(3);
  read_rsp(5);
  read_reg(CONTROL_REG);
end
endtask

// Keep the old status in the head slot while the capture goes from primed to
// stopped. A younger data fetch must proceed independently of that held reply.
task automatic test_held_control_read;
  integer start_r, end_cycle;
begin
  read_addr(CONTROL_REG);
  wait (l_rvalid);
  @(negedge clk);
  trigger_input = 1;
  end_cycle = cycle;
  while (!stopped) begin
    @(negedge clk);
    if (cycle-end_cycle > ACCESS_BOUND) die("completion with an old status response held");
  end
  trigger_input = 0;
  start_r = r_count;
  read_addr(DATA_REG);
  // RREADY may now resume. Compare both words to their independent snapshots.
  read_rsp(11);
  read_rsp();
  if (r_count != start_r+1) die("stopped data-read fetch count");
  read_reg(CONTROL_REG);
  check_ring();
end
endtask

reg [31:0] top_base;
integer j, initial_irqs;
initial begin
  repeat (4) @(negedge clk);
  rst_n = 1;
  repeat (4) @(negedge clk);

  phase = 1;
  read_reg(CONTROL_REG);
  read_reg(BASE_REG);
  read_reg(CAPACITY_REG);
  write_cmd(CONTROL_REG, 0, 15);
  write_cmd(BASE_REG, 32'hfacecafe, 4'b1101);
  write_cmd(BASE_REG, 32'h10203040, 4'b0110, 2, 0);
  read_reg(BASE_REG);
  write_cmd(CAPACITY_REG, 32'hfacecafe, 4'b1011);
  write_cmd(CAPACITY_REG, 32'h01020304, 4'b0101, 0, 3);
  read_reg(CAPACITY_REG);
  stage(1, 8);
  write_cmd(CONTROL_REG, 0, 15); // misalignment
  stage(0, 32'hffffffff);
  write_cmd(CONTROL_REG, 0, 15); // wide shift overflows address domain
  top_base = LIMIT-BYTES_PER_BEAT*9;
  stage(top_base, 10);
  write_cmd(CONTROL_REG, 0, 15);
  stage(top_base, 9);
  write_cmd(CONTROL_REG, 9, 15);
  stage(0, 7);
  write_cmd(CONTROL_REG, 20'h80000, 15); // high unsigned holdoff, H>N
  for (j = 0; j < AXI_SIZE; j = j+1) begin
    stage(1<<j, 3);
    write_cmd(CONTROL_REG, 1, 15);
  end
  stage(0, FULL_CAPACITY+1);
  write_cmd(CONTROL_REG, 0, 15);
  if (A < 32) begin
    stage(LIMIT, 1);
    write_cmd(CONTROL_REG, 0, 15);
  end

  phase = 2;
  // Continuous ADC, several full bursts, a 4K split, a short ring tail and wraps.
  source_enabled = 1;
  trigger_input = 1;
  stage(SPLIT_BASE, SPLIT_CAPACITY);
  write_cmd(CONTROL_REG, 32'h04000007, 15, 0, 4, 4);
  wait_count(SPLIT_CAPACITY*2+21);
  read_reg(CONTROL_REG, 5);
  stage(1, 12);
  write_cmd(CONTROL_REG, 0, 15, 4, 0);
  read_reg(CONTROL_REG); // invalid ARM from running
  stage(top_base, 9);
  wait_count(SPLIT_CAPACITY*3+43);

  phase = 3;
  // Suppress and drain the old prefix; old data/addresses and flags remain
  // observable until the B ownership boundary. Stream stays continuously valid.
  stall_memory_b = 1;
  fork
    write_cmd(CONTROL_REG, 32'h08000002, 15, 0, 5, 18);
    begin
      wait (arm_pending);
      read_reg(CONTROL_REG, 3);
      delay_edges(17);
      stall_memory_b = 0;
      read_reg(CONTROL_REG);
    end
  join
  wait_stop();
  check_ring();
  for (j = 0; j<9*LANES+5; j = j+1) read_reg(DATA_REG, j%3);
  write_cmd(DATA_REG, 0, 15, 2, 0);
  for (j = 0; j<2*LANES+1; j = j+1) read_reg(DATA_REG);

  phase = 4;
  // Accept a stopped read and keep its Lite response stalled across ARM.
  // Delay both DDR address and return, then finish reading the old snapshot
  // while samples already transfer to the new allocation.
  stage(ACTIVE_BASE, 17);
  stall_memory_ar = 1;
  stall_memory_r = 1;
  trigger_input = 0;
  fork
    begin
      read_addr(DATA_REG);
      write_cmd(CONTROL_REG, 32'h00000006, 15);
    end
    begin
      delay_edges(25);
      stall_memory_ar = 0;
      delay_edges(19);
      stall_memory_r = 0;
    end
  join
  wait_count(10);
  read_rsp(6);

  phase = 5;
  ideal = 0;
  source_gaps = 1;
  wait_count(33);
  // First disable external trigger, then exercise MANUAL and DISABLE on the
  // execution edge. Old flags govern any stream beat on that edge.
  write_cmd(CONTROL_REG, 32'h840fffff, 15);
  trigger_input = 1;
  wait_count(55);
  read_reg(CONTROL_REG);
  write_cmd(CONTROL_REG, 32'h8c0fffff, 15, 4, 0);
  wait (triggered);
  wait (capture_closed);
  stall_memory_b = 1;
  @(negedge clk);
  read_reg(CONTROL_REG, 4);
  initial_irqs = irqs;
  stage(REPLACE_BASE, REPLACE_CAPACITY);
  fork
    write_cmd(CONTROL_REG, 32'h00000003, 15, 5, 0);
    begin
      wait (arm_pending);
      read_reg(CONTROL_REG);
      delay_edges(18);
      stall_memory_b = 0;
    end
  join
  if (irqs != initial_irqs) die("cancelled interrupt");
  trigger_input = 0;
  // Randomized AXI and sparse source, independent prolonged W / AW stalls,
  // multiple outstanding bursts and non-power-of-two wrap.

  phase = 6;
  stall_memory_w = 1;
  delay_edges(40);
  stall_memory_w = 0;
  wait_count(MAX_BURST+9);
  stall_memory_aw = 1;
  delay_edges(57);
  stall_memory_aw = 0;
  wait_count(MAX_BURST*8+47);
  trigger_input = 1;
  wait (triggered);
  read_reg(CONTROL_REG);
  wait_stop();
  check_ring();
  for (j = 0; j<LANES*7+1; j = j+1) read_reg(DATA_REG, j%5);
  write_cmd(DATA_REG, 32'hdeadbeef, 15);
  read_reg(DATA_REG);

  phase = 7;
  // One-beat bursts at one sample per clock force adjacent AW and WLAST
  // handshakes and one B per clock. Use H=0, then idle the source on closure.
  source_gaps = 0;
  ideal = 1;
  trigger_input = 0;
  stage(16'h0100, 1);
  write_cmd(CONTROL_REG, 0, 15);
  wait_count(2048);
  read_reg(CONTROL_REG);
  @(negedge clk);
  trigger_input = 1;
  source_enabled = 0;
  wait_stop();
  check_ring();
  for (j = 0; j<LANES*3; j = j+1) read_reg(DATA_REG);

  phase = 8;
  // Test the largest legal capacity and holdoff at each width. For the
  // one-page address space, pass through the top address, wrap and complete
  // the entire H=N-1 tail while the source stays continuously valid.
  stage(0, FULL_CAPACITY);
  source_enabled = 1;
  trigger_input = 0;
  write_cmd(CONTROL_REG, FULL_HOLDOFF, 15);
  read_reg(CONTROL_REG);
  if (A == 12) begin
    wait_count(2*FULL_CAPACITY+3);
    trigger_input = 1;
    wait_count(3*FULL_CAPACITY+3);
    wait_stop();
    check_ring();
    for (j = 0; j<2*LANES+3; j = j+1) read_reg(DATA_REG);
    trigger_input = 0;
    write_cmd(CONTROL_REG, FULL_HOLDOFF, 15);
  end
  wait_count(31);
  // Force W and AW to freeze after a short prefix, release to an idle producer.
  stall_memory_aw = 1;
  stall_memory_w = 1;
  source_enabled = 0;
  delay_edges(3);
  stage(16'h700, 5);
  fork
    write_cmd(CONTROL_REG, 32'h08000000, 15, 0, 2);
    begin
      wait (arm_pending);
      delay_edges(33);
      stall_memory_w = 0;
      delay_edges(9);
      stall_memory_aw = 0;
    end
  join
  source_enabled = 1;
  wait_stop();
  check_ring();
  for (j = 0; j<5*LANES+2; j = j+1) read_reg(DATA_REG);
  repeat (7) @(negedge clk);
  if (expected_bursts != checked_bursts || w_count != total_samples ||
      expected_reads != ar_count || expected_reads != r_count || b_count != expected_bursts ||
      single_slot_rate < 2048 || irqs != (A == 12 ? 5 : 4))
    die("incomplete transfers");
  phase = 9;
  for (j=0; j<8; j=j+1) test_stage_during_arm(j);
  if (max_outstanding_writes != 3 || issued_writes != accepted_responses || checked_bursts != expected_bursts ||
      w_count != total_samples || expected_reads != r_count || b_count != expected_bursts)
    die("unfinished work after pipelined writes");
  phase = 10;
  for (j=0; j<4; j=j+1) test_two_reads_at_handoff(j);
  test_held_control_read();
  if (paired_reads < 5 || issued_reads != returned_reads || checked_reads != expected_reads ||
      expected_reads != ar_count || expected_reads != r_count) die("unfinished two-read trial");
  $display("READS issued=%0d returned=%0d second_ARs=%0d max_second_AR_offer_clocks=%0d",
           issued_reads, returned_reads, paired_reads, longest_second_ar);
  $display("OVERLAP staged_while_pending=%0d staged_on_ARM_B=%0d max_outstanding=%0d max_staging_accept_clocks=%0d writes=%0d responses=%0d",
           overlapped_writes, handoff_collisions, max_outstanding_writes, longest_staging_acceptance,
           issued_writes, accepted_responses);
  $display("PASS A=%0d D=%0d G=%0d cycles=%0d samples=%0d bursts=%0d DDR_reads=%0d valid_ARMs=%0d invalid_ARMs=%0d natural_IRQ=%0d ideal_beats=%0d one_slot_ideal_beats=%0d max_observed_ARM_clocks=%0d",
           A, D, G, cycle, total_samples, expected_bursts, expected_reads, good_arms, bad_arms, irqs, checked_rate, single_slot_rate, max_drain);
  $finish;
end
endmodule
PIN_LEVEL_BENCH

# Exercise the generated hierarchy through its public pins.
workdir=$(mktemp -d)
trap 'rm -rf -- "$workdir"' EXIT
cd /app/out/rtl
/usr/bin/iverilog -g2012 -s tb -o "$workdir/capture.vvp" *.v ../sim/tb.sv
/usr/bin/vvp "$workdir/capture.vvp"
