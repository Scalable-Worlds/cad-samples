`timescale 1ns/1ps
`default_nettype none
`ifdef COMPRESSED_CLIENT
`define SUBJECT memscopc
`define BASELINE baseline_memscopc
`else
`define SUBJECT memscope
`define BASELINE baseline_memscope
`endif

// One copy drives both libraries with the same legal bus handshakes.
module compare_legacy #(
    parameter DEFAULTS=0, A=32, D=32, I=1, IDLE=1, G=9, L=16,
              ID=0, H=20, DEFAULT_H=0, CYCLES=2300
) (input wire clk, output reg done=0);
    // A wider output binding also observes the public port's sign extension.
    typedef struct packed {
        logic [1:0] tready;
        logic [1:0] l_awready;
        logic [1:0] l_wready;
        logic [1:0] l_bvalid;
        logic [2:0] l_bresp;
        logic [1:0] l_arready;
        logic [1:0] l_rvalid;
        logic [32:0] l_rdata;
        logic [2:0] l_rresp;
        logic [1:0] m_awvalid;
        logic [I:0] m_awid;
        logic [A:0] m_awaddr;
        logic [8:0] m_awlen;
        logic [3:0] m_awsize;
        logic [2:0] m_awburst;
        logic [1:0] m_awlock;
        logic [4:0] m_awcache;
        logic [3:0] m_awprot;
        logic [4:0] m_awqos;
        logic [1:0] m_wvalid;
        logic [D:0] m_wdata;
        logic [D/8:0] m_wstrb;
        logic [1:0] m_wlast;
        logic [1:0] m_bready;
        logic [1:0] m_arvalid;
        logic [I:0] m_arid;
        logic [A:0] m_araddr;
        logic [8:0] m_arlen;
        logic [3:0] m_arsize;
        logic [2:0] m_arburst;
        logic [1:0] m_arlock;
        logic [4:0] m_arcache;
        logic [3:0] m_arprot;
        logic [4:0] m_arqos;
        logic [1:0] m_rready;
        logic [1:0] irq;
    } pins_t;
    wire pins_t original, revised;
    reg resetn=0;
    reg trigger=0;
    reg tvalid=0;
    reg [D-1:0] tdata=0;
    reg l_awvalid=0;
    reg [4-1:0] l_awaddr=0;
    reg [2:0] l_awprot=0;
    reg l_wvalid=0;
    reg [32-1:0] l_wdata=0;
    reg [32/8-1:0] l_wstrb=0;
    reg l_bready=0;
    reg l_arvalid=0;
    reg [4-1:0] l_araddr=0;
    reg [2:0] l_arprot=0;
    reg l_rready=0;
    reg m_awready=0;
    reg m_wready=0;
    reg m_bvalid=0;
    reg [I-1:0] m_bid=0;
    reg [1:0] m_bresp=0;
    reg m_arready=0;
    reg m_rvalid=0;
    reg [I-1:0] m_rid=0;
    reg [D-1:0] m_rdata=0;
    reg m_rlast=0;
    reg [1:0] m_rresp=0;
`ifdef COMPRESSED_CLIENT
    wire [D-2:0] sample = tdata[D-2:0];
`else
    wire [D-1:0] sample = tdata;
`endif
`define NAMED_PINS(p) \
    .S_AXI_ACLK(clk), .S_AXI_ARESETN(resetn), \
    .i_trigger(trigger), .S_AXIS_TVALID(tvalid), \
    .S_AXIS_TREADY(p.tready), .S_AXIS_TDATA(sample), \
    .S_AXIL_AWVALID(l_awvalid), .S_AXIL_AWREADY(p.l_awready), \
    .S_AXIL_AWADDR(l_awaddr), .S_AXIL_AWPROT(l_awprot), \
    .S_AXIL_WVALID(l_wvalid), .S_AXIL_WREADY(p.l_wready), \
    .S_AXIL_WDATA(l_wdata), .S_AXIL_WSTRB(l_wstrb), \
    .S_AXIL_BVALID(p.l_bvalid), .S_AXIL_BREADY(l_bready), \
    .S_AXIL_BRESP(p.l_bresp), .S_AXIL_ARVALID(l_arvalid), \
    .S_AXIL_ARREADY(p.l_arready), .S_AXIL_ARADDR(l_araddr), \
    .S_AXIL_ARPROT(l_arprot), .S_AXIL_RVALID(p.l_rvalid), \
    .S_AXIL_RREADY(l_rready), .S_AXIL_RDATA(p.l_rdata), \
    .S_AXIL_RRESP(p.l_rresp), .M_AXI_AWVALID(p.m_awvalid), \
    .M_AXI_AWREADY(m_awready), .M_AXI_AWID(p.m_awid), \
    .M_AXI_AWADDR(p.m_awaddr), .M_AXI_AWLEN(p.m_awlen), \
    .M_AXI_AWSIZE(p.m_awsize), .M_AXI_AWBURST(p.m_awburst), \
    .M_AXI_AWLOCK(p.m_awlock), .M_AXI_AWCACHE(p.m_awcache), \
    .M_AXI_AWPROT(p.m_awprot), .M_AXI_AWQOS(p.m_awqos), \
    .M_AXI_WVALID(p.m_wvalid), .M_AXI_WREADY(m_wready), \
    .M_AXI_WDATA(p.m_wdata), .M_AXI_WSTRB(p.m_wstrb), \
    .M_AXI_WLAST(p.m_wlast), .M_AXI_BVALID(m_bvalid), \
    .M_AXI_BREADY(p.m_bready), .M_AXI_BID(m_bid), \
    .M_AXI_BRESP(m_bresp), .M_AXI_ARVALID(p.m_arvalid), \
    .M_AXI_ARREADY(m_arready), .M_AXI_ARID(p.m_arid), \
    .M_AXI_ARADDR(p.m_araddr), .M_AXI_ARLEN(p.m_arlen), \
    .M_AXI_ARSIZE(p.m_arsize), .M_AXI_ARBURST(p.m_arburst), \
    .M_AXI_ARLOCK(p.m_arlock), .M_AXI_ARCACHE(p.m_arcache), \
    .M_AXI_ARPROT(p.m_arprot), .M_AXI_ARQOS(p.m_arqos), \
    .M_AXI_RVALID(m_rvalid), .M_AXI_RREADY(p.m_rready), \
    .M_AXI_RID(m_rid), .M_AXI_RDATA(m_rdata), \
    .M_AXI_RLAST(m_rlast), .M_AXI_RRESP(m_rresp), \
    .o_int(p.irq)
`define ORDERED_PINS(p) \
    clk, resetn, trigger, tvalid, p.tready, sample, \
    l_awvalid, p.l_awready, l_awaddr, l_awprot, l_wvalid, p.l_wready, \
    l_wdata, l_wstrb, p.l_bvalid, l_bready, p.l_bresp, l_arvalid, \
    p.l_arready, l_araddr, l_arprot, p.l_rvalid, l_rready, p.l_rdata, \
    p.l_rresp, p.m_awvalid, m_awready, p.m_awid, p.m_awaddr, p.m_awlen, \
    p.m_awsize, p.m_awburst, p.m_awlock, p.m_awcache, p.m_awprot, p.m_awqos, \
    p.m_wvalid, m_wready, p.m_wdata, p.m_wstrb, p.m_wlast, m_bvalid, \
    p.m_bready, m_bid, m_bresp, p.m_arvalid, m_arready, p.m_arid, \
    p.m_araddr, p.m_arlen, p.m_arsize, p.m_arburst, p.m_arlock, p.m_arcache, \
    p.m_arprot, p.m_arqos, m_rvalid, p.m_rready, m_rid, m_rdata, \
    m_rlast, m_rresp, p.irq
    generate if (DEFAULTS) begin : binding
        `SUBJECT dut(`NAMED_PINS(revised));
        `BASELINE ref_dut(`NAMED_PINS(original));
`ifndef FIXED_LEGACY_NETLIST
`ifndef SUPPLIED_RAW_CLIENT
`ifndef COMPRESSED_CLIENT
        initial if (dut.OPT_WINDOWED !== 0) $fatal(1,"OPT_WINDOWED must default to zero");
`endif
`endif
`endif
    end else begin : binding
`ifdef FIXED_LEGACY_NETLIST
        `SUBJECT dut(`ORDERED_PINS(revised));
`else
        `SUBJECT #(A,D,I,IDLE,G,L,ID,H,DEFAULT_H) dut(`ORDERED_PINS(revised));
`endif
        `BASELINE #(A,D,I,IDLE,G,L,ID,H,DEFAULT_H) ref_dut(`ORDERED_PINS(original));
    end endgenerate
`undef NAMED_PINS
`undef ORDERED_PINS

    function automatic [31:0] advance(input [31:0] v);
        begin
            v = v ^ (v << 13);
            v = v ^ (v >> 17);
            advance = v ^ (v << 5);
        end
    endfunction

    task compare_pins(input integer cycle_no);
        integer bit_no;
        pins_t expected, actual;
        begin
            expected = original;
            actual = revised;
`ifdef FIXED_LEGACY_NETLIST
            // A gate network may propagate unspecified startup RAM bits
            // into idle payloads. Compare the protocol values while offered.
            if (!original.l_rvalid) begin expected.l_rdata=0; actual.l_rdata=0; end
            if (!original.m_awvalid) begin
                expected.m_awaddr=0; actual.m_awaddr=0;
                expected.m_awlen=0; actual.m_awlen=0;
            end
            if (!original.m_wvalid) begin expected.m_wdata=0; actual.m_wdata=0; end
            if (!original.m_arvalid) begin expected.m_araddr=0; actual.m_araddr=0; end
`endif
            if (actual !== expected) begin
                // A bit the original leaves uninitialized can legally be defined.
                for (bit_no=0; bit_no < $bits(pins_t); bit_no=bit_no+1) begin
                    if ((expected[bit_no] === 1'b0 || expected[bit_no] === 1'b1)
                            && actual[bit_no] !== expected[bit_no])
                        $fatal(1, "Legacy pin difference at cycle %0d (A=%0d D=%0d G=%0d), output bit %0d: reference %h observed %h",
                              cycle_no, A, D, G, bit_no, expected, actual);
                end
            end
        end
    endtask

    integer n, lane;
    reg [31:0] random_word=32'h93f9689d ^ A ^ (D<<16);
    reg was_sample=0, was_b=0, was_r=0;
    reg write_active=0, need_aw=0, need_w=0, read_active=0;
    reg want_read=0;
    integer writes=0, aw_total=0, last_total=0, b_total=0;
    integer ar_total=0, r_total=0;
    reg [A-1:0] read_address;

    initial begin
        @(posedge clk);
        // The first sampled reset has taken effect on each subsequent edge.
        for (n=-3; n<0; n=n+1) begin
            @(posedge clk);
            compare_pins(n);
        end
        for (n=0; n<CYCLES; n=n+1) begin
            @(negedge clk);
            resetn=1;
            random_word=advance(random_word);
            if (!tvalid || was_sample) begin
                tvalid=random_word[0] || random_word[1];
                // Alternate stretches of identical data with incompressible words.
                if ((n % 257 < 90) || (random_word[8:6]==0))
                    for (lane=0; lane<D/32; lane=lane+1)
                        tdata[32*lane +: 32]=advance(random_word ^ (lane*32'h9c58ad71));
            end
            trigger = random_word[18:13] == 0;
            l_bready=random_word[3] || random_word[4];
            l_rready=random_word[5] || random_word[6];
            m_awready=random_word[9] || random_word[10];
            m_wready=random_word[11] || random_word[12];
            m_arready=random_word[17] || random_word[18];
            if (!m_bvalid || was_b) begin
                m_bvalid=(b_total<aw_total && b_total<last_total && random_word[19]);
                m_bid=ID;
            end
            if (!m_rvalid || was_r) begin
                m_rvalid=(r_total<ar_total && random_word[20]);
                m_rid=ID;
                m_rlast=1;
                for (lane=0; lane<D/32; lane=lane+1)
                    m_rdata[32*lane +: 32]=advance(32'hf0a98741 ^ read_address ^ lane);
            end
            // A and W are launched independently. Both remain asserted until accepted.
            if (!write_active && (n<5 || random_word[26:21]==0)) begin
                write_active=1;
                need_aw=1;
                need_w=1;
                l_awaddr=(writes%6<3) ? 0 : {random_word[28:27],2'b0};
                l_awprot=random_word[10:8];
                l_wdata={1'b1,3'b0,random_word[28:27],6'b0,random_word[19:0]};
                if (writes==0 || writes%31==0) l_wdata=32'd9;
                l_wstrb=(writes%3==0) ? 4'hf : random_word[3:0];
                writes=writes+1;
            end
            l_awvalid=need_aw && (l_awvalid || random_word[23]);
            l_wvalid=need_w && (l_wvalid || random_word[24]);
            if (!read_active && random_word[6:4]==0) begin
                read_active=1;
                want_read=1;
                l_araddr={random_word[18:17],2'b0};
                l_arprot=random_word[15:13];
            end
            l_arvalid=want_read;
            @(posedge clk);
            compare_pins(n);
            was_sample=tvalid && original.tready;
            was_b=m_bvalid && original.m_bready;
            was_r=m_rvalid && original.m_rready;
            if (l_awvalid && original.l_awready) need_aw=0;
            if (l_wvalid && original.l_wready) need_w=0;
            if (l_bready && original.l_bvalid) write_active=0;
            if (l_arvalid && original.l_arready) want_read=0;
            if (l_rready && original.l_rvalid) read_active=0;
            if (original.m_awvalid && m_awready) aw_total=aw_total+1;
            if (original.m_wvalid && m_wready && original.m_wlast) last_total=last_total+1;
            if (was_b) b_total=b_total+1;
            if (original.m_arvalid && m_arready) begin
                ar_total=ar_total+1;
                read_address=original.m_araddr;
            end
            if (was_r) r_total=r_total+1;
        end
        // A later common reset cancels all host and memory transactions.
        @(negedge clk);
        resetn=0;
        tvalid=0;
        l_awvalid=0; l_wvalid=0; l_arvalid=0;
        m_bvalid=0; m_rvalid=0;
        repeat (4) begin
            @(posedge clk);
            compare_pins(n);
            n=n+1;
        end
        @(negedge clk);
        resetn=1;
        repeat (4) begin
            @(posedge clk);
            compare_pins(n);
            n=n+1;
        end
        done=1;
    end
endmodule

module tb;
    reg clk=0;
    always #5 clk=!clk;
`ifdef FIXED_LEGACY_NETLIST
    wire done;
`ifdef POSITIONAL_LEGACY_NETLIST
    compare_legacy #(.A(16),.D(64),.I(3),.IDLE(0),.G(5),.L(15),.ID(5),
                     .H(12),.DEFAULT_H(37)) fixed(clk,done);
`else
    compare_legacy #(.DEFAULTS(1)) fixed(clk,done);
`endif
`elsif LONG_LEGACY_CAPTURE
    wire done;
    compare_legacy #(.A(16),.D(64),.I(3),.IDLE(0),.G(5),.L(15),.ID(5),
                     .H(12),.DEFAULT_H(37),.CYCLES((1<<16)+512)) spanning(clk,done);
`else
    wire [3:0] done;
    compare_legacy #(.DEFAULTS(1)) defaults(clk,done[0]);
    compare_legacy #(.A(16),.D(64),.I(3),.IDLE(0),.G(5),.L(15),.ID(5),.H(12),.DEFAULT_H(37)) wide(clk,done[1]);
    compare_legacy #(.A(32),.D(128),.I(2),.G(9),.L(16),.ID(2),.H(17),.DEFAULT_H(129)) widest(clk,done[2]);
    compare_legacy #(.A(16),.D(32),.I(2),.IDLE(0),.G(5),.L(15),.ID(1),.H(1),.DEFAULT_H(1)) one_bit_holdoff(clk,done[3]);
`endif
    initial begin
        wait (&done);
        $display("BENCH_OK");
        $finish;
    end
endmodule
`undef SUBJECT
`undef BASELINE
