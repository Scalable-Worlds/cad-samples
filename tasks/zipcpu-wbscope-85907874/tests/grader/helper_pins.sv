`timescale 1ns/1ps
module tb;
    reg clk=0, reset=1;
    always #5 clk=!clk;
    integer cycle=0;
    initial begin
        repeat (3) @(negedge clk);
        reset=0;
        for (cycle=0; cycle<180; cycle=cycle+1) @(negedge clk);
        $display("BENCH_OK");
        $finish;
    end

    genvar mode;
    generate for (mode=0; mode<3; mode=mode+1) begin : compare
        localparam WIDTH = 8 + mode;
        reg [WIDTH-1:0] word=0;
        reg offer=0, accept=0;
        wire input_ready;
`ifdef COMPARE_sfifo
        localparam LG=2;
        // During each period, first fill, then exchange while full, drain,
        // then exchange while empty. The option values separate both policies.
        wire wr=(cycle%20 < 8 || cycle%20 >= 16) && !reset;
        wire rd=(cycle%20 >= 5) && !reset;
        wire [WIDTH+LG+2:0] expected, actual;
        wire full0, empty0, full1, empty1;
        wire [LG:0] fill0, fill1;
        wire [WIDTH-1:0] data0, data1;
        baseline_sfifo #(WIDTH,LG,mode!=1,mode==0,mode==1)
            ref_dut(clk,reset,wr,word,full0,fill0,rd,data0,empty0);
        sfifo #(WIDTH,LG,mode!=1,mode==0,mode==1)
            dut(clk,reset,wr,word,full1,fill1,rd,data1,empty1);
        assign expected={full0,empty0,fill0,data0};
        assign actual={full1,empty1,fill1,data1};
        assign input_ready=1;
`elsif COMPARE_skidbuffer
        wire in0, out0, in1, out1;
        wire [WIDTH-1:0] data0, data1;
        wire [WIDTH+1:0] expected={in0,out0,data0}, actual={in1,out1,data1};
        baseline_skidbuffer #(mode==0,mode==1,mode==2,WIDTH,mode!=1)
            ref_dut(clk,reset,offer,in0,word,out0,accept,data0);
        skidbuffer #(mode==0,mode==1,mode==2,WIDTH,mode!=1)
            dut(clk,reset,offer,in1,word,out1,accept,data1);
        assign input_ready=in0;
`else
        wire in0, out0, trig0, in1, out1, trig1;
        wire [WIDTH-1:0] data0, data1;
        wire [WIDTH+2:0] expected={in0,out0,trig0,data0}, actual={in1,out1,trig1,data1};
        wire trigger=(cycle%29==10), encode=(cycle%31 < 19);
        baseline_axisrle #(WIDTH)
            ref_dut(clk,!reset,offer,in0,word[WIDTH-2:0],out0,accept,data0,trigger,encode,trig0);
        axisrle #(WIDTH)
            dut(clk,!reset,offer,in1,word[WIDTH-2:0],out1,accept,data1,trigger,encode,trig1);
        assign input_ready=in0;
`endif
        integer b;
        always @(posedge clk) begin
            if (!reset) begin
                for (b=0; b<$bits(expected); b=b+1)
                    if ((expected[b]===1'b0 || expected[b]===1'b1) && actual[b]!==expected[b])
                        $fatal(1,"Positional helper pin differs in mode %0d at cycle %0d",mode,cycle);
                accept <= cycle%11 > 3;
                if (!offer || input_ready) begin
                    offer <= cycle%13 != 4;
                    if (cycle%19 > 6) word <= word*5 + 3;
                end
            end
        end
    end endgenerate
endmodule
