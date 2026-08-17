//`timescale 1ns/1ps

module tb_mrsn_omega16_p19;

  // Parameters
  parameter N = 19;
  
   parameter CLK_PERIOD = 10;

   parameter C = 162504;

   parameter S = 42613;

   parameter DELAY = 3;
   

  // Clock and reset signals
  logic clk;
  logic rst_n;
   logic en;

  // DUT signals
  logic [N-1:0] a, a_q;
  logic [N-1:0] a_c, a_s;

  logic mode;

  longint prime;

  // Expected sum for verification
  logic [N-1:0] expected_a_c, expected_a_c_q;
  logic [N-1:0] expected_a_s, expected_a_s_q;
   

  // Instantiate the DUT
  mrsn_omega16_p19 #(
      .W(N)
  ) dut (
      .clk_i(clk),
      .rst_ni(rst_n),
      .en_i(en),
      .x_i(a),
      .x_c_o(a_c),
      .x_s_o(a_s)
  );

  // Clock generation

  initial clk = 1'b0;
  always #(CLK_PERIOD / 2) clk = ~clk;

  initial prime = (1 << N) - 1;

   pipe_reg #(.WIDTH(N), .DEPTH(DELAY)) i_a (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(a), .output_o(a_q));
   pipe_reg #(.WIDTH(N), .DEPTH(DELAY)) i_a_c (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(expected_a_c), .output_o(expected_a_c_q));
   pipe_reg #(.WIDTH(N), .DEPTH(DELAY)) i_a_s (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(expected_a_s), .output_o(expected_a_s_q));

   
  

  initial begin

    rst_n = 1;
    @(posedge clk);
    rst_n = 0;
    @(posedge clk);
    rst_n = 1;

     @(posedge clk);
     en = 1'b1;
     a = 5623;
  
     expected_a_c = (a*C) % prime;
     expected_a_s = (a*S) % prime;

     
     @(posedge clk);
     a = 1234;
  
     expected_a_c = (a*C) % prime;
     expected_a_s = (a*S) % prime;

     
     @(posedge clk);
     a = 8190;
  
     expected_a_c = (a*C) % prime;
     expected_a_s = (a*S) % prime;

     
     @(posedge clk);
     a = prime-1;
  
     expected_a_c = (a*C) % prime;
     expected_a_s = (a*S) % prime;

     
     @(posedge clk);
     a = 0;
  
     expected_a_c = (a*C) % prime;
     expected_a_s = (a*S) % prime;
     
     repeat (4) @(posedge clk);

    $display("=== TEST COMPLETE ===");
    $finish;

  end  // initial begin


  always @(posedge clk) begin
    $display("[%0t]  Input: %0d Prod0: %0d Expected: %0d, Prod1: %0d Expected: %0d  Error: %0d %0d ", $time, a_q, a_c, expected_a_c_q, a_s, expected_a_s_q, $signed(a_c-expected_a_c_q), $signed(a_s-expected_a_s_q) ) ;
  end



endmodule
