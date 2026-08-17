//`timescale 1ns/1ps

module tb_mrsn_add_sub;

  // Parameters
  parameter N = 19;
  parameter CLK_PERIOD = 10;

  // Clock and reset signals
  logic clk;
  logic rst_n;
  logic en;
   

  // DUT signals
  logic [N-1:0] a, a_q;
  logic [N-1:0] b, b_q;
  logic [N-1:0] sum;

  logic mode, mode_q;

  longint prime;

  // Expected sum for verification
  logic [N-1:0] expected_sum, expected_sum_q;

  // Instantiate the DUT
  mrsn_add_sub #(
      .N(N)
  ) dut (
      .clk_i(clk),
      .rst_ni(rst_n),
      .en_i(en),
      .mode_i(mode),
      .a_i(a),
      .b_i(b),
      .sum_o(sum)
  );

  // Clock generation

  initial clk = 1'b0;
  always #(CLK_PERIOD / 2) clk = ~clk;

  initial prime = (1 << N) - 1;

  always @(posedge clk) begin
    a_q <= a;

    b_q <= b;

    expected_sum_q <= expected_sum;

     mode_q = mode;
     
  end


  initial begin

    rst_n = 1;
    @(posedge clk);
    rst_n = 0;
    @(posedge clk);
    rst_n = 1;

     
     @(posedge clk);
     en = 1;
     a = 0;
     b = 0;
     mode = 0;
    expected_sum =  (longint'(a) + longint'(b)) % prime;
    @(posedge clk);
    a = 35;
    b = 120;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 74762;
    b = 349525;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 362096;
    b = 131248;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 82;
    b = 31;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 362096;
    b = 431248;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 35;
    b = 120;
    mode = 1;
    expected_sum =  (longint'(a) + prime - longint'(b)) % prime;

    @(posedge clk);
    a = 74762;
    b = 349525;
    mode = 1;
    expected_sum =  (longint'(a) + prime - longint'(b)) % prime;

    @(posedge clk);
    a = 362096;
    b = 131248;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 82;
    b = 31;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 362096;
    b = 431248;
    mode = 1;
    expected_sum =  (longint'(a) + prime - longint'(b)) % prime;

    @(posedge clk);
    a = 362096;
    b = 362096;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 0;
    b = 0;
    mode = 1;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 362096;
    b = 0;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 0;
    b = 362096;
    mode = 1;
    expected_sum =  (longint'(a) + prime - longint'(b)) % prime;

    @(posedge clk);
    a = 524287;
    b = 0;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 0;
    b = 524287;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    @(posedge clk);
    a = 524287;
    b = 0;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 0;
    b = 524287;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 524286;
    b = 0;
    mode = 1;
    expected_sum = (longint'(a) - longint'(b)) % prime;

    @(posedge clk);
    a = 0;
    b = 1;
    mode = 1;
    expected_sum =  (longint'(a) + prime - longint'(b)) % prime;


     
    @(posedge clk);
    a = 524287;
    b = 524287;
    mode = 0;
     expected_sum = (longint'(a) + longint'(b)) % prime;


     
     @(posedge clk);
    a = 524286;
     b = 524287;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;
     
    @(posedge clk);
    a = 0;
    b = 0;
    mode = 0;
    expected_sum = (longint'(a) + longint'(b)) % prime;

    repeat (1) @(posedge clk);

    $display("=== TEST COMPLETE ===");
    $finish;

  end  // initial begin


  always @(posedge clk) begin
     if ($signed(longint'(sum) - longint'(expected_sum_q)) != 0) begin
     
	$display("[%0t]  Input: %0d %0d %0d Sum: %0d Expected: %0d Error: %0d", $time, a_q, b_q, mode_q, sum,
		 expected_sum_q, $signed(longint'(sum) - longint'(expected_sum_q)));
     end
  end



endmodule
