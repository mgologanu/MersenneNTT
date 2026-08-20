//`timescale 1ns/1ps

module tb_bf2;

  import modular_complex_pkg::*;

  localparam CLK_PERIOD = 10;

  localparam W = 19;

  localparam P = (1 << W) - 1;

  localparam FACTOR = (1 << ((W - 1) / 2));

  localparam LATENCY = 1;

  logic clk;

  logic rst_n;

  logic en;

  longint prime;

  logic   [W-1:0] a_re_i;
  logic   [W-1:0] a_im_i;
  logic   [W-1:0] b_re_i;
  logic   [W-1:0] b_im_i;

  logic   [W-1:0] a_re_o;
  logic   [W-1:0] a_im_o;
  logic   [W-1:0] b_re_o;
  logic   [W-1:0] b_im_o;

  bf2 #(
      .W(W)
  ) dut (
      .clk_i (clk),
      .rst_ni(rst_n),
      .en_i  (en),
      .a_re_i(a_re_i),
      .a_im_i(a_im_i),
      .b_re_i(b_re_i),
      .b_im_i(b_im_i),
      .a_re_o(a_re_o),
      .a_im_o(a_im_o),
      .b_re_o(b_re_o),
      .b_im_o(b_im_o)
  );


  initial clk = 1'b0;

  always #(CLK_PERIOD / 2) clk = ~clk;

  initial prime = (1 << W) - 1;


  function automatic void compute_expected(
      input logic [W-1:0] ar, input logic [W-1:0] ai, input logic [W-1:0] br,
      input logic [W-1:0] bi, output logic [W-1:0] ar_exp, output logic [W-1:0] ai_exp,
      output logic [W-1:0] br_exp, output logic [W-1:0] bi_exp, input longint p);

    complex_mod_t a, b, a_exp, b_exp, z1, z2;

    a.re = ar;
    a.im = ai;
    b.re = br;
    b.im = bi;

    z1.re = 0;
    z1.im = 1;


    z2    = cmpl_mult_mod(b, z1, p);

    a_exp = cmpl_add_mod(a, z2, p);

    b_exp = cmpl_sub_mod(a, z2, p);

    ar_exp = a_exp.re;
    ai_exp = a_exp.im;
    br_exp = b_exp.re;
    bi_exp = b_exp.im;
  endfunction


   
  // --------------------------------------------------------------------------
  // Tasks 
  // --------------------------------------------------------------------------

  task automatic reset_dut();
    en     = 1'b0;
    a_re_i = '0;
    a_im_i = '0;
    b_re_i = '0;
    b_im_i = '0;
    rst_n  = 1'b0;
    repeat (2) @(posedge clk);
    rst_n = 1'b1;
    repeat (1) @(posedge clk);
  endtask

  task automatic run_vector(input logic [W-1:0] ar, 
			    input logic [W-1:0]	ai, 
			    input logic [W-1:0]	br,
                            input logic [W-1:0]	bi);
    logic [W-1:0] ar_exp;
    logic [W-1:0] ai_exp;
    logic [W-1:0] br_exp;
    logic [W-1:0] bi_exp;

    compute_expected(ar, ai, br, bi, ar_exp, ai_exp, br_exp, bi_exp, prime);

    // Apply inputs on a falling edge so they are stable before the next posedge.
    @(negedge clk);
    en     = 1'b1;
    a_re_i = ar;
    a_im_i = ai;
    b_re_i = br;
    b_im_i = bi;

    @(posedge clk);  // sample edge
    #1;

    if (LATENCY > 1) begin
      repeat (LATENCY - 1) @(posedge clk);
      #1;
    end

    if (a_re_o !== ar_exp || a_im_o !== ai_exp || b_re_o !== br_exp || b_im_o !== bi_exp) begin
      $display(   "Mismatch at time %0t:", $time);
      $display("  Input : a = (%0d,%0d), b = (%0d,%0d)", ar, ai, br, bi);
      $display("  Actual: a = (%0d,%0d), b = (%0d,%0d)", a_re_o, a_im_o, b_re_o, b_im_o);
      $display("  Expect: a = (%0d,%0d), b = (%0d,%0d)\n", ar_exp, ai_exp, br_exp, bi_exp);
      $fatal(1, "bf2 output mismatch");
    end

    en = 1'b0;
  endtask

  // --------------------------------------------------------------------------
  // Stimulus
  // --------------------------------------------------------------------------
  logic [W-1:0] ar, ai, br, bi;
  integer i;

  initial begin
    reset_dut();


    run_vector(1'd0, 1'd0, 1'd0, 1'd0);
    run_vector(1'd1, 1'd0, 1'd0, 1'd0);
    run_vector(1'd0, 1'd1, 1'd0, 1'd0);

    run_vector(1'd0, 1'd0, 1'd1, 1'd0);
    run_vector(1'd0, 1'd0, 1'd0, 1'd1);

    run_vector(P - 1, P - 1, P - 1, P - 1);

    for (i = 0; i < 2000; i++) begin
      ar = $urandom_range(0, P - 1);
      ai = $urandom_range(0, P - 1);
      br = $urandom_range(0, P - 1);
      bi = $urandom_range(0, P - 1);
      run_vector(ar, ai, br, bi);
    end


    $display("bf2 testbench PASSED for W=%0d, %0d random vectors.", W, 2000);
    $finish;
  end

endmodule
