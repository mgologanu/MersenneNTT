module tb_mrsn_complex_multiply;

  import modular_complex_pkg::*;

  // Parameters

   parameter WIDTH = 32;
   
   parameter W0 = 13;

   parameter W1 = 19;

   parameter W_NOTUSED = WIDTH-W0-W1;
   
   parameter CLK_PERIOD = 10;

   localparam DELAY = 5;
   

   // Clock and reset signals
   logic    clk;
   logic rst_n;
   logic en;
   
   
   // DUT signals
   logic [WIDTH-1:0]	a_re, a_im, b_re, b_im, z_re, z_im;

   logic [W0-1:0]	a0_re, a0_im, b0_re, b0_im, z0_re, z0_im, z0_re_expected, z0_im_expected;
   logic [W0-1:0]	a0_re_q, a0_im_q, b0_re_q, b0_im_q, z0_re_expected_q, z0_im_expected_q;

   logic [W1-1:0]	a1_re, a1_im, b1_re, b1_im, z1_re, z1_im, z1_re_expected, z1_im_expected;
   logic [W1-1:0]	a1_re_q, a1_im_q, b1_re_q, b1_im_q, z1_re_expected_q, z1_im_expected_q;
   
   complex_mod_t  a0, b0, a1, b1;
   
      
   longint		prime0, prime1;

   logic [31:0]		rand_val1, rand_val2, rand_val3, rand_val4;
   
  // Instantiate the DUT
  mrsn_complex_multiply2 
    #(
      .WIDTH(WIDTH),
      .W0(W0),
      .W1(W1)
      ) 
   dut (
	.clk_i(clk),
	.rst_ni(rst_n),
	.en_i(en),
	.a_re_i({{W_NOTUSED{1'b0}}, a1_re, a0_re}),
	.a_im_i({{W_NOTUSED{1'b0}}, a1_im, a0_im}),
	.b_re_i({{W_NOTUSED{1'b0}}, b1_re, b0_re}),
	.b_im_i({{W_NOTUSED{1'b0}}, b1_im, b0_im}),
	.z_re_o(z_re),
	.z_im_o(z_im)
	);


   // Delayed values

   
   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_a0_re (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(a0_re), .output_o(a0_re_q));
   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_a0_im (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(a0_im), .output_o(a0_im_q));

   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_b0_re (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(b0_re), .output_o(b0_re_q));
   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_b0_im (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(b0_im), .output_o(b0_im_q));
   
   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_z0_re (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(z0_re_expected), .output_o(z0_re_expected_q));
   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_z0_im (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(z0_im_expected), .output_o(z0_im_expected_q));

      
   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_a1_re (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(a1_re), .output_o(a1_re_q));
   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_a1_im (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(a1_im), .output_o(a1_im_q));

   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_b1_re (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(b1_re), .output_o(b1_re_q));
   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_b1_im (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(b1_im), .output_o(b1_im_q));
   
   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_z1_re (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(z1_re_expected), .output_o(z1_re_expected_q));
   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_z1_im (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(z1_im_expected), .output_o(z1_im_expected_q));

   function automatic integer cmpl_mul_re (integer a_re, integer a_im, integer b_re, integer b_im, longint p);
      complex_mod_t a, b, z;
      a.re = a_re;
      a.im = a_im;
      b.re = b_re;
      b.im = b_im;
      z = cmpl_mult_mod(a, b, p);
      
      return z.re;
   endfunction

   function automatic integer cmpl_mul_im (integer a_re, integer a_im, integer b_re, integer b_im, longint p);
      complex_mod_t a, b, z;
      a.re = a_re;
      a.im = a_im;
      b.re = b_re;
      b.im = b_im;
      z = cmpl_mult_mod(a, b, p);
      
      return z.im;
   endfunction

   
   
  // Clock generation


   initial clk = 1'b0;
   always #(CLK_PERIOD / 2) clk = ~clk;
   
   initial prime0 = (1 << W0) - 1;
   initial prime1 = (1 << W1) - 1;

   assign z0_re = z_re[W0-1:0];
   assign z0_im = z_im[W0-1:0];
   
   assign z1_re = z_re[W1+W0-1:W0];
   assign z1_im = z_im[W1+W0-1:W0];
   
   initial begin
      
      rst_n = 1;
      @(posedge clk);
      rst_n = 0;
      @(posedge clk);
      rst_n = 1;


      @(posedge clk);
      en = 1'b1;
      
      a0_re = 2345;
      a0_im = 4576;

      b0_re = 1245;
      b0_im = 5555;
     
      a1_re = 16991;
      a1_im = 37634;
 
      b1_re = 12456;
      b1_im = 105557;

      z0_re_expected = cmpl_mul_re(a0_re, a0_im, b0_re, b0_im, prime0);
      z0_im_expected = cmpl_mul_im(a0_re, a0_im, b0_re, b0_im, prime0);
      
      z1_re_expected = cmpl_mul_re(a1_re, a1_im, b1_re, b1_im, prime1);
      z1_im_expected = cmpl_mul_im(a1_re, a1_im, b1_re, b1_im, prime1);

      @(posedge clk);
      a0_re = 1345;
      a0_im = 2576;

      b0_re = 245;
      b0_im = 6555;
     
      a1_re = 26991;
      a1_im = 47634;
 
      b1_re = 112456;
      b1_im = 5557;

      z0_re_expected = cmpl_mul_re(a0_re, a0_im, b0_re, b0_im, prime0);
      z0_im_expected = cmpl_mul_im(a0_re, a0_im, b0_re, b0_im, prime0);
      
      z1_re_expected = cmpl_mul_re(a1_re, a1_im, b1_re, b1_im, prime1);
      z1_im_expected = cmpl_mul_im(a1_re, a1_im, b1_re, b1_im, prime1);


      rand_val1    = 32'hDEAD_BEEF;
      rand_val2    = 32'hBEEF_DEAD;
      rand_val3    = 32'hBEEF_BEEF;
      rand_val4    = 32'hDEAD_DEAD;
      
      for (integer i = 0; i < 100; i = i + 1) begin

	 @(posedge clk);
	 // Simple LFSR-style pseudo-random
         rand_val1 = rand_val1 ^ (rand_val1 << 13);
         rand_val1 = rand_val1 ^ (rand_val1 >> 17);
         rand_val1 = rand_val1 ^ (rand_val1 << 5);

	 rand_val2 = rand_val2 ^ (rand_val2 << 13);
         rand_val2 = rand_val2 ^ (rand_val2 >> 17);
         rand_val2 = rand_val2 ^ (rand_val2 << 5);

	 rand_val3 = rand_val3 ^ (rand_val3 << 13);
         rand_val3 = rand_val3 ^ (rand_val3 >> 17);
         rand_val3 = rand_val3 ^ (rand_val3 << 5);
	 rand_val4 = rand_val4 ^ (rand_val4 << 13);
         rand_val4 = rand_val4 ^ (rand_val4 >> 17);
         rand_val4 = rand_val4 ^ (rand_val4 << 5);
	 
         a0_re = rand_val1[W0-1:0];
	 a1_re = rand_val1[W1+W0-1:W0];
	 
	 a0_im = rand_val2[W0-1:0];
	 a1_im = rand_val2[W1+W0-1:W0];

	 b0_re = rand_val3[W0-1:0];
	 b1_re = rand_val3[W1+W0-1:W0];

	 b0_im = rand_val4[W0-1:0];
	 b1_im = rand_val4[W1+W0-1:W0];

	 

	 if (a0_re >= prime0) a0_re = a0_re - prime0;
	 if (a0_im >= prime0) a0_im = a0_im - prime0;
	 if (b0_re >= prime0) b0_re = b0_re - prime0;
	 if (b0_im >= prime0) b0_im = b0_im - prime0;

	 if (a1_re >= prime1) a1_re = a1_re - prime1;
	 if (a1_im >= prime1) a1_im = a1_im - prime1;
	 if (b1_re >= prime1) b1_re = b1_re - prime1;
	 if (b1_im >= prime1) b1_im = b1_im - prime1;

	 z0_re_expected = cmpl_mul_re(a0_re, a0_im, b0_re, b0_im, prime0);
	 z0_im_expected = cmpl_mul_im(a0_re, a0_im, b0_re, b0_im, prime0);
	 
	 z1_re_expected = cmpl_mul_re(a1_re, a1_im, b1_re, b1_im, prime1);
	 z1_im_expected = cmpl_mul_im(a1_re, a1_im, b1_re, b1_im, prime1);
	 
      
      end // for (i = 0; i < 10; i = i + 1)
      
      
      
      
      repeat (DELAY) @(posedge clk);

      @(posedge clk);
      
    $display("=== TEST COMPLETE ===");
    $finish;

   end  // initial begin

   
   always @(posedge clk) begin
      $display("[%0t]  Input: %0d+j%0d , %0d+j%0d, Output: %0d+j%0d, Expected: %0d+j%0d, Input: %0d+j%0d , %0d+j%0d, Output: %0d+j%0d, Expected: %0d+j%0d, Errors: %0d %0d %0d %0d",
	       $time, 
	       a0_re_q, a0_im_q, b0_re_q, b0_im_q, z0_re, z0_im, z0_re_expected_q, z0_im_expected_q,
	       a1_re_q, a1_im_q, b1_re_q, b1_im_q, z1_re, z1_im, z1_re_expected_q, z1_im_expected_q,
	       ($signed(z0_re)-$signed(z0_re_expected_q)), ($signed(z0_im)-$signed(z0_im_expected_q)),
	       ($signed(z1_re)-$signed(z1_re_expected_q)), ($signed(z1_im)-$signed(z1_im_expected_q))	       
	       );
    end

 
 

endmodule
