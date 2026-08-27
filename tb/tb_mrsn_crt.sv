module tb_mrsn_crt;

  // Parameters

   parameter WIDTH = 32;
   
   parameter W0 = 13;

   parameter W1 = 19;

   parameter CLK_PERIOD = 10;

   localparam DELAY = 5;
   

   // Clock and reset signals
   logic    clk;
   logic    rst_n;
   logic    en;
   
   
   // DUT signals
  
   logic [W0-1:0] r0, r0_q;

   logic [W1-1:0] r1, r1_q;

   
   logic signed [W0-1:0] u0;

   longint u0_expected,  u0_expected_q;
   
   logic signed [W1-1:0] u1;


   longint		 u1_expected,  u1_expected_q;
      
   longint		 prime0, prime1;

   logic [31:0]		 rand_val1;
   
   
  // Instantiate the DUT
  mrsn_crt
    #(
      .WIDTH(WIDTH)
      ) 
   dut (
	.clk_i(clk),
	.rst_ni(rst_n),
	.en_i(en),
	.r0(r0),
	.r1(r1),
	.v0(u0),
	.v1(u1)
	);


   // Delayed values

   
   pipe_reg #(.WIDTH(W0), .DEPTH(DELAY)) i_r0 (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(r0), .output_o(r0_q));
   pipe_reg #(.WIDTH(W1), .DEPTH(DELAY)) i_r1 (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(r1), .output_o(r1_q));

   
   
   pipe_reg #(.WIDTH(64), .DEPTH(DELAY)) i_r2 (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(u0_expected), .output_o(u0_expected_q));
   pipe_reg #(.WIDTH(64), .DEPTH(DELAY)) i_r3 (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(u1_expected), .output_o(u1_expected_q));

   
   function automatic longint signed_newton0 (longint r0, longint r1, longint p0, longint p1);

      longint		z;
      
      if (r0 > (p0-1)/2) begin
	 z = r0 - p0;
      end else begin
	 z = r0;
      end
      return z;
   endfunction

   
   
   function automatic longint signed_newton1 (longint r0, longint r1, longint p0, longint p1);

      longint zz, z1, z2, z3, z, t;

   

      z1 = (r1-r0) % p1;

 

      t = z1 * 8321;

      
      z2 = t % p1;

      if (r0 > (p0-1)/2) begin
	 z3 = (z2 + 1) % p1;
      end else begin
	 z3 = z2;
      end
      
      if (z3 > (p1-1)/2) begin
	 z = z3 -p1;
      end else begin
	 z = z3;
      end
      
      return z3;
      
   endfunction


   
   
  // Clock generation


   initial clk = 1'b0;
   always #(CLK_PERIOD / 2) clk = ~clk;
   
   initial prime0 = (1 << W0) - 1;
   initial prime1 = (1 << W1) - 1;

   
   initial begin
      
      rst_n = 1;
      @(posedge clk);
      rst_n = 0;
      @(posedge clk);
      rst_n = 1;


      @(posedge clk);
      en = 1'b1;
      
/* -----\/----- EXCLUDED -----\/-----
      r0 = 2345;
      r1 = 105557;

      u0_expected =  signed_newton0(r0, r1, prime0, prime1);
      u1_expected =  signed_newton1(r0, r1, prime0, prime1);

      @(posedge clk);
      en = 1'b1;
      
      r0 = 7345;
      r1 = 205557;

      u0_expected =  signed_newton0(r0, r1, prime0, prime1);
      u1_expected =  signed_newton1(r0, r1, prime0, prime1);
 -----/\----- EXCLUDED -----/\----- */
/* -----\/----- EXCLUDED -----\/-----

      rand_val1    = 32'hDEAD_BEEF;
      rand_val2    = 32'hBEEF_DEAD;
       
      for (integer i = 0; i < 100; i = i + 1) begin

	 @(posedge clk);


	 // Simple LFSR-style pseudo-random
         rand_val1 = rand_val1 ^ (rand_val1 << 13);
         rand_val1 = rand_val1 ^ (rand_val1 >> 17);
         rand_val1 = rand_val1 ^ (rand_val1 << 5);

	 r0 = rand_val1[W0-1:0];
	 r1 = rand_val1[W1+W0-1:W0];
	 
	 u0_expected =  signed_newton0(r0, r1, prime0, prime1);
	 u1_expected =  signed_newton1(r0, r1, prime0, prime1);

	 
      end // for (i = 0; i < 10; i = i + 1)
      
 -----/\----- EXCLUDED -----/\----- */
      
      
      
      @(posedge clk);
      en = 1'b1;
      
      r0 = 8191;
      r1 = 524287;

      u0_expected =  signed_newton0({51'b0,r0}, {45'b0,r1}, prime0, prime1);
      u1_expected =  signed_newton1({51'b0,r0}, {45'b0,r1}, prime0, prime1);
      

      
      
      
      repeat (DELAY) @(posedge clk);

      @(posedge clk);
      
    $display("=== TEST COMPLETE ===");
    $finish;

   end  // initial begin

   
   always @(posedge clk) begin
      $display("[%0t]  Input: %0d, %0d , Output: %0d, %0d, Expected: %0d, %0d,  Errors: %0d %0d",
	       $time, r0_q, r1_q, u0, u1, u0_expected_q, u1_expected_q,
	               $signed(u0) - u0_expected[13:0], $signed(u1) - u1_expected[31:0]
	      	       );
    end

 
 

endmodule
