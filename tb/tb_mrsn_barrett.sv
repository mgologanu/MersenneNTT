module tb_mrsn_crt;

  // Parameters

   parameter WIDTH = 32;

   parameter W0 = 13;

   parameter W1 = 19;


   parameter CLK_PERIOD = 10;

   localparam DELAY = 7;
   

   // Clock and reset signals
   logic    clk;
   logic    rst_n;
   logic    en;
   
   
   // DUT signals
  

   
   longint  prime0, prime1, ql, Rql, rl;
   

   longint zq_expected,  zq_expected_q;
 
   logic [31:0]		 rand_val1;

   logic [WIDTH-1:0]	 q, Rq;
   
   logic signed [WIDTH-1:0] zq, r, r_q;
   
   
   logic [W0-1:0] t0;
   
   logic [W1-1:0] t1;

   
  // Instantiate the DUT
  mrsn_barrett
    #(
      .WIDTH(WIDTH)
      ) 
   dut (
	.clk_i(clk),
	.rst_ni(rst_n),
	.en_i(en),
	.q(q),
	.Rq(Rq),
	.r(r),
	.zq(zq)
	);


   // Delayed values

   
   pipe_reg #(.WIDTH(32), .DEPTH(DELAY)) i_r4 (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(r), .output_o(r_q));
   pipe_reg #(.WIDTH(64), .DEPTH(DELAY)) i_r5 (.clk_i(clk), .rst_ni(rst_n), .en_i(en), .input_i(zq_expected), .output_o(zq_expected_q));
   
   
   function automatic longint signed_reduction (longint rl,  longint ql);

      longint signed z, z1, z2;


      z1 = rl % ql;
      

      if (z1 > (ql-1)/2) begin
	 z2 = z1 - ql;
      end else begin
	 z2 = z1;
      end

      
      if (z2 < -(ql-1)/2) begin
	 z = z2 + ql;
      end else begin
	 z = z2;
      end

      
      return z;
      
   endfunction


   

   
   
  // Clock generation


   initial clk = 1'b0;
   always #(CLK_PERIOD / 2) clk = ~clk;
   
   initial prime0 = (1 << W0) - 1;

   initial prime1 = (1 << W1) - 1;

   initial q = 32'd3329;

   initial Rq = 32'd2580335;

   initial ql = 3329;

   initial Rql = 2580335;
   

   
   initial begin
      
      rst_n = 1;
      @(posedge clk);
      rst_n = 0;
      @(posedge clk);
      rst_n = 1;


      @(posedge clk);
      en = 1'b1;

      
  
      @(posedge clk);
      en = 1'b1;
      
      
      t0 = 7345;
      t1 = 405557;

      rl = {51'b0,t0} + {45'b0,t1} * prime0;

      
      if (rl > (prime0*prime1-1)/2) begin
	rl = rl - prime0*prime1;
      end

      r = rl[32-1:0];
      
      zq_expected = signed_reduction(rl, ql);


      
      rand_val1    = 32'hDEAD_BEEF;
    
       
      for (integer i = 0; i < 100; i = i + 1) begin

	 @(posedge clk);


	 // Simple LFSR-style pseudo-random
         rand_val1 = rand_val1 ^ (rand_val1 << 13);
         rand_val1 = rand_val1 ^ (rand_val1 >> 17);
         rand_val1 = rand_val1 ^ (rand_val1 << 5);

	 t0 = rand_val1[W0-1:0];
	 t1 = rand_val1[W1+W0-1:W0];
	 
	 rl = {51'b0,t0} + {45'b0,t1} * prime0;

      
	 if (rl > (prime0*prime1-1)/2) begin
	    rl = rl - prime0*prime1;
	 end
	 
	 r = rl[32-1:0];
	 
	 zq_expected = signed_reduction(rl, ql);
	 
      end // for (i = 0; i < 10; i = i + 1)
      
      
      repeat (DELAY) @(posedge clk);

      @(posedge clk);
      
    $display("=== TEST COMPLETE ===");
    $finish;

   end  // initial begin

   
   always @(posedge clk) begin
      $display("[%0t]  Input: %0d, Output: %0d, Expected: %0d, Error: %0d",
	       $time, 
	       r_q, zq, zq_expected_q,
	       zq - zq_expected_q[32-1:0]
	       );
    end

 
 

endmodule
