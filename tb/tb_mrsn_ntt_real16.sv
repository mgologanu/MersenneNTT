module tb_mrsn_ntt_real16;

  

   parameter WIDTH = 32;


   parameter LEN          = 16;
   
   parameter CLK_PERIOD = 10;

   localparam DELAY = 11;

   parameter  W0 = 13;

   parameter  W1  = 19;
   

   logic     clk;
   logic     rst_n;
   logic     en;

   logic signed [WIDTH-1:0] a [0:LEN-1];

   logic signed [WIDTH-1:0] x_array [0:3][0:LEN-1];
   

   logic  [WIDTH-1:0] c [0:LEN-1];

   logic  [W0-1:0] c0 [0:LEN-1];

   logic  [W1-1:0] c1 [0:LEN-1];
				      

   genvar		    i;
   
   mrsn_ntt_real16 
     #(
       .WIDTH(WIDTH),
       .LEN(LEN)
       ) 
   dut
     (
      .clk_i(clk),
      .rst_ni(rst_n),
      .en_i(en),
      .a_i(a),
      .c_o(c)
      );
   
    

  generate
     for (i = 0; i < LEN; i++) begin
	assign	c0[i] = c[i][W0-1:0];
	assign  c1[i] = c[i][W1+W0-1:W0];
     end
  endgenerate
   // Clock generation
   
   
   initial clk = 1'b0;
   always #(CLK_PERIOD / 2) clk = ~clk;

   
   initial begin
      x_array[0] = '{
		     32'sh00AA962B, 32'sh00951398, 32'sh00858979, 32'shFFF5F2F9,
		     32'shFFFDF071, 32'sh0090BEA5, 32'shFF17023A, 32'shFF434766,
		     32'sh004F87B6, 32'sh00A9F864, 32'sh00213BBF, 32'shFF2B54C6,
		     32'shFF3A4BF0, 32'shFFB8E78C, 32'shFFD241D8, 32'sh0030B112
		     };

      /* Expected:
	  Input:  11179563  9769880  8751481  -658695  -135055  9485989  -15269318  -12368026  5212086  11139172  2177983  -13937466  -12956688  -4659316  -2997800  3191058 
          Output_13: 4089  4018  3734  707  3373  643  151  5761  5498  6292  4664  2146  1970  2369  69  7135 
	  Output_19: 121822  304199  123546  67899  95186  106530  127579  70039  477129  47365  169259  196621  211214  164583  30553  369353
       */
      
      x_array[1] = '{
		     32'sh000BF5D0, 32'sh004A50F6, 32'sh008100AD, 32'shFF434743,
		     32'sh00F691FA, 32'shFF00585E, 32'sh00677828, 32'sh00077EF7,
		     32'shFF96434F, 32'sh00784DD7, 32'shFF976CDF, 32'shFF55D2DC,
		     32'shFF1A48E9, 32'sh00AA5F66, 32'sh00F8C844, 32'sh00FC30E4
		     };


      /* Expected:
	    Input:  783824  4870390  8454317  -12368061  16159226  -16754594  6780968  491255  -6929585  7884247  -6853409  -11152676  -15054615  11165542  16304196  16527588 
	    Output_13: 2114  883  5006  1088  1660  1246  2609  7501  3506  1010  6492  4395  6413  7032  1250  1426 
	    Output_19: 106450  95730  17372  138771  71486  77368  475298  226212  142664  337043  41444  11913  491633  455481  205662  368085
       */
      
        // Reset & enable
      rst_n  = 1'b0;
      en = 1'b0;
      repeat (2) @(posedge clk);
      rst_n  = 1'b1;
      en = 1'b1;
      
      @(posedge clk);
      a = x_array[0];

      @(posedge clk);
      a = x_array[1];


      
      repeat (DELAY) @(posedge clk);

      
      @(posedge clk);
      
      $display("=== TEST COMPLETE ===");
      $finish;
      
   end  // initial begin

  always @(posedge clk) begin
     $display("[%0t]  Input:  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n Output_13: %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n Output_19: %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d\n",
	      $time,
	      a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7],
	      a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15],
	      c0[0], c0[1], c0[2], c0[3], c0[4], c0[5], c0[6], c0[7],
	      c0[8], c0[9], c0[10], c0[11], c0[12], c0[13], c0[14], c0[15],
	      c1[0], c1[1], c1[2], c1[3], c1[4], c1[5], c1[6], c1[7],
	      c1[8], c1[9], c1[10], c1[11], c1[12], c1[13], c1[14], c1[15]);

    // $display("[%0t]  %0d %0d %0d %0d", $time, dut.c1_re_s1_q[0], dut.tmp1_re_s1[0], dut.c1_re_s2a[0], {dut.tmp1_re_s1[0][9:0], dut.tmp1_re_s1[0][19-1:10]});
     
   end

 
 

endmodule
