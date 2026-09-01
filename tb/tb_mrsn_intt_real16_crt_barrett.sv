module tb_mrsn_ntt_real16_crt_barrett;

   parameter WIDTH = 32;

   parameter LEN          = 16;
   
   parameter CLK_PERIOD = 10;

   localparam DELAY = 23;

   parameter  W0 = 13;

   parameter  W1  = 19;

   localparam W_NOTUSED = WIDTH - W0 - W1;

   logic     clk;
   logic     rst_n;
   logic     en;


   logic [WIDTH-1:0] q, Rq;
   
   logic  [WIDTH-1:0] c [LEN-1:0];

   logic  [W0-1:0] c0 [LEN-1:0];

   logic  [W1-1:0] c1 [LEN-1:0];


   logic signed [WIDTH-1:0] zq  [LEN-1 : 0];
   

	 
   logic signed [W0-1:0] c0_array [0:3][LEN-1:0];
   logic signed [W1-1:0] c1_array [0:3][LEN-1:0];
   
   genvar		    i;
   
   mrsn_intt_real16_crt_barrett
     #(
       .WIDTH(WIDTH),
       .LEN(LEN)
       ) 
   dut
     (
      .clk_i(clk),
      .rst_ni(rst_n),
      .en_i(en),
      .q(q),
      .Rq(Rq),
      .c_i(c),
      .zq(zq)
      );
   
    

  generate
     for (i = 0; i < LEN; i++) begin
	assign	c[i] =  {{W_NOTUSED{1'b0}}, c1[i], c0[i]};
     end
  endgenerate
   

   // Clock generation
   
   
   initial clk = 1'b0;
   always #(CLK_PERIOD / 2) clk = ~clk;

   initial q = 32'd3329;

   initial Rq = 32'd2580335;


   
   initial begin

      //initialize in reverse order
      
      c0_array[0] = '{7135,69,2369,1970,2146,4664,6292,5498,5761,151,643,3373,707,3734,4018,4089};
      
      c1_array[0] = '{369353,30553,164583,211214,196621,169259,47365,477129,70039,127579,106530,95186,67899,123546,304199,121822};
      
		     
      /*Expected:
       Input_13: 4089  4018  3734  707  3373  643  151  5761  5498  6292  4664  2146  1970  2369  69  7135 
       Input_19: 121822  304199  123546  67899  95186  106530  127579  70039  477129  47365  169259  196621  211214  164583  30553  369353


       Output_13:  7166  518  3371  5444  772  6488  6030  3072  4498  3487  1607  4355  3601  2713  848  5308 
       Output_19:  307714  40277  281677  497597  492421  390584  4327  146035  278015  508873  122393  173403  155322  474136  134802  362688
       */

      
      c0_array[1] = '{1426,1250,7032,6413,4395,6492,1010,3506,7501,2609,1246,1660,1088,5006,883,2114};
      
      c1_array[1] = '{368085,205662,455481,491633,11913,41444,337043,142664,226212,475298,77368,71486,138771,17372,95730,106450};

      /* Expected
       Input_13: 2114  883  5006  1088  1660  1246  2609  7501  3506  1010  6492  4395  6413  7032  1250  1426 
       Input_19: 106450  95730  17372  138771  71486  77368  475298  226212  142664  337043  41444  11913  491633  455481  205662  368085
   
       Output_13:  4477  6724  1449  2792  3446  772  6942  6551  8  3276  3282  3155  3544  1481  84  1582 
       Output_19:  503435  165882  1513  145755  299206  180720  246183  260031  137742  159536  222863  431669  149090  195546  410392  100380
      */

      
      // Reset & enable
      rst_n  = 1'b0;
      en = 1'b0;
      repeat (2) @(posedge clk);
      rst_n  = 1'b1;
      en = 1'b1;
      
      @(posedge clk);
      c0 = c0_array[0];
      c1 = c1_array[0];

      @(posedge clk);
      c0 = c0_array[1];
      c1 = c1_array[1];

 
      repeat (DELAY) @(posedge clk);

      @(posedge clk);
      
      $display("=== TEST COMPLETE ===");
      $finish;
      
   end  // initial begin

  always @(posedge clk) begin
     $display("[%0t]  \nInput_13: %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n Input_19: %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d\n Output:  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n ",
	      $time,
	      c0[0], c0[1], c0[2], c0[3], c0[4], c0[5], c0[6], c0[7],
	      c0[8], c0[9], c0[10], c0[11], c0[12], c0[13], c0[14], c0[15],
	      c1[0], c1[1], c1[2], c1[3], c1[4], c1[5], c1[6], c1[7],
	      c1[8], c1[9], c1[10], c1[11], c1[12], c1[13], c1[14], c1[15],
	      zq[0], zq[1], zq[2], zq[3], zq[4], zq[5], zq[6], zq[7],
	      zq[8], zq[9], zq[10], zq[11], zq[12], zq[13], zq[14], zq[15]
	    
	      );

    // $display("[%0t]  %0d %0d %0d %0d", $time, dut.c0_re_s2[7], dut.c0_im_s2[7], dut.c0_re_s2a[7], dut.c0_im_s2a[7]);
     
     
   end

 
 

endmodule
