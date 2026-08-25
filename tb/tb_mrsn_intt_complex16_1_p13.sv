module tb_mrsn_intt_complex16_1_p13;

  
   parameter LEN          = 16;
   
   parameter CLK_PERIOD = 10;

   localparam DUT_LATENCY = 8;

   parameter  W = 13;

   

   logic     clk;
   logic     rst_n;
   logic     en;

   logic [W-1:0] a_re[LEN-1:0];
   
   logic [W-1:0] a_im[LEN-1:0];
   
   logic [W-1:0] c_re[LEN-1:0];
   
   logic [W-1:0] c_im[LEN-1:0];
   
   
   logic [W-1:0] x_re_array [0:3][LEN-1:0];
   logic [W-1:0] x_im_array [0:3][LEN-1:0];
   
   
   genvar	 i;
   
   mrsn_intt_complex16_1
     #(
       .W(W),
       .LEN(LEN)
       ) 
   dut
     (
      .clk_i(clk),
      .rst_ni(rst_n),
      .en_i(en),
      .a_re_i(a_re),
      .a_im_i(a_im),
      .c_re_o(c_re),
      .c_im_o(c_im)
      );
   
    
   // Clock generation
   
   
   initial clk = 1'b0;
   always #(CLK_PERIOD / 2) clk = ~clk;

   
   initial begin 

      //initialize in reverse order
      
      x_re_array[0] = '{4898,3181,7252,5835,2494,456,3647,266,342,4179,2783,6605,4229,1802,7392,2621};
      
      x_im_array[0] = '{27,1754,4458,1260,6342,3203,5278,6168,3669,2648,4585,7129,791,2050,5954,1732};

      /*
       
       Input:
       2621  7392  1802  4229  6605  2783  4179  342  266  3647  456  2494  5835  7252  3181  4898
       1732  5954  2050  791  7129  4585  2648  3669  6168  5278  3203  6342  1260  4458  1754  27

       Expected Output:
       645   6725   7298   5718   7436   4102   2170   3819   1924   7264   3269   7075   5025    625   3191    223
       7902   3046   1783   1154   6799   4976   1443   7784     68   1747   6364   6874   2766   7624   1599   6738
 
      
        */
      x_re_array[1] = '{462,6638,6700,1182,5103,114,1782,3224,447,1335,1646,7913,7123,4813,7545,7698};
      
      x_im_array[1] = '{6099,715,1994,3691,4231,2721,4474,1354,6984,813,4598,3359,549,2993,575,2344};

      /*
       
       Input:
       7698  7545  4813  7123  7913  1646  1335  447  3224  1782  114  5103  1182  6700  6638  462
       2344  575  2993  549  3359  4598  813  6984  1354  4474  2721  4231  3691  1994  715  6099

       Expected Output:
       6388   3172   7567   3475   7928    472   7878   1659   5124   4166   7558   1342   7366   2578   6607    742
       6539   7898   1061   6284   3657   3600   1932   3777   5127   2401   6712   1788   2330   2417   4185   2369

	    */
      // Reset & enable
      rst_n  = 1'b0;
      en = 1'b0;
      
      repeat (2) @(posedge clk);
      
      rst_n  = 1'b1;
      en = 1'b1;
      
      @(posedge clk);
      a_re = x_re_array[0];
      a_im = x_im_array[0];

      @(posedge clk);
      a_re = x_re_array[1];
      a_im = x_im_array[1];


      
      repeat (DUT_LATENCY) @(posedge clk);

      
      @(posedge clk);
      
      $display("=== TEST COMPLETE ===");
      $finish;
      
   end  // initial begin

  always @(posedge clk) begin
     $display("[%0t]  Input:\n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d\nOutput:\n %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d\n",
	      $time,
	      a_re[0], a_re[1], a_re[2], a_re[3], a_re[4], a_re[5], a_re[6], a_re[7],
	      a_re[8], a_re[9], a_re[10], a_re[11], a_re[12], a_re[13], a_re[14], a_re[15],
	      a_im[0], a_im[1], a_im[2], a_im[3], a_im[4], a_im[5], a_im[6], a_im[7],
	      a_im[8], a_im[9], a_im[10], a_im[11], a_im[12], a_im[13], a_im[14], a_im[15],
	      c_re[0], c_re[1], c_re[2], c_re[3], c_re[4], c_re[5], c_re[6], c_re[7],
	      c_re[8], c_re[9], c_re[10], c_re[11], c_re[12], c_re[13], c_re[14], c_re[15],
	      c_im[0], c_im[1], c_im[2], c_im[3], c_im[4], c_im[5], c_im[6], c_im[7],
	      c_im[8], c_im[9], c_im[10], c_im[11], c_im[12], c_im[13], c_im[14], c_im[15]);

    // $display("[%0t]  %0d %0d %0d %0d\n  %0d %0d %0d %0d\n\n", $time, dut.c_re_tmp[12], dut.c_im_tmp[12], dut.c_re_tmp[14], dut.c_im_tmp[14],   dut.c_re_s4[12], dut.c_im_s4[12], dut.c_re_s4[14], dut.c_im_s4[14]);
     
     
   end

 
 

endmodule
