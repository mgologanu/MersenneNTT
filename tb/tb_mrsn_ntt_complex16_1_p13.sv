module tb_mrsn_ntt_complex16_1_p13;

  
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
   
   mrsn_ntt_complex16_1
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
      
      x_re_array[0] = '{7693,3783,551,826,1978,2764,454,2168,5870,5255,3328,6608,3429,1480,2980,2600};
      
      x_im_array[0] = '{1445,7779,4572,7340,5549,6541,1645,2052,4582,1626,311,8104,1096,3695,3262,7661};
      
      /*
       
       Input:
  
       2600  2980  1480  3429  6608  3328  5255  5870  2168  454  2764  1978  826  551  3783  7693     
       7661  3262  3695  1096  8104  311  1626  4582  2052  1645  6541  5549  7340  4572  7779  1445
  
       
       Expected output:
       
       2621  7392  1802  4229  6605  2783  4179   342   266  3647   456  2494  5835  7252  3181  4898
       1732  5954  2050   791  7129  4585  2648  3669  6168  5278  3203  6342  1260  4458  1754    27

	*/

      
      x_re_array[1] = '{3118,8092,1185,3532,7251,3544,3332,2368,5735,3564,4125,4591,1753,8152,2246,2447};
      
      x_im_array[1] = '{660,4869,663,5265,6255,4515,662,3904,748,6264,225,4836,6536,2626,5613,6040};
            
      
      /*
       
       Input:
       
       2447  2246  8152  1753  4591  4125  3564  5735  2368  3332  3544  7251  3532  1185  8092  3118
       6040  5613  2626  6536  4836  225  6264  748  3904  662  4515  6255  5265  663  4869  660
        
       Expected ouput:
       
       7698  7545  4813  7123  7913  1646  1335   447  3224  1782   114  5103  1182  6700  6638   462
       2344   575  2993   549  3359  4598   813  6984  1354  4474  2721  4231  3691  1994   715  6099
       
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
     $display("[%0t]  \nInput:\n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d\nOutput:\n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d \n%0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d  %0d\n",
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
