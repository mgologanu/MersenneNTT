module tb_mrsn_ntt_complex16_1_p13;

  
   parameter LEN          = 16;
   
   parameter CLK_PERIOD = 10;

   localparam DUT_LATENCY = 8;

   parameter  W = 13;

   

   logic     clk;
   logic     rst_n;
   logic     en;

   logic [W-1:0] a_re[0:LEN-1];
   
   logic [W-1:0] a_im[0:LEN-1];
   
   logic [W-1:0] c_re[0:LEN-1];
   
   logic [W-1:0] c_im[0:LEN-1];
   
   
   logic [W-1:0] x_re_array [0:3][0:LEN-1];
   logic [W-1:0] x_im_array [0:3][0:LEN-1];
   
   
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
      x_re_array[0] = '{
			13'd2600, 13'd2980, 13'd1480, 13'd3429, 
			13'd6608, 13'd3328, 13'd5255, 13'd5870, 
			13'd2168, 13'd454, 13'd2764, 13'd1978, 
			13'd826, 13'd551, 13'd3783, 13'd7693
			};
      x_im_array[0] = '{
			13'd7661, 13'd3262, 13'd3695, 13'd1096, 
			13'd8104, 13'd311, 13'd1626, 13'd4582, 
			13'd2052, 13'd1645, 13'd6541, 13'd5549, 
			13'd7340, 13'd4572, 13'd7779, 13'd1445
			};

      /*
       
       Expected output:
       
       X_re =  2621  7392  1802  4229  6605  2783  4179   342   266  3647   456  2494  5835  7252  3181  4898

       X_im =  1732  5954  2050   791  7129  4585  2648  3669  6168  5278  3203  6342  1260  4458  1754    27

	*/

      
      x_re_array[1] = '{13'd2447, 13'd2246, 13'd8152, 13'd1753, 13'd4591, 13'd4125, 13'd3564, 13'd5735, 13'd2368, 13'd3332, 13'd3544, 13'd7251, 13'd3532, 13'd1185, 13'd8092, 13'd3118};
      
      x_im_array[1] = '{13'd6040, 13'd5613, 13'd2626, 13'd6536, 13'd4836, 13'd225, 13'd6264, 13'd748, 13'd3904, 13'd662, 13'd4515, 13'd6255, 13'd5265, 13'd663, 13'd4869, 13'd660};
      
      /*
       Expected ouput:
       
       X_re =  7698  7545  4813  7123  7913  1646  1335   447  3224  1782   114  5103  1182  6700  6638   462

       X_im =  2344   575  2993   549  3359  4598   813  6984  1354  4474  2721  4231  3691  1994   715  6099
       
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
