module tb_mrsn_intt_complex16_1_p13;

  
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

      x_re_array[0] = '{13'd2621, 13'd7392, 13'd1802, 13'd4229, 13'd6605, 13'd2783, 13'd4179, 13'd342, 13'd266, 13'd3647, 13'd456, 13'd2494, 13'd5835, 13'd7252, 13'd3181, 13'd4898};
      
      x_im_array[0] = '{13'd1732, 13'd5954, 13'd2050, 13'd791, 13'd7129, 13'd4585, 13'd2648, 13'd3669, 13'd6168, 13'd5278, 13'd3203, 13'd6342, 13'd1260, 13'd4458, 13'd1754, 13'd27};

      x_re_array[1] = '{13'd7698, 13'd7545, 13'd4813, 13'd7123, 13'd7913, 13'd1646, 13'd1335, 13'd447, 13'd3224, 13'd1782, 13'd114, 13'd5103, 13'd1182, 13'd6700, 13'd6638, 13'd462};
      
      x_im_array[1] = '{13'd2344, 13'd575, 13'd2993, 13'd549, 13'd3359, 13'd4598, 13'd813, 13'd6984, 13'd1354, 13'd4474, 13'd2721, 13'd4231, 13'd3691, 13'd1994, 13'd715, 13'd6099};
      
      
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
