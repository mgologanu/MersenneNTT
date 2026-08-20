module tb_mrsn_ntt_complex16_1_p19;

  
   parameter LEN          = 16;
   
   parameter CLK_PERIOD = 10;

   localparam DUT_LATENCY = 8;

   parameter  W = 19;

   

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
      
      x_re_array[0] = '{19'd420818, 19'd206611, 19'd257994, 19'd193509, 19'd192120, 19'd211700, 19'd2439, 19'd402625, 19'd122800, 19'd323063, 19'd439217, 19'd283031, 19'd422202, 19'd484880, 19'd420806, 19'd520007};

      x_im_array[0] = '{19'd39363, 19'd315709, 19'd233576, 19'd373911, 19'd448368, 19'd13442, 19'd294119, 19'd251937, 19'd29978, 19'd270204, 19'd390509, 19'd498381, 19'd14876, 19'd162307, 19'd20081, 19'd263837};

      /* expected:
       X_re =   185239  177257  139601  459654  106358  181038   51792  426570   33129  456988   46125  494723  245646  113172  383554   86520
       X_im =   474876  369716  469956  291505  276828    7585  369814  391596  449122  253325  421179   90734  326658  327709  291360   12141
       */
      
      x_re_array[1] = '{19'd118327, 19'd279406, 19'd20280, 19'd444060, 19'd430019, 19'd458429, 19'd391437, 19'd190834, 19'd427684, 19'd148463, 19'd405405, 19'd49244, 19'd489331, 19'd149907, 19'd275117, 19'd50279};
      
      x_im_array[1] = '{19'd467429, 19'd236682, 19'd506017, 19'd152194, 19'd513846, 19'd490114, 19'd452285, 19'd486043, 19'd108961, 19'd13482, 19'd169527, 19'd108919, 19'd114176, 19'd403571, 19'd30384, 19'd243113};
      

      /* expected
       X_re =  133926  262691  219542    2415  479951  484769  360691  326955  191880  161170  252041  228810   79018   22238   116165  143831
       X_im =  302447  228507  347987  268698   62725  400847  403494  499267  344156  343576  326576  182316  366205   45778  339214  395636
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
