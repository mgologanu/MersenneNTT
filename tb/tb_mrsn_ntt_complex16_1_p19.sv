module tb_mrsn_ntt_complex16_1_p19;

  
   parameter LEN          = 16;
   
   parameter CLK_PERIOD = 10;

   localparam DUT_LATENCY = 8;

   parameter  W = 19;

   

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
      
      x_re_array[0] = '{520007,420806,484880,422202,283031,439217,323063,122800,402625,2439,211700,192120,193509,257994,206611,420818};
      
      x_im_array[0] = '{263837,20081,162307,14876,498381,390509,270204,29978,251937,294119,13442,448368,373911,233576,315709,39363};
      
      /* 
       
       Input:
       
       420818  206611  257994  193509  192120  211700  2439  402625  122800  323063  439217  283031  422202  484880  420806  520007
       39363  315709  233576  373911  448368  13442  294119  251937  29978  270204  390509  498381  14876  162307  20081  263837
       
       Expected Output:

       185239  177257  139601  459654  106358  181038   51792  426570   33129  456988   46125  494723  245646  113172  383554   86520
       474876  369716  469956  291505  276828    7585  369814  391596  449122  253325  421179   90734  326658  327709  291360   12141
       */

      
      x_re_array[1] = '{50279,275117,149907,489331,49244,405405,148463,427684,190834,391437,458429,430019,444060,20280,279406,118327};
      
      x_im_array[1] = '{243113,30384,403571,114176,108919,169527,13482,108961,486043,452285,490114,513846,152194,506017,236682,467429};
      
      /* 
       
       Input:
       
       118327  279406  20280  444060  430019  458429  391437  190834  427684  148463  405405  49244  489331  149907  275117  50279
       467429  236682  506017  152194  513846  490114  452285  486043  108961  13482  169527  108919  114176  403571  30384  243113
       
       Expected Output:
            
       133926  262691  219542    2415  479951  484769  360691  326955  191880  161170  252041  228810   79018   22238   116165  143831
       302447  228507  347987  268698   62725  400847  403494  499267  344156  343576  326576  182316  366205   45778  339214  395636
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
