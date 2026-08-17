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
      x_array[1] = '{
		     32'sh000BF5D0, 32'sh004A50F6, 32'sh008100AD, 32'shFF434743,
		     32'sh00F691FA, 32'shFF00585E, 32'sh00677828, 32'sh00077EF7,
		     32'shFF96434F, 32'sh00784DD7, 32'shFF976CDF, 32'shFF55D2DC,
		     32'shFF1A48E9, 32'sh00AA5F66, 32'sh00F8C844, 32'sh00FC30E4
		     };

      
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
