`include "mrsn_ntt.svh"
	      
module mrsn_ntt_complex16_1 #(
    parameter  W = 19,
    parameter  LEN   = 16

) (
    input logic		     clk_i,
    input logic		     rst_ni,
    input logic		     en_i,
    input logic  [W-1:0] a_re_i [LEN - 1 : 0],
    input logic  [W-1:0] a_im_i [LEN - 1 : 0],
    output logic [W-1:0] c_re_o [LEN - 1 : 0],
    output logic [W-1:0] c_im_o [LEN - 1 : 0]
);

  localparam P8 = (W - 1) / 2;
   
  localparam LAT_SEC_HALF = 4;
  
  logic [W-1:0] c_re_s0[LEN];
  logic [W-1:0] c_im_s0[LEN];
  
  logic [W-1:0] c_re_s1[LEN];
  logic [W-1:0] c_im_s1[LEN];

  logic [W-1:0] c_re_s2[LEN];
  logic [W-1:0] c_im_s2[LEN];

  logic [W-1:0] c_re_s3[LEN];
  logic [W-1:0] c_im_s3[LEN];

  logic [W-1:0] c_re_s4[LEN];
  logic [W-1:0] c_im_s4[LEN];

  logic [W-1:0] c_re_tmp[LEN];
  logic [W-1:0] c_im_tmp[LEN];
   
  logic [W-1:0] tmp0_re_c_s2, tmp0_re_s_s2, tmp0_im_c_s2, tmp0_im_s_s2,
		tmp1_re_c_s2, tmp1_re_s_s2, tmp1_im_c_s2, tmp1_im_s_s2,
		tmp2_re_c_s2, tmp2_re_s_s2, tmp2_im_c_s2, tmp2_im_s_s2,
		tmp3_re_c_s2, tmp3_re_s_s2, tmp3_im_c_s2, tmp3_im_s_s2;

  logic [W-1:0] tmp_re_10, tmp_im_10,  tmp_re_14, tmp_im_14;
  
  genvar i,j;

   generate
      for (i=0; i<LEN; i++) begin
	 assign c_re_s0[i] = a_re_i[i];
         assign c_im_s0[i] = a_im_i[i];
	 assign c_re_o[i] = c_re_s4[i];
         assign c_im_o[i] = c_im_s4[i];
	 
      end
   endgenerate


   /* ******************************************************************************
    
    Stage 0  -  1 cycle
    
    (x^16 - 1) -> (x^8 - 1) and (x^8 + 1) 
    
    c[0:7] = c[0:7] + c[8:15]
    
    c[8:15] = c[0:7] - c[8:15]
    
    */
   
   
   generate
      for (i = 0; i < LEN/2; i++) begin
	 
	 bf1 #(.W(W)) bf_s0 
	       (
		.clk_i, 
		.rst_ni, 
		.en_i, 
		.a_re_i(c_re_s0[i]), 
		.a_im_i(c_im_s0[i]), 
		.b_re_i(c_re_s0[i+LEN/2]), 
		.b_im_i(c_im_s0[i+LEN/2]),
		.a_re_o(c_re_s1[i]), 
		.a_im_o(c_im_s1[i]), 
		.b_re_o(c_re_s1[i+LEN/2]), 
		.b_im_o(c_im_s1[i+LEN/2])
		);
      end
   endgenerate

   
   /* ******************************************************************************
    
    Stage 1  -  1 cycle
    
    (x^8 - 1) -> (x^4 - 1) and (x^4 + 1) 
    (x^8 + 1) -> (x^4 - i) and (x^4 + i)
    
    c[0:3] = c[0:3] + c[4:7]
    
    c[0:3] = c[0:3] - c[4:7]
    
    c[8:11] = c[8:11] + i * c[12:15] 
    
    c[12:15] = c[8:11] - i * c[12:15] 
    
    */
       
   generate
      for (i = 0; i < LEN/4; i++) begin

	 bf1 #(.W(W)) bf_0_s1 
	       (
		.clk_i, 
		.rst_ni, 
		.en_i, 
		.a_re_i(c_re_s1[i]), 
		.a_im_i(c_im_s1[i]), 
		.b_re_i(c_re_s1[i+LEN/4]), 
		.b_im_i(c_im_s1[i+LEN/4]),
		.a_re_o(c_re_s2[i]), 
		.a_im_o(c_im_s2[i]), 
		.b_re_o(c_re_s2[i+LEN/4]), 
		.b_im_o(c_im_s2[i+LEN/4])
		);
	 
	 bf2 #(.W(W)) bf_1_s1 
	   (
	    .clk_i, 
	    .rst_ni, 
	    .en_i, 
	    .a_re_i(c_re_s1[i+LEN/2]), 
	    .a_im_i(c_im_s1[i+LEN/2]), 
	    .b_re_i(c_re_s1[i+LEN/2+LEN/4]), 
	    .b_im_i(c_im_s1[i+LEN/2+LEN/4]),
	    .a_re_o(c_re_s2[i+LEN/2]), 
	    .a_im_o(c_im_s2[i+LEN/2]), 
	    .b_re_o(c_re_s2[i+LEN/2+LEN/4]), 
	    .b_im_o(c_im_s2[i+LEN/2+LEN/4])
	    );

      end 

      
  endgenerate

   

   /* ******************************************************************************
   
    Stage 2  for first half (8 complex) - 2 cycles
    
   
    (x^4 - 1) -> (x^2 - 1) and (x^2 + 1) 
    (x^4 + 1) -> (x^2 - i) and (x^2 + i)
    
    c[0:1] = c[0:1] + c[2:3]
  
    c[2:3] = c[0:1] - c[2:3]
  
    c[4:5] = c[4:5] + i * c[6:7] 
   
    c[6:7] = c[4:5] - i * c[6:7] 
       
  */

   generate
      for (i = 0; i < 2; i++) begin

	 bf1 #(.W(W)) bf_0_s2 
	       (
		.clk_i, 
		.rst_ni, 
		.en_i, 
		.a_re_i(c_re_s2[i]), 
		.a_im_i(c_im_s2[i]), 
		.b_re_i(c_re_s2[i+2]), 
		.b_im_i(c_im_s2[i+2]),
		.a_re_o(c_re_s3[i]), 
		.a_im_o(c_im_s3[i]), 
		.b_re_o(c_re_s3[i+2]), 
		.b_im_o(c_im_s3[i+2])
		);

	 bf2 #(.W(W)) bf_1_s2 
	       (
		.clk_i, 
		.rst_ni, 
		.en_i, 
		.a_re_i(c_re_s2[4+i]), 
		.a_im_i(c_im_s2[4+i]), 
		.b_re_i(c_re_s2[4+i+2]), 
		.b_im_i(c_im_s2[4+i+2]),
		.a_re_o(c_re_s3[4+i]), 
		.a_im_o(c_im_s3[4+i]), 
		.b_re_o(c_re_s3[4+i+2]), 
		.b_im_o(c_im_s3[4+i+2])
		);
	 	 

      end // for (i = 0; i < 2; i++)
   endgenerate

   
   /*  ******************************************************************************

    Stage 3 for first 4 values - 1 cycle
   
    x^2 - 1 =>  x - 1 and x + 1
    x^2 + 1 =>  x - i and x + i
    
    */

   bf1 #(.W(W)) bf_0_s3
     (
      .clk_i, 
      .rst_ni, 
      .en_i, 
      .a_re_i(c_re_s3[0]), 
      .a_im_i(c_im_s3[0]), 
      .b_re_i(c_re_s3[1]), 
      .b_im_i(c_im_s3[1]),
      .a_re_o(c_re_tmp[0]), 
      .a_im_o(c_im_tmp[0]), 
      .b_re_o(c_re_tmp[1]), 
      .b_im_o(c_im_tmp[1])
      );
   
   bf2 #(.W(W)) bf_1_s3 
     (
      .clk_i, 
      .rst_ni, 
      .en_i, 
      .a_re_i(c_re_s3[2]), 
      .a_im_i(c_im_s3[2]), 
      .b_re_i(c_re_s3[3]), 
      .b_im_i(c_im_s3[3]),
      .a_re_o(c_re_tmp[2]), 
      .a_im_o(c_im_tmp[2]), 
      .b_re_o(c_re_tmp[3]), 
      .b_im_o(c_im_tmp[3])
      );


   //save [0:3]
   generate
      for (i = 0; i < 4; i++) begin
	 pipe_reg #(.WIDTH(W), .DEPTH(LAT_SEC_HALF)) pipe_reg0_s4 (.clk_i, .rst_ni, .en_i, .input_i(c_re_tmp[i]), .output_o(c_re_s4[i]));
	 pipe_reg #(.WIDTH(W), .DEPTH(LAT_SEC_HALF)) pipe_reg1_s4 (.clk_i, .rst_ni, .en_i, .input_i(c_im_tmp[i]), .output_o(c_im_s4[i]));
      end	 
   endgenerate


   /* *******************************************************************************
    Stage 3 for c[4:7]  - 2 cycles
    
    c[4:5]
    x^2 - i => x - sqrt(i) and    x + sqrt(i)   
    
    with  sqrt(i)  = (1+i)*2^((W-1)/2)
    
    c[6:7]
    x^2 + i => x - sqrt(-i) and   x + sqrt(-i) 

    with sqrt(-i) = (-1+i)*2^((W-1)/2)
    
    */

   bf3 #(.W(W)) bf_2_s3
     (
      .clk_i, 
      .rst_ni, 
      .en_i, 
      .a_re_i(c_re_s3[4]), 
      .a_im_i(c_im_s3[4]), 
      .b_re_i(c_re_s3[5]), 
      .b_im_i(c_im_s3[5]),
      .a_re_o(c_re_tmp[4]), 
      .a_im_o(c_im_tmp[4]), 
      .b_re_o(c_re_tmp[5]), 
      .b_im_o(c_im_tmp[5])
      );

   bf4 #(.W(W)) bf_3_s3
     (
      .clk_i, 
      .rst_ni, 
      .en_i, 
      .a_re_i(c_re_s3[6]), 
      .a_im_i(c_im_s3[6]), 
      .b_re_i(c_re_s3[7]), 
      .b_im_i(c_im_s3[7]),
      .a_re_o(c_re_tmp[6]), 
      .a_im_o(c_im_tmp[6]), 
      .b_re_o(c_re_tmp[7]), 
      .b_im_o(c_im_tmp[7])
      );

   
  generate
     for (i = 4; i < 8; i++) begin
   	pipe_reg #(.WIDTH(W), .DEPTH(LAT_SEC_HALF-1)) pipe_reg2_s4 (.clk_i, .rst_ni, .en_i, .input_i(c_re_tmp[i]), .output_o(c_re_s4[i]));
	pipe_reg #(.WIDTH(W), .DEPTH(LAT_SEC_HALF-1)) pipe_reg3_s4 (.clk_i, .rst_ni, .en_i, .input_i(c_im_tmp[i]), .output_o(c_im_s4[i]));
     end
  endgenerate

   
   
   /*  *******************************************************************************
    Stage 2a  - second half (8 complex) - 4 cycles
    
    (x^4 - i) * w_16.^[0:3] -> x^4 -1   
    
    with  w_16.^[0:3]       = [(1 0)   (c s)  (1+i)*2^((W-1)/2)   (s c)]
    
    (x^4 + i) * conj(w).^[0:3] -> x^4 -1 
    
    with conj(w^16).^[0:3]  = [(1 0)   (c -s) (1-i)*2^((W-1)/2)  (s -c)]
        
    */

   
   pipe_reg #(.WIDTH(W), .DEPTH(4)) pipe_reg0_s2a (.clk_i, .rst_ni, .en_i, .input_i(c_re_s2[8]), .output_o(c_re_tmp[8]));
   pipe_reg #(.WIDTH(W), .DEPTH(4)) pipe_reg1_s2a (.clk_i, .rst_ni, .en_i, .input_i(c_im_s2[8]), .output_o(c_im_tmp[8]));

   pipe_reg #(.WIDTH(W), .DEPTH(4)) pipe_reg2_s2a (.clk_i, .rst_ni, .en_i, .input_i(c_re_s2[12]), .output_o(c_re_tmp[12]));
   pipe_reg #(.WIDTH(W), .DEPTH(4)) pipe_reg3_s2a (.clk_i, .rst_ni, .en_i, .input_i(c_im_s2[12]), .output_o(c_im_tmp[12]));
   
   
   //tmp = (b_re + i b_im)(1+i) = (b_re + i b_im) + i (b_re + i b_im)   = b_re - b_im + i(b_re + b_im)

   mrsn_add_sub #(.N(W)) adder_1_s2a  (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i({c_re_s2[10][P8:0], c_re_s2[10][W-1:P8+1]}), .b_i({c_im_s2[10][P8:0], c_im_s2[10][W-1:P8+1]}), .sum_o(tmp_re_10));
   mrsn_add_sub #(.N(W)) adder_2_s2a  (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i({c_re_s2[10][P8:0], c_re_s2[10][W-1:P8+1]}), .b_i({c_im_s2[10][P8:0], c_im_s2[10][W-1:P8+1]}), .sum_o(tmp_im_10));

   //tmp = (b_re + i b_im)(1-i) = (b_re + i b_im) - i (b_re + i b_im)   = b_re + b_im + i(-b_re + b_im)

   mrsn_add_sub #(.N(W)) adder_3_s2a  (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i({c_re_s2[14][P8:0], c_re_s2[14][W-1:P8+1]}), .b_i({c_im_s2[14][P8:0], c_im_s2[14][W-1:P8+1]}), .sum_o(tmp_re_14));
   mrsn_add_sub #(.N(W)) adder_4_s2a  (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i({c_im_s2[14][P8:0], c_im_s2[14][W-1:P8+1]}), .b_i({c_re_s2[14][P8:0], c_re_s2[14][W-1:P8+1]}), .sum_o(tmp_im_14));

   pipe_reg #(.WIDTH(W), .DEPTH(3)) pipe_reg4_s2a (.clk_i, .rst_ni, .en_i, .input_i(tmp_re_10), .output_o(c_re_tmp[10]));
   pipe_reg #(.WIDTH(W), .DEPTH(3)) pipe_reg5_s2a (.clk_i, .rst_ni, .en_i, .input_i(tmp_im_10), .output_o(c_im_tmp[10]));

   
   pipe_reg #(.WIDTH(W), .DEPTH(3)) pipe_reg6_s2a (.clk_i, .rst_ni, .en_i, .input_i(tmp_re_14), .output_o(c_re_tmp[14]));
   pipe_reg #(.WIDTH(W), .DEPTH(3)) pipe_reg7_s2a (.clk_i, .rst_ni, .en_i, .input_i(tmp_im_14), .output_o(c_im_tmp[14]));


   
   //(re,im)*(c,s) = (re*c - im*s, re*s + im*c) with (c,s) = (c16 s16) 
   generate
      if (W == 13) begin
	 mrsn_omega16_p13 #(.W(W)) mult_0 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[9]), .x_c_o(tmp0_re_c_s2), .x_s_o(tmp0_re_s_s2));
	 mrsn_omega16_p13 #(.W(W)) mult_1 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[9]), .x_c_o(tmp0_im_c_s2), .x_s_o(tmp0_im_s_s2));
      end else if (W == 19) begin
	 mrsn_omega16_p19 #(.W(W)) mult_0 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[9]), .x_c_o(tmp0_re_c_s2), .x_s_o(tmp0_re_s_s2));
	 mrsn_omega16_p19 #(.W(W)) mult_1 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[9]), .x_c_o(tmp0_im_c_s2), .x_s_o(tmp0_im_s_s2));
      end
   endgenerate
   mrsn_add_sub #(.N(W)) adder_5_s2a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp0_re_c_s2), .b_i(tmp0_im_s_s2), .sum_o(c_re_tmp[9]));
   mrsn_add_sub #(.N(W)) adder_6_s2a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp0_re_s_s2), .b_i(tmp0_im_c_s2), .sum_o(c_im_tmp[9]));

   
   //(re,im)*(s,c) = (re*s - im*c, re*c + im*s)  (c,s) = (c16 s16)
   generate
      if (W == 13) begin
	 mrsn_omega16_p13 #(.W(W)) mult_2 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[11]), .x_c_o(tmp1_re_c_s2), .x_s_o(tmp1_re_s_s2));
	 mrsn_omega16_p13 #(.W(W)) mult_3 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[11]), .x_c_o(tmp1_im_c_s2), .x_s_o(tmp1_im_s_s2));
      end else if (W == 19) begin
	 mrsn_omega16_p19 #(.W(W)) mult_2 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[11]), .x_c_o(tmp1_re_c_s2), .x_s_o(tmp1_re_s_s2));
	 mrsn_omega16_p19 #(.W(W)) mult_3 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[11]), .x_c_o(tmp1_im_c_s2), .x_s_o(tmp1_im_s_s2));
      end
   endgenerate
   mrsn_add_sub #(.N(W)) adder_7_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp1_re_s_s2), .b_i(tmp1_im_c_s2), .sum_o(c_re_tmp[11]));
   mrsn_add_sub #(.N(W)) adder_8_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp1_im_s_s2), .b_i(tmp1_re_c_s2), .sum_o(c_im_tmp[11]));
   

   //(re,im)*(c,-s) = (re*c + im*s, -re*s + im*c) 
   generate
      if (W == 13) begin
	 mrsn_omega16_p13 #(.W(W)) mult_4 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[13]), .x_c_o(tmp2_re_c_s2), .x_s_o(tmp2_re_s_s2));
	 mrsn_omega16_p13 #(.W(W)) mult_5 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[13]), .x_c_o(tmp2_im_c_s2), .x_s_o(tmp2_im_s_s2));
      end else if (W == 19) begin
	 mrsn_omega16_p19 #(.W(W)) mult_4 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[13]), .x_c_o(tmp2_re_c_s2), .x_s_o(tmp2_re_s_s2));
	 mrsn_omega16_p19 #(.W(W)) mult_5 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[13]), .x_c_o(tmp2_im_c_s2), .x_s_o(tmp2_im_s_s2));
      end
   endgenerate
   mrsn_add_sub #(.N(W)) adder_9_s2a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp2_re_c_s2), .b_i(tmp2_im_s_s2), .sum_o(c_re_tmp[13]));
   mrsn_add_sub #(.N(W)) adder_10_s2a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp2_im_c_s2), .b_i(tmp2_re_s_s2), .sum_o(c_im_tmp[13]));

 
   //(re,im)*(s,-c) = (re*s + im*c, -re*c + im*s)  (c,s) = (c16 s16)
   generate
      if (W == 13) begin
	 mrsn_omega16_p13 #(.W(W)) mult_6 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[15]), .x_c_o(tmp3_re_c_s2), .x_s_o(tmp3_re_s_s2));
	 mrsn_omega16_p13 #(.W(W)) mult_7 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[15]), .x_c_o(tmp3_im_c_s2), .x_s_o(tmp3_im_s_s2));
      end else if (W == 19) begin
	 mrsn_omega16_p19 #(.W(W)) mult_6 (.clk_i, .rst_ni, .en_i, .x_i(c_re_s2[15]), .x_c_o(tmp3_re_c_s2), .x_s_o(tmp3_re_s_s2));
	 mrsn_omega16_p19 #(.W(W)) mult_7 (.clk_i, .rst_ni, .en_i, .x_i(c_im_s2[15]), .x_c_o(tmp3_im_c_s2), .x_s_o(tmp3_im_s_s2));
      end
   endgenerate
   mrsn_add_sub #(.N(W)) adder_11_s2a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp3_re_s_s2), .b_i(tmp3_im_c_s2), .sum_o(c_re_tmp[15]));
   mrsn_add_sub #(.N(W)) adder_12_s2a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp3_im_s_s2), .b_i(tmp3_re_c_s2), .sum_o(c_im_tmp[15]));


   /* ******************************************************************************
   
    Stage 2b for second half  - 1 cycle
    
   
    (x^4 - 1) -> (x^2 - 1) and (x^2 + 1)  
    
    twice for c[8:11] and c[12:15]
   
    */

   
   generate
      for (j = 0; j < 2; j++) begin
	 for (i = 0; i < 2; i++) begin
	 
	 bf1 #(.W(W)) bf_2_s2 
	       (
		.clk_i, 
		.rst_ni, 
		.en_i, 
		.a_re_i(c_re_tmp[8+4*j+i]), 
		.a_im_i(c_im_tmp[8+4*j+i]), 
		.b_re_i(c_re_tmp[8+4*j+i+2]), 
		.b_im_i(c_im_tmp[8+4*j+i+2]),
		.a_re_o(c_re_s3[8+4*j+i]), 
		.a_im_o(c_im_s3[8+4*j+i]), 
		.b_re_o(c_re_s3[8+4*j+i+2]), 
		.b_im_o(c_im_s3[8+4*j+i+2])
		);

	 end // for (i = 0; i < 2; i++)
      end // for (j = 0; j < 2; j++)
   endgenerate
   

   /*  ******************************************************************************
    Stage 3 for second half - 1 cycle
    
    
    
    x^2 - 1 =>  x - 1 and x + 1
    x^2 + 1 =>  x - i and x + i
  
    twice for c[8:11] and c[12:15]
    */

   
   
   generate
      for (j = 0; j < 2; j++) begin
	 bf1 #(.W(W)) bf_0_s3
	       (
		.clk_i, 
		.rst_ni, 
		.en_i, 
		.a_re_i(c_re_s3[8+4*j+0]), 
		.a_im_i(c_im_s3[8+4*j+0]), 
		.b_re_i(c_re_s3[8+4*j+1]), 
		.b_im_i(c_im_s3[8+4*j+1]),
		.a_re_o(c_re_s4[8+4*j+0]), 
		.a_im_o(c_im_s4[8+4*j+0]), 
		.b_re_o(c_re_s4[8+4*j+1]), 
		.b_im_o(c_im_s4[8+4*j+1])
		);
	 
	 bf2 #(.W(W)) bf_1_s3 
	   (
	    .clk_i, 
	    .rst_ni, 
	    .en_i, 
	    .a_re_i(c_re_s3[8+4*j+2]), 
	    .a_im_i(c_im_s3[8+4*j+2]), 
	    .b_re_i(c_re_s3[8+4*j+3]), 
	    .b_im_i(c_im_s3[8+4*j+3]),
	    .a_re_o(c_re_s4[8+4*j+2]), 
	    .a_im_o(c_im_s4[8+4*j+2]), 
	    .b_re_o(c_re_s4[8+4*j+3]), 
	    .b_im_o(c_im_s4[8+4*j+3])
	    );
	 
      end 
   endgenerate

endmodule
