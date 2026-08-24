/*

 Inverse NTT16 for 8 complex values with output 16 real values
 
 
 NTT is evaluated modulo Mersenne primes P0 = 2^13-1 and P1 = 2^19-1
 (or P1 = 2^17-1 - Work in progress)

*/


`define ADD  1'b0
`define SUB  1'b1

`define W0   13
`define W1   19

	      
module mrsn_intt_real16 #(
    parameter  WIDTH = 32,
    parameter  LEN   = 16
) (
    input logic			          clk_i             ,
    input logic			          rst_ni            ,
    input logic			          en_i              ,
    input logic  [WIDTH-1:0]  c_i [LEN - 1 : 0] ,
    output logic [WIDTH-1:0]	a_o [LEN - 1 : 0] 
);

   
  localparam W0 = `W0;

  localparam W1 = `W1;

   
  localparam W_NOTUSED = WIDTH - W0 - W1;

  localparam P0_om8 = (W0 - 1) / 2;
  localparam P1_om8 = (W1 - 1) / 2;

   
  logic [W0-1:0] a0_s0   [  LEN];
  logic [W1-1:0] a1_s0   [  LEN];

  logic [W0-1:0] c0_re_s1[LEN/2];
  logic [W0-1:0] c0_im_s1[LEN/2];

  logic [W1-1:0] c1_re_s1[LEN/2];
  logic [W1-1:0] c1_im_s1[LEN/2];

  logic [W0-1:0] c0_re_s1_q[LEN/4], tmp0_re_s1[LEN/4];
  logic [W0-1:0] c0_im_s1_q[LEN/4], tmp0_im_s1[LEN/4];

  logic [W1-1:0] c1_re_s1_q[LEN/4], tmp1_re_s1[LEN/4];
  logic [W1-1:0] c1_im_s1_q[LEN/4], tmp1_im_s1[LEN/4];
   
   logic [W0-1:0] tmp0_re_c_s2a, tmp0_re_s_s2a, tmp0_im_c_s2a, tmp0_im_s_s2a,
		 tmp2_re_c_s2a, tmp2_re_s_s2a, tmp2_im_c_s2a, tmp2_im_s_s2a,
		 tmp4_re_c_s2a, tmp4_re_s_s2a, tmp4_im_c_s2a, tmp4_im_s_s2a,
		 tmp6_re_c_s2a, tmp6_re_s_s2a, tmp6_im_c_s2a, tmp6_im_s_s2a,
		 tmp8_re_c_s2a, tmp8_re_s_s2a, tmp8_im_c_s2a, tmp8_im_s_s2a,
		 tmp10_re_c_s2a, tmp10_re_s_s2a, tmp10_im_c_s2a, tmp10_im_s_s2a;

  logic [W1-1:0] tmp1_re_c_s2a, tmp1_re_s_s2a, tmp1_im_c_s2a, tmp1_im_s_s2a,
		 tmp3_re_c_s2a, tmp3_re_s_s2a, tmp3_im_c_s2a, tmp3_im_s_s2a,
		 tmp5_re_c_s2a, tmp5_re_s_s2a, tmp5_im_c_s2a, tmp5_im_s_s2a,
		 tmp7_re_c_s2a, tmp7_re_s_s2a, tmp7_im_c_s2a, tmp7_im_s_s2a,
		 tmp9_re_c_s2a, tmp9_re_s_s2a, tmp9_im_c_s2a, tmp9_im_s_s2a,
		 tmp11_re_c_s2a, tmp11_re_s_s2a, tmp11_im_c_s2a, tmp11_im_s_s2a;

   
  logic [W0-1:0] c0_re_s2a[LEN/2];
  logic [W0-1:0] c0_im_s2a[LEN/2];

  logic [W1-1:0] c1_re_s2a[LEN/2];
  logic [W1-1:0] c1_im_s2a[LEN/2];

  logic [W0-1:0] c0_re_s2[LEN/2];
  logic [W0-1:0] c0_im_s2[LEN/2];

  logic [W1-1:0] c1_re_s2[LEN/2];
  logic [W1-1:0] c1_im_s2[LEN/2];
   
  logic [W0-1:0] c0_re_s3[LEN/2];
  logic [W0-1:0] c0_im_s3[LEN/2];

  logic [W1-1:0] c1_re_s3[LEN/2];
  logic [W1-1:0] c1_im_s3[LEN/2];

  logic [W0-1:0] c0_re_s4[LEN/2];
  logic [W0-1:0] c0_im_s4[LEN/2];

  logic [W1-1:0] c1_re_s4[LEN/2];
  logic [W1-1:0] c1_im_s4[LEN/2];


   logic [W0-1:0] cucu;
   

   
  genvar i,j;

   generate
      for (i=0; i<LEN/2; i++) begin
	 
	 assign c0_re_s4[i] = c_i[2*i]  [W0-1:0];
	 assign c0_im_s4[i] = c_i[2*i+1][W0-1:0];
	 
	 assign c1_re_s4[i] = c_i[2*i]  [W1+W0-1:W0];
	 assign c1_im_s4[i] = c_i[2*i+1][W1+W0-1:W0];
      end
   endgenerate
      
 /* ******************************************************************************
   
   
   Stage 3  - 1 cycle 

  direct:   
    c[0] = c[0] + c[1]
    c[1] = c[0] - c[1]
    
    c[2] = c[2] + i c[3]
    c[3] = c[2] - i c[3]
    
    c[4] = c[4] + c[5]
    c[5] = c[4] - c[5]
    
    c[6] = c[6] + i c[7]
    c[7] = c[6] - i c[7]
  
  inverse:
  
  c[0] = c[0] + c[1]
  c[1] = c[0] - c[1]
  
  c[2] = c[2] + c[3]
  c[3] = (c[2] - c[3])/i = (c[3] - c[2])*i   = (c_re[3] - c_re[2] + i (c_im[3] - c_im[2]))*i 
       = c_im[2] - c_im[3] + i (c_re[3] - c_re[2])
  
  same for 4,5,6,7

      */

   
   generate
     for (j = 0; j<2; j++) begin

	mrsn_add_sub #(.N(W0)) adder0_0_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_re_s4[0+4*j]), .b_i(c0_re_s4[1+4*j]), .sum_o(c0_re_s3[0+4*j]));
	mrsn_add_sub #(.N(W0)) adder0_1_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_re_s4[0+4*j]), .b_i(c0_re_s4[1+4*j]), .sum_o(c0_re_s3[1+4*j]));
	
	mrsn_add_sub #(.N(W0)) adder0_2_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_im_s4[0+4*j]), .b_i(c0_im_s4[1+4*j]), .sum_o(c0_im_s3[0+4*j]));
	mrsn_add_sub #(.N(W0)) adder0_3_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_im_s4[0+4*j]), .b_i(c0_im_s4[1+4*j]), .sum_o(c0_im_s3[1+4*j]));
	
	mrsn_add_sub #(.N(W0)) adder0_4_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_re_s4[2+4*j]), .b_i(c0_re_s4[3+4*j]), .sum_o(c0_re_s3[2+4*j]));
	mrsn_add_sub #(.N(W0)) adder0_5_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_im_s4[2+4*j]), .b_i(c0_im_s4[3+4*j]), .sum_o(c0_re_s3[3+4*j]));

	mrsn_add_sub #(.N(W0)) adder0_6_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_im_s4[2+4*j]), .b_i(c0_im_s4[3+4*j]), .sum_o(c0_im_s3[2+4*j]));
	mrsn_add_sub #(.N(W0)) adder0_7_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_re_s4[3+4*j]), .b_i(c0_re_s4[2+4*j]), .sum_o(c0_im_s3[3+4*j]));


	
	mrsn_add_sub #(.N(W1)) adder1_0_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_re_s4[0+4*j]), .b_i(c1_re_s4[1+4*j]), .sum_o(c1_re_s3[0+4*j]));
	mrsn_add_sub #(.N(W1)) adder1_1_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_re_s4[0+4*j]), .b_i(c1_re_s4[1+4*j]), .sum_o(c1_re_s3[1+4*j]));
	
	mrsn_add_sub #(.N(W1)) adder1_2_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_im_s4[0+4*j]), .b_i(c1_im_s4[1+4*j]), .sum_o(c1_im_s3[0+4*j]));
	mrsn_add_sub #(.N(W1)) adder1_3_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_im_s4[0+4*j]), .b_i(c1_im_s4[1+4*j]), .sum_o(c1_im_s3[1+4*j]));

	mrsn_add_sub #(.N(W1)) adder1_4_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_re_s4[2+4*j]), .b_i(c1_re_s4[3+4*j]), .sum_o(c1_re_s3[2+4*j]));
	mrsn_add_sub #(.N(W1)) adder1_5_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_im_s4[2+4*j]), .b_i(c1_im_s4[3+4*j]), .sum_o(c1_re_s3[3+4*j]));

	mrsn_add_sub #(.N(W1)) adder1_6_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_im_s4[2+4*j]), .b_i(c1_im_s4[3+4*j]), .sum_o(c1_im_s3[2+4*j]));
	mrsn_add_sub #(.N(W1)) adder1_7_s3 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_re_s4[3+4*j]), .b_i(c1_re_s4[2+4*j]), .sum_o(c1_im_s3[3+4*j]));
  
     end // for (j = 0; j<2; j++)
  endgenerate

   /* ******************************************************************************
   
   Stage 2  - 1 cycle 

   For c[0:3] and c[4:7] do
   
      (x^4-1)_complex => 
            (x^2-1)_complex and 
            (x^2+1)_complex

   c[0] = c[0] + c[2]
   c[2] = c[0] - c[2]
   
   c[1] = c[1] + c[3]
   c[3] = c[1] - c[3]
   
   c[4] = c[4] + c[6]
   c[6] = c[4] - c[6]
   
   c[5] = c[5] + c[7]
   c[7] = c[5] - c[7]
   
   For the last 2, correct sign for the imaginary part of c[7]

  */
   
  generate
     for (j = 0; j< 2; j++) begin
       for (i = 0; i < 2; i++) begin
	  mrsn_add_sub #(.N(W0)) adder0_0_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_re_s3[i+4*j]), .b_i(c0_re_s3[i+4*j+2]), .sum_o(c0_re_s2[i+4*j]));
	  mrsn_add_sub #(.N(W0)) adder0_1_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_re_s3[i+4*j]), .b_i(c0_re_s3[i+4*j+2]), .sum_o(c0_re_s2[i+4*j+2]));
	  
	  mrsn_add_sub #(.N(W1)) adder1_0_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_re_s3[i+4*j]), .b_i(c1_re_s3[i+4*j+2]), .sum_o(c1_re_s2[i+4*j]));
	  mrsn_add_sub #(.N(W1)) adder1_1_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_re_s3[i+4*j]), .b_i(c1_re_s3[i+4*j+2]), .sum_o(c1_re_s2[i+4*j+2]));
       end
     end
  endgenerate

   generate
      for (i = 0; i < 2; i++) begin
	mrsn_add_sub #(.N(W0)) adder0_2_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_im_s3[i]), .b_i(c0_im_s3[i+2]), .sum_o(c0_im_s2[i]));
	mrsn_add_sub #(.N(W0)) adder0_3_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_im_s3[i]), .b_i(c0_im_s3[i+2]), .sum_o(c0_im_s2[i+2]));
	 
        mrsn_add_sub #(.N(W1)) adder1_2_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_im_s3[i]), .b_i(c1_im_s3[i+2]), .sum_o(c1_im_s2[i]));
	mrsn_add_sub #(.N(W1)) adder1_3_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_im_s3[i]), .b_i(c1_im_s3[i+2]), .sum_o(c1_im_s2[i+2]));
      end
   endgenerate


   mrsn_add_sub #(.N(W0)) adder0_4_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_im_s3[4]), .b_i(c0_im_s3[6]), .sum_o(c0_im_s2[4]));
   mrsn_add_sub #(.N(W0)) adder0_5_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_im_s3[4]), .b_i(c0_im_s3[6]), .sum_o(c0_im_s2[6]));
   
   mrsn_add_sub #(.N(W1)) adder1_4_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_im_s3[4]), .b_i(c1_im_s3[6]), .sum_o(c1_im_s2[4]));
   mrsn_add_sub #(.N(W1)) adder1_5_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_im_s3[4]), .b_i(c1_im_s3[6]), .sum_o(c1_im_s2[6]));

   /* direct with sign correction for imaginary part: 
    c_im[5] = c_im[5] - c_im[7]
    c_im[7] = c_im[5] + c_im[7]
 
    inverse with sign correction for c_im[7]
    c_im[5] = c_im[5] + c_im[7]
    c_im[7] = c_im[7] - c_im[5] 
   */
 
   mrsn_add_sub #(.N(W0)) adder0_6_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_im_s3[5]), .b_i(c0_im_s3[7]), .sum_o(c0_im_s2[5]));
   mrsn_add_sub #(.N(W0)) adder0_7_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_im_s3[7]), .b_i(c0_im_s3[5]), .sum_o(c0_im_s2[7]));
   
   mrsn_add_sub #(.N(W1)) adder1_6_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_im_s3[5]), .b_i(c1_im_s3[7]), .sum_o(c1_im_s2[5]));
   mrsn_add_sub #(.N(W1)) adder1_7_s2 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_im_s3[7]), .b_i(c1_im_s3[5]), .sum_o(c1_im_s2[7]));

 /* ******************************************************************************
    
  Stage 1b - 4 cycles where 3 cycles used to multiply by c and s using
  only adds and rotations and 1 cycle for final add/sub
    
  Pointwise multiplication 
        
  c[0:3] .* w32 .^ [0:3]
  c[4:7] .* conj(w32^3) .^ [0:3]
  
  w32 .^ [0:3]              = [(1,0)  w32   w16  w32^3]                     = [(1 0), (c32 s32)      (c16  s16) (c32_3 s32_3)] 

  (conj(w32// )^3) .^ [0:3] = [(1,0)  conj(w32)^3 conj(w16)^3  conj(w32)^9] = [(1 0), (c32_3 -s32_3) (s16 -c16) (-s32 -c32)] 
   
 */
   
   //c[0:3] .* [(1 0)  (c32 s32)  (c16 s16) (c32_3 s32_3)] 
   
   //(1 0)
   
   pipe_reg #(.WIDTH(W0), .DEPTH(3+1)) pipe_reg0_s1 (.clk_i, .rst_ni, .en_i, .input_i(c0_re_s2[0]), .output_o(c0_re_s2a[0]));
   
   pipe_reg #(.WIDTH(W0), .DEPTH(3+1)) pipe_reg1_s1 (.clk_i, .rst_ni, .en_i, .input_i(c0_im_s2[0]), .output_o(c0_im_s2a[0]));

   pipe_reg #(.WIDTH(W1), .DEPTH(3+1)) pipe_reg4_s1 (.clk_i, .rst_ni, .en_i, .input_i(c1_re_s2[0]), .output_o(c1_re_s2a[0]));
   
   pipe_reg #(.WIDTH(W1), .DEPTH(3+1)) pipe_reg5_s1 (.clk_i, .rst_ni, .en_i, .input_i(c1_im_s2[0]), .output_o(c1_im_s2a[0]));

   //direct:   (re,im)*(c,s) = (re*c - im*s, re*s + im*c) with (c,s) = (c32  s32)

   //inverse: (re,im)*(c,-s) = (re*c + im*s, -re*s + im*c) with (c,s) = (c32  s32)
   
   mrsn_omega32_p13 #(.W(W0)) mult0_0 (.clk_i, .rst_ni, .en_i, .x_i(c0_re_s2[1]), .x_c_o(tmp0_re_c_s2a), .x_s_o(tmp0_re_s_s2a));
   mrsn_omega32_p13 #(.W(W0)) mult0_1 (.clk_i, .rst_ni, .en_i, .x_i(c0_im_s2[1]), .x_c_o(tmp0_im_c_s2a), .x_s_o(tmp0_im_s_s2a));

   mrsn_add_sub #(.N(W0)) adder0_0_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp0_re_c_s2a), .b_i(tmp0_im_s_s2a), .sum_o(c0_re_s2a[1]));
   mrsn_add_sub #(.N(W0)) adder0_1_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp0_im_c_s2a), .b_i(tmp0_re_s_s2a), .sum_o(c0_im_s2a[1]));

   mrsn_omega32_p19 #(.W(W1)) mult1_0 (.clk_i, .rst_ni, .en_i, .x_i(c1_re_s2[1]), .x_c_o(tmp1_re_c_s2a), .x_s_o(tmp1_re_s_s2a));
   mrsn_omega32_p19 #(.W(W1)) mult1_1 (.clk_i, .rst_ni, .en_i, .x_i(c1_im_s2[1]), .x_c_o(tmp1_im_c_s2a), .x_s_o(tmp1_im_s_s2a));

   mrsn_add_sub #(.N(W1)) adder1_0_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp1_re_c_s2a), .b_i(tmp1_im_s_s2a), .sum_o(c1_re_s2a[1]));
   mrsn_add_sub #(.N(W1)) adder1_1_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp1_im_c_s2a), .b_i(tmp1_re_s_s2a), .sum_o(c1_im_s2a[1]));

   //direct: (re,im)*(c,s) = (re*c - im*s, re*s + im*c) with (c,s) = (c16  s16)

   //inverse: (re,im)*(c,-s) = (re*c + im*s, -re*s + im*c) with (c,s) = (c16  s16) 
  
   mrsn_omega16_p13 #(.W(W0)) mult0_2 (.clk_i, .rst_ni, .en_i, .x_i(c0_re_s2[2]), .x_c_o(tmp2_re_c_s2a), .x_s_o(tmp2_re_s_s2a));
   mrsn_omega16_p13 #(.W(W0)) mult0_3 (.clk_i, .rst_ni, .en_i, .x_i(c0_im_s2[2]), .x_c_o(tmp2_im_c_s2a), .x_s_o(tmp2_im_s_s2a));

   mrsn_add_sub #(.N(W0)) adder0_2_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp2_re_c_s2a), .b_i(tmp2_im_s_s2a), .sum_o(c0_re_s2a[2]));
   mrsn_add_sub #(.N(W0)) adder0_3_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp2_im_c_s2a), .b_i(tmp2_re_s_s2a), .sum_o(c0_im_s2a[2]));

   mrsn_omega16_p19 #(.W(W1)) mult1_2 (.clk_i, .rst_ni, .en_i, .x_i(c1_re_s2[2]), .x_c_o(tmp3_re_c_s2a), .x_s_o(tmp3_re_s_s2a));
   mrsn_omega16_p19 #(.W(W1)) mult1_3 (.clk_i, .rst_ni, .en_i, .x_i(c1_im_s2[2]), .x_c_o(tmp3_im_c_s2a), .x_s_o(tmp3_im_s_s2a));

   mrsn_add_sub #(.N(W1)) adder1_2_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp3_re_c_s2a), .b_i(tmp3_im_s_s2a), .sum_o(c1_re_s2a[2]));
   mrsn_add_sub #(.N(W1)) adder1_3_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp3_im_c_s2a), .b_i(tmp3_re_s_s2a), .sum_o(c1_im_s2a[2]));

   
   //direct: (re,im)*(c,s) = (re*c - im*s, re*s + im*c)    (c,s) = (c32,   s32)^3 = (c32_3,   s32_3)

   //inverse: (re,im)*(c,s) = (re*c + im*s, -re*s + im*c)    (c,s) = (c32,   s32)^3 = (c32_3,   s32_3)
  
   mrsn_omega32_3_p13 #(.W(W0)) mult0_4 (.clk_i, .rst_ni, .en_i, .x_i(c0_re_s2[3]), .x_c_o(tmp4_re_c_s2a), .x_s_o(tmp4_re_s_s2a));
   mrsn_omega32_3_p13 #(.W(W0)) mult0_5 (.clk_i, .rst_ni, .en_i, .x_i(c0_im_s2[3]), .x_c_o(tmp4_im_c_s2a), .x_s_o(tmp4_im_s_s2a));

   mrsn_add_sub #(.N(W0)) adder0_4_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp4_re_c_s2a), .b_i(tmp4_im_s_s2a), .sum_o(c0_re_s2a[3]));
   mrsn_add_sub #(.N(W0)) adder0_5_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp4_im_c_s2a), .b_i(tmp4_re_s_s2a), .sum_o(c0_im_s2a[3]));

   mrsn_omega32_3_p19 #(.W(W1)) mult1_4 (.clk_i, .rst_ni, .en_i, .x_i(c1_re_s2[3]), .x_c_o(tmp5_re_c_s2a), .x_s_o(tmp5_re_s_s2a));
   mrsn_omega32_3_p19 #(.W(W1)) mult1_5 (.clk_i, .rst_ni, .en_i, .x_i(c1_im_s2[3]), .x_c_o(tmp5_im_c_s2a), .x_s_o(tmp5_im_s_s2a));

   mrsn_add_sub #(.N(W1)) adder1_4_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp5_re_c_s2a), .b_i(tmp5_im_s_s2a), .sum_o(c1_re_s2a[3]));
   mrsn_add_sub #(.N(W1)) adder1_5_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp5_im_c_s2a), .b_i(tmp5_re_s_s2a), .sum_o(c1_im_s2a[3]));

   //c[4:7] .* [(1 0)  (c32_3 -s32_3)  (s16 -c16)  (-s32   -c32) ]
   
   //(1 0)
   
   pipe_reg #(.WIDTH(W0), .DEPTH(3+1)) pipe_reg2_s1 (.clk_i, .rst_ni, .en_i, .input_i(c0_re_s2[4]), .output_o(c0_re_s2a[4]));
   
   pipe_reg #(.WIDTH(W0), .DEPTH(3+1)) pipe_reg3_s1 (.clk_i, .rst_ni, .en_i, .input_i(c0_im_s2[4]), .output_o(c0_im_s2a[4]));

   pipe_reg #(.WIDTH(W1), .DEPTH(3+1)) pipe_reg6_s1 (.clk_i, .rst_ni, .en_i, .input_i(c1_re_s2[4]), .output_o(c1_re_s2a[4]));

   pipe_reg #(.WIDTH(W1), .DEPTH(3+1)) pipe_reg7_s1 (.clk_i, .rst_ni, .en_i, .input_i(c1_im_s2[4]), .output_o(c1_im_s2a[4]));
   
   //direct: (re,im)*(c,-s) = (re*c + im*s, -re*s + im*c)    (c,s) = (c32,   s32)^3 = (c32_3,   s32_3)

   //inverse: (re,im)*(c,s) = (re*c - im*s, re*s + im*c)     (c,s) = (c32,   s32)^3 = (c32_3,   s32_3)
  
   mrsn_omega32_3_p13 #(.W(W0)) mult0_6 (.clk_i, .rst_ni, .en_i, .x_i(c0_re_s2[5]), .x_c_o(tmp6_re_c_s2a), .x_s_o(tmp6_re_s_s2a));
   mrsn_omega32_3_p13 #(.W(W0)) mult0_7 (.clk_i, .rst_ni, .en_i, .x_i(c0_im_s2[5]), .x_c_o(tmp6_im_c_s2a), .x_s_o(tmp6_im_s_s2a));

   mrsn_add_sub #(.N(W0)) adder0_6_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp6_re_c_s2a), .b_i(tmp6_im_s_s2a), .sum_o(c0_re_s2a[5]));
   mrsn_add_sub #(.N(W0)) adder0_7_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp6_im_c_s2a), .b_i(tmp6_re_s_s2a), .sum_o(c0_im_s2a[5]));

   mrsn_omega32_3_p19 #(.W(W1)) mult1_6 (.clk_i, .rst_ni, .en_i, .x_i(c1_re_s2[5]), .x_c_o(tmp7_re_c_s2a), .x_s_o(tmp7_re_s_s2a));
   mrsn_omega32_3_p19 #(.W(W1)) mult1_7 (.clk_i, .rst_ni, .en_i, .x_i(c1_im_s2[5]), .x_c_o(tmp7_im_c_s2a), .x_s_o(tmp7_im_s_s2a));

   mrsn_add_sub #(.N(W1)) adder1_6_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp7_re_c_s2a), .b_i(tmp7_im_s_s2a), .sum_o(c1_re_s2a[5]));
   mrsn_add_sub #(.N(W1)) adder1_7_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp7_im_c_s2a), .b_i(tmp7_re_s_s2a), .sum_o(c1_im_s2a[5]));

   //direct: (re,im)*(s,-c) = (re*s + im*c, -re*c + im*s)  (c,s) = (c16,s16)

   //inverse: (re,im)*(s,c) = (re*s - im*c, re*c + im*s)  (c,s) = (c16,s16)

   mrsn_omega16_p13 #(.W(W0)) mult0_8 (.clk_i, .rst_ni, .en_i, .x_i(c0_re_s2[6]), .x_c_o(tmp8_re_c_s2a), .x_s_o(tmp8_re_s_s2a));
   mrsn_omega16_p13 #(.W(W0)) mult0_9 (.clk_i, .rst_ni, .en_i, .x_i(c0_im_s2[6]), .x_c_o(tmp8_im_c_s2a), .x_s_o(tmp8_im_s_s2a));

   mrsn_add_sub #(.N(W0)) adder0_8_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp8_re_s_s2a), .b_i(tmp8_im_c_s2a), .sum_o(c0_re_s2a[6]));
   mrsn_add_sub #(.N(W0)) adder0_9_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp8_im_s_s2a), .b_i(tmp8_re_c_s2a), .sum_o(c0_im_s2a[6]));

   mrsn_omega16_p19 #(.W(W1)) mult1_8 (.clk_i, .rst_ni, .en_i, .x_i(c1_re_s2[6]), .x_c_o(tmp9_re_c_s2a), .x_s_o(tmp9_re_s_s2a));
   mrsn_omega16_p19 #(.W(W1)) mult1_9 (.clk_i, .rst_ni, .en_i, .x_i(c1_im_s2[6]), .x_c_o(tmp9_im_c_s2a), .x_s_o(tmp9_im_s_s2a));

   mrsn_add_sub #(.N(W1)) adder1_8_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp9_re_s_s2a), .b_i(tmp9_im_c_s2a), .sum_o(c1_re_s2a[6]));
   mrsn_add_sub #(.N(W1)) adder1_9_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp9_im_s_s2a), .b_i(tmp9_re_c_s2a), .sum_o(c1_im_s2a[6]));

   //direct:  (re,im)*(-s,-c) = (-re*s + im*c, - re*c - im*s) with (c,s) = (c32,   s32)

   //inverse: (re,im)*(-s,c) = (-re*s - im*c,   re*c - im*s) with (c,s) = (c32,   s32) 
   
   // NOTE direct: for the imaginary part -re*c -im*s we calculate re*c+im*s and do a sign correction at next stage

   // NOTE inverse: the imaginary part s is sign corrected to -s so that we need (re,im)*(s,c) = (re*s - im*c,   re*c + im*s)

   mrsn_omega32_p13 #(.W(W0)) mult0_10 (.clk_i, .rst_ni, .en_i, .x_i(c0_re_s2[7]), .x_c_o(tmp10_re_c_s2a), .x_s_o(tmp10_re_s_s2a));
   mrsn_omega32_p13 #(.W(W0)) mult0_11 (.clk_i, .rst_ni, .en_i, .x_i(c0_im_s2[7]), .x_c_o(tmp10_im_c_s2a), .x_s_o(tmp10_im_s_s2a));

   mrsn_add_sub #(.N(W0)) adder0_10_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp10_im_c_s2a), .b_i(tmp10_re_s_s2a), .sum_o(c0_re_s2a[7]));
   mrsn_add_sub #(.N(W0)) adder0_11_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp10_re_c_s2a), .b_i(tmp10_im_s_s2a), .sum_o(c0_im_s2a[7]));
   
   mrsn_omega32_p19 #(.W(W1)) mult1_10 (.clk_i, .rst_ni, .en_i, .x_i(c1_re_s2[7]), .x_c_o(tmp11_re_c_s2a), .x_s_o(tmp11_re_s_s2a));
   mrsn_omega32_p19 #(.W(W1)) mult1_11 (.clk_i, .rst_ni, .en_i, .x_i(c1_im_s2[7]), .x_c_o(tmp11_im_c_s2a), .x_s_o(tmp11_im_s_s2a));

   mrsn_add_sub #(.N(W1)) adder1_10_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(tmp11_im_c_s2a), .b_i(tmp11_re_s_s2a), .sum_o(c1_re_s2a[7]));
   mrsn_add_sub #(.N(W1)) adder1_11_s1 (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(tmp11_re_c_s2a), .b_i(tmp11_im_s_s2a), .sum_o(c1_im_s2a[7]));


   /*
    Stage 1a - 2 cycles
    
    Calculate (x^4 - sqrt(i)) and (x^4 + sqrt(i))
    
    c[0:3] = c[0:3] + c[4:7] * sqrt(i)
    c[4:7] = c[0:3] - c[4:7] * sqrt(i)
    
    Here sqrt(i) = (1+i)*2^((W-1)/2) for p = 2^W-1
                 = (1+i)*2^6         for p = 2^13-1
                  = (1+i)*2^9         for p = 2^19-1
    

    Inverse:
    
    c[0:3] = c[0:3] + c[4:7]
    
    c[4:7] = (c[0:3] - c[4:7]) / sqrt(i) =  (c[0:3] + c[4:7])*(1-i)*2^((W-1)/2)
    
    (re + i*im) * (1 - i) = re + im  + i (im - re)  
    */
  generate
    for (i = 0; i < LEN/4; i++) begin

      mrsn_add_sub #(.N(W0)) adder0_0_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_re_s2a[i]), .b_i(c0_re_s2a[i+LEN/4]), .sum_o(c0_re_s1_q[i]));
       
      mrsn_add_sub #(.N(W0)) adder0_1_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c0_im_s2a[i]), .b_i(c0_im_s2a[i+LEN/4]), .sum_o(c0_im_s1_q[i]));

      mrsn_add_sub #(.N(W1)) adder1_0_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_re_s2a[i]), .b_i(c1_re_s2a[i+LEN/4]), .sum_o(c1_re_s1_q[i]));
       
      mrsn_add_sub #(.N(W1)) adder1_1_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(c1_im_s2a[i]), .b_i(c1_im_s2a[i+LEN/4]), .sum_o(c1_im_s1_q[i]));
       

      mrsn_add_sub #(.N(W0)) adder0_2_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_re_s2a[i]), .b_i(c0_re_s2a[i+LEN/4]), .sum_o(tmp0_re_s1[i]));
       
      mrsn_add_sub #(.N(W0)) adder0_3_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c0_im_s2a[i]), .b_i(c0_im_s2a[i+LEN/4]), .sum_o(tmp0_im_s1[i]));

      mrsn_add_sub #(.N(W1)) adder1_2_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_re_s2a[i]), .b_i(c1_re_s2a[i+LEN/4]), .sum_o(tmp1_re_s1[i]));
       
      mrsn_add_sub #(.N(W1)) adder1_3_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(c1_im_s2a[i]), .b_i(c1_im_s2a[i+LEN/4]), .sum_o(tmp1_im_s1[i]));

       

      pipe_reg #(.WIDTH(W0), .DEPTH(1)) pipe_reg0_s1a (.clk_i, .rst_ni, .en_i, .input_i(c0_re_s1_q[i]), .output_o(c0_re_s1[i]));

      pipe_reg #(.WIDTH(W0), .DEPTH(1)) pipe_reg1_s1a (.clk_i, .rst_ni, .en_i, .input_i(c0_im_s1_q[i]), .output_o(c0_im_s1[i]));

      pipe_reg #(.WIDTH(W1), .DEPTH(1)) pipe_reg2_s1a (.clk_i, .rst_ni, .en_i, .input_i(c1_re_s1_q[i]), .output_o(c1_re_s1[i]));

      pipe_reg #(.WIDTH(W1), .DEPTH(1)) pipe_reg3_s1a (.clk_i, .rst_ni, .en_i, .input_i(c1_im_s1_q[i]), .output_o(c1_im_s1[i]));
       
 
      mrsn_add_sub #(.N(W0)) adder0_4_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i({tmp0_re_s1[i][P0_om8:0], tmp0_re_s1[i][W0-1:P0_om8+1]}), .b_i({tmp0_im_s1[i][P0_om8:0], tmp0_im_s1[i][W0-1:P0_om8+1]}), .sum_o(c0_re_s1[i+LEN/4]));

      mrsn_add_sub #(.N(W0)) adder0_5_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i({tmp0_im_s1[i][P0_om8:0], tmp0_im_s1[i][W0-1:P0_om8+1]}), .b_i({tmp0_re_s1[i][P0_om8:0], tmp0_re_s1[i][W0-1:P0_om8+1]}), .sum_o(c0_im_s1[i+LEN/4]));
  
      mrsn_add_sub #(.N(W1)) adder1_4_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i({tmp1_re_s1[i][P1_om8:0], tmp1_re_s1[i][W1-1:P1_om8+1]}), .b_i({tmp1_im_s1[i][P1_om8:0], tmp1_im_s1[i][W1-1:P1_om8+1]}), .sum_o(c1_re_s1[i+LEN/4]));

      mrsn_add_sub #(.N(W1)) adder1_5_s1a (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i({tmp1_im_s1[i][P1_om8:0], tmp1_im_s1[i][W1-1:P1_om8+1]}), .b_i({tmp1_re_s1[i][P1_om8:0], tmp1_re_s1[i][W1-1:P1_om8+1]}), .sum_o(c1_im_s1[i+LEN/4]));
    
    end
  endgenerate

   
 /* ******************************************************************************
   
  Stage 0  -  0 cycles
  
 */

 generate
   for (i = 0; i < LEN / 2; i++) begin
       assign a0_s0[i]       = c0_re_s1[i];
       assign a0_s0[i+LEN/2] = c0_im_s1[i];
       assign a1_s0[i]       = c1_re_s1[i];
       assign a1_s0[i+LEN/2] = c1_im_s1[i];
   end
 endgenerate

 /* ******************************************************************************
  
  Calculate remainder mod (p_13 * p_19) in Z using Chinese Remainder Theorem
  
  
  r =  mod(r0 * t0, p0) * p1 +  mod(r1 * t1, p1) * p0
   
  where t1 = 8321, t0 = 8061.
  
  Multiplication by constants t0 and t1 - 2 cycles
  
  Product with p0 and p1 - 1 cycle
  
  Sum                    - 1 cycle - 33 bits
  
  Correction 1 - 1 cycle  - 32 bits
  
  Correction 2 - 1 cycle  - 32 bits
 
   
  if (r >= p0*p1) r = r - p0*p1
  
  if (r > (p0*p1-1)/2 r = r - p0*p1
    
 */

   
   
   // register output 
  generate
    for (i = 0; i < LEN; i++) begin
       always_ff @(posedge clk_i or negedge rst_ni) begin
	  if (~rst_ni) begin
	     a_o[i]   <= '0;
	  end else begin
	     if (en_i) begin
		a_o[i]  <= {{W_NOTUSED{1'b0}}, a1_s0[i], a0_s0[i]};
	     end
	  end
	  
       end
    end
  endgenerate

endmodule


