/* Modular butterfly type 3 for Mersenne prime modulus p = 2^W-1
 
 
 Input: complex
 
    a
    b
 
 Output:  
    
    a + b * sqrt(i)
    a - b * sqrt(i)
 
    Here   sqrt(i) = (1+i)*2^((W-1)/2)
 
    and multiplication with power of 2^alpha is equivalent to bit
    rotation with alpha positions to the left
 
 */

`define ADD 1'b0
`define SUB 1'b1


module bf3 #(
    parameter W = 13
) (
    input logic clk_i,
    input logic rst_ni,
    input logic en_i,

    input logic [W-1:0] a_re_i,
    input logic [W-1:0] a_im_i,

    input logic [W-1:0] b_re_i,
    input logic [W-1:0] b_im_i,

    output logic [W-1:0] a_re_o,
    output logic [W-1:0] a_im_o,

    output logic [W-1:0] b_re_o,
    output logic [W-1:0] b_im_o

);

  localparam P8 = (W - 1) / 2;

  //  a = a + b * sqrt(i)
  //  b = a - b * sqrt(i)

  //with sqrt(i)  = (1+i)*2^((W-1)/2)

  //tmp = (b_re + i b_im)(1+i) = (b_re + i b_im) + i (b_re + i b_im)   = b_re - b_im + i(b_re + b_im)



  //  a_re + rot(tmp_re) + i (a_im + rot(tmp_im))
  //  a_re - rot(tmp_re) + i (a_im - rot(tmp_im))

  logic [W-1:0] a_re_q, a_im_q, tmp_re, tmp_im;


  pipe_reg #(
      .WIDTH(W),
      .DEPTH(1)
  ) pipe0 (
      .clk_i,
      .rst_ni,
      .en_i,
      .input_i (a_re_i),
      .output_o(a_re_q)
  );
  pipe_reg #(
      .WIDTH(W),
      .DEPTH(1)
  ) pipe1 (
      .clk_i,
      .rst_ni,
      .en_i,
      .input_i (a_im_i),
      .output_o(a_im_q)
  );


  mrsn_add_sub #(
      .N(W)
  ) adder1 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(b_re_i),
      .b_i(b_im_i),
      .sum_o(tmp_re)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder2 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(b_re_i),
      .b_i(b_im_i),
      .sum_o(tmp_im)
  );

  mrsn_add_sub #(
      .N(W)
  ) adder3 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(a_re_q),
      .b_i({tmp_re[P8:0], tmp_re[W-1:P8+1]}),
      .sum_o(a_re_o)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder4 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(a_re_q),
      .b_i({tmp_re[P8:0], tmp_re[W-1:P8+1]}),
      .sum_o(b_re_o)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder5 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(a_im_q),
      .b_i({tmp_im[P8:0], tmp_im[W-1:P8+1]}),
      .sum_o(a_im_o)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder6 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(a_im_q),
      .b_i({tmp_im[P8:0], tmp_im[W-1:P8+1]}),
      .sum_o(b_im_o)
  );


endmodule

