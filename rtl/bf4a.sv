/* Modular butterfly inverse of type 4 for Mersenne prime modulus p = 2^W-1
 
 
 Input: complex
 
    a
    b
 
 Output: complex  
 
     a + b
    (a - b) * conj(sqrt(-i))
    
 
  Here sqrt(-i) = (-1+i)*2^((W-1)/2)
 
  and multiplication with power of 2^alpha is equivalent to bit
  rotation with alpha positions to the left
 
 */

`define ADD 1'b0
`define SUB 1'b1

module bf4a #(
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

  /*
   
   Direct:
   
    a = a + b * sqrt(-i)
    b = a - b * sqrt(-i)
    
    with sqrt(i)  = (-1+i)*2^((W-1)/2)
    
    tmp = (b_re + i b_im)(-1+i) 
        = (-b_re - i b_im) + i (b_re + i b_im)   
     
    tmp = -b_re - b_im + i(b_re - b_im)
        
    a = a_re + rot(tmp_re) + i (a_im + rot(tmp_im))
    b = a_re - rot(tmp_re) + i (a_im - rot(tmp_im))
   
   Inverse:
   
    a = a + b
  
    b*sqrt(-i) = a - b => b = conj(sqrt(-i))*(a-b) = (a-b)*(-1-i)*2^((W-1)/2)
   
    a =  a + b
    b = (b - a)*(1+i)*2^((W-1)/2)
   
    tmp = b - a
   
    b =  (tmp_re + i tmp_im)(1+i)   
      = (tmp_re + i tmp_im) + i(tmp_re + i tmp_im) 
      = tmp_re - tmp_im + i(tmp_im + tmp_re)
   
    a   = a_re + b_re + i (a_im + b_im)
    tmp = b_re - a_re + i (b_im - a_im)
    b   = tmp_re - tmp_im + i(tmp_im + tmp_re)
    */


  logic [W-1:0] a_re_q, a_im_q, tmp_re, tmp_im;



  mrsn_add_sub #(
      .N(W)
  ) adder1 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(a_re_i),
      .b_i(b_re_i),
      .sum_o(a_re_q)
  );

  mrsn_add_sub #(
      .N(W)
  ) adder2 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(a_im_i),
      .b_i(b_im_i),
      .sum_o(a_im_q)
  );



  pipe_reg #(
      .WIDTH(W),
      .DEPTH(1)
  ) pipe0 (
      .clk_i,
      .rst_ni,
      .en_i,
      .input_i (a_re_q),
      .output_o(a_re_o)
  );

  pipe_reg #(
      .WIDTH(W),
      .DEPTH(1)
  ) pipe1 (
      .clk_i,
      .rst_ni,
      .en_i,
      .input_i (a_im_q),
      .output_o(a_im_o)
  );



  mrsn_add_sub #(
      .N(W)
  ) adder3 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(b_re_i),
      .b_i(a_re_i),
      .sum_o(tmp_re)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder4 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(b_im_i),
      .b_i(a_im_i),
      .sum_o(tmp_im)
  );


  mrsn_add_sub #(
      .N(W)
  ) adder5 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i({tmp_re[P8:0], tmp_re[W-1:P8+1]}),
      .b_i({tmp_im[P8:0], tmp_im[W-1:P8+1]}),
      .sum_o(b_re_o)
  );

  mrsn_add_sub #(
      .N(W)
  ) adder6 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i({tmp_im[P8:0], tmp_im[W-1:P8+1]}),
      .b_i({tmp_re[P8:0], tmp_re[W-1:P8+1]}),
      .sum_o(b_im_o)
  );



endmodule

