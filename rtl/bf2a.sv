/* Modular butterfly inverse of type 2 for Mersenne prime modulus p = 2^W-1
 
 
 Input: complex 
 
    a
    b
 
 Output: complex

    a + b
   (a - b) * -i
  
 */


`define ADD 1'b0
`define SUB 1'b1


module bf2a #(
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

  /*
    
   Direct:
     a = a + i b =  a_re + i a_im + (-b_im + i b_re)
     a = a - i b =  a_re + i a_im - (-b_im + i b_re)
   
     a_re - b_im + i (a_im + b_re)
     a_re + b_im + i (a_im - b_re)

   Invers (neglecting the 1/2 factor):

     a   = a + b 
         = a_re + b_re + i ( a_im + b_im)

     i b = (a - b)  => 
    
      b  = -i (a-b) 
         = i(b-a) 
         = i b_re - b_im - (i a_re - a_im) 
         = a_im - b_im + i (b_re - a_re)

    
    */

  mrsn_add_sub #(
      .N(W)
  ) adder1 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(a_re_i),
      .b_i(b_re_i),
      .sum_o(a_re_o)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder2 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(a_im_i),
      .b_i(b_im_i),
      .sum_o(b_re_o)
  );

  mrsn_add_sub #(
      .N(W)
  ) adder3 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`ADD),
      .a_i(a_im_i),
      .b_i(b_im_i),
      .sum_o(a_im_o)
  );
  mrsn_add_sub #(
      .N(W)
  ) adder4 (
      .clk_i,
      .rst_ni,
      .en_i,
      .mode_i(`SUB),
      .a_i(b_re_i),
      .b_i(a_re_i),
      .sum_o(b_im_o)
  );


endmodule

