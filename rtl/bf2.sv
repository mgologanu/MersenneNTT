`define ADD  1'b0
`define SUB  1'b1


module bf2 #(
  parameter  W = 13
) (
   input logic		clk_i,
   input logic		rst_ni,
   input logic		en_i,

   input logic [W-1:0]	a_re_i,
   input logic [W-1:0]	a_im_i,

   input logic [W-1:0]	b_re_i,
   input logic [W-1:0]	b_im_i,

   output logic [W-1:0]	a_re_o,
   output logic [W-1:0]	a_im_o,

   output logic [W-1:0]	b_re_o,
   output logic [W-1:0]	b_im_o
   
   );

   //  a = a + i b =  a_re + i a_im + (-b_im + i b_re)
   //  a = a - i b =  a_re + i a_im - (-b_im + i b_re)
   
   //  a_re - b_im + i (a_im + b_re)
   //  a_re + b_im + i (a_im - b_re)
   

   mrsn_add_sub #(.N(W)) adder1  (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(a_re_i), .b_i(b_im_i), .sum_o(a_re_o));
   mrsn_add_sub #(.N(W)) adder2  (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(a_re_i), .b_i(b_im_i), .sum_o(b_re_o));

   mrsn_add_sub #(.N(W)) adder3  (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(a_im_i), .b_i(b_re_i), .sum_o(a_im_o));
   mrsn_add_sub #(.N(W)) adder4  (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(a_im_i), .b_i(b_re_i), .sum_o(b_im_o));
   

endmodule
	     
