`define ROT6 6
`define ROT7 7


module mrsn_mult_8321_p19 #(
    parameter W = 19
) (
   input logic		clk_i,
   input logic		rst_ni,
   input logic		en_i,
   input logic [W-1:0]	x_i,
   output logic [W-1:0]	p_o
);

   logic [W-1:0]	x_q, w64, w65, w8320;
   

   
  always_ff @(posedge clk_i or negedge rst_ni) begin
     if (~rst_ni) begin
	x_q <= '0;
     end else begin
	if (en_i) begin
           x_q <= x_i;
	end
     end
  end

   
   assign w64   = {x_i[W-1-`ROT6 : 0], x_i[W-1 : W-`ROT6]};  // w64 = w1 << 6; 
   
   mrsn_add_sub #(
      .N(W)
  ) adder1 (
      .clk_i(clk_i),
      .rst_ni(rst_ni),
      .en_i(en_i),
      .mode_i(`ADD),
      .a_i(x_i),
      .b_i(w64),
      .sum_o(w65)
  );  // w65 = w1 + w64; 
   

  assign w8320   = {w65[W-1-`ROT7 : 0], w65[W-1 : W-`ROT7]};  // w8320 = w65 << 7; 

   mrsn_add_sub #(
      .N(W)
  ) adder2 (
      .clk_i(clk_i),
      .rst_ni(rst_ni),
      .en_i(en_i),
      .mode_i(`ADD),
      .a_i(x_q),
      .b_i(w8320),
      .sum_o(p_o)
  );  // w8321 = w1 + w8320;


   
endmodule
