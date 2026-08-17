//Multiplication of two complex numbers modulo Mersenne primes 13 and 17 (or 13 and 19). 

`define ADD 1'b0
`define SUB 1'b1

module mrsn_complex_multiply #(
    parameter WIDTH = 32,	       
    parameter W0 = 13,
    parameter W1 = 19
) (
    input logic		     clk_i,
    input logic		     rst_ni,
    input logic               en_i,
    input logic [WIDTH-1:0]  a_re_i,
    input logic [WIDTH-1:0]  a_im_i,
    input logic [WIDTH-1:0]  b_re_i,
    input logic [WIDTH-1:0]  b_im_i,
    output logic [WIDTH-1:0] z_re_o,
    output logic [WIDTH-1:0] z_im_o
   
);

   localparam		     W_NOTUSED = WIDTH-W0-W1;
   
   logic [W0-1:0]	     c0, s0, x0, y0;
   logic [W1-1:0]	     c1, s1, x1, y1;
   
   logic [2*W0-1:0]	     prod0_xc, prod0_xs, prod0_yc, prod0_ys;
   
   logic [W0-1:0]	     prod0_mod_xc, prod0_mod_xs, prod0_mod_ys, prod0_mod_yc, prod0_mod_x, prod0_mod_y;

   logic [2*W1-1:0]	     prod1_xc, prod1_xs, prod1_yc, prod1_ys;
   
   logic [W1-1:0]	     prod1_mod_xc, prod1_mod_xs, prod1_mod_ys, prod1_mod_yc, prod1_mod_x, prod1_mod_y;

  
   logic [WIDTH-1:0]	     z_re;
   logic [WIDTH-1:0]	     z_im;
   
   
   assign z_re_o = z_re;
   assign z_im_o = z_im;

  always_ff @(posedge clk_i or negedge rst_ni) begin
    if (~rst_ni) begin
       x0 <= '0 ;
       x1 <= '0;
       y0 <= '0;
       y1 <= '0;
       c0 <= '0;
       c1 <= '0;
       s0 <= '0;
       s1 <= '0;
    end else begin
       if(en_i) begin
	  x0 <= a_re_i[W0-1     :  0]    ;
	  x1 <= a_re_i[W0+W1-1  : W0];
	  y0 <= a_im_i[W0-1     :  0];
	  y1 <= a_im_i[W0+W1-1  : W0];
	  c0 <= b_re_i[W0-1     :  0];
	  c1 <= b_re_i[W0+W1-1  : W0];
	  s0 <= b_im_i[W0-1     :  0];
	  s1 <= b_im_i[W0+W1-1  : W0];
       end
    end
  end
   
   

   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 prod0_xc <= '0;
	 prod0_xs <= '0;
	 prod0_yc <= '0;
	 prod0_ys <= '0;
	 prod1_xc <= '0;
	 prod1_xs <= '0;
	 prod1_yc <= '0;
	 prod1_ys <= '0;
      end else begin
	 if(en_i) begin
	    prod0_xc <= x0 * c0;
	    prod0_xs <= x0 * s0;
	    prod0_yc <= y0 * c0;
	    prod0_ys <= y0 * s0;
	    prod1_xc <= x1 * c1;
	    prod1_xs <= x1 * s1;
	    prod1_yc <= y1 * c1;
	    prod1_ys <= y1 * s1;
	 end
      end
   end

   
  mrsn_add_sub #(.N(W0)) adder0_xc (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod0_xc[W0-1:0]), .b_i(prod0_xc[2*W0-1:W0]), .sum_o(prod0_mod_xc));

  mrsn_add_sub #(.N(W0)) adder0_ys (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod0_ys[W0-1:0]), .b_i(prod0_ys[2*W0-1:W0]), .sum_o(prod0_mod_ys));
    
  mrsn_add_sub #(.N(W0)) adder0_x  (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(prod0_mod_xc),     .b_i(prod0_mod_ys),        .sum_o(prod0_mod_x));

  mrsn_add_sub #(.N(W0)) adder0_xs (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod0_xs[W0-1:0]), .b_i(prod0_xs[2*W0-1:W0]), .sum_o(prod0_mod_xs));

  mrsn_add_sub #(.N(W0)) adder0_yc (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod0_yc[W0-1:0]), .b_i(prod0_yc[2*W0-1:W0]), .sum_o(prod0_mod_yc));
    
  mrsn_add_sub #(.N(W0)) adder0_y  (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod0_mod_xs),     .b_i(prod0_mod_yc),        .sum_o(prod0_mod_y));


   
  mrsn_add_sub #(.N(W1)) adder1_xc (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod1_xc[W1-1:0]), .b_i(prod1_xc[2*W1-1:W1]), .sum_o(prod1_mod_xc));

  mrsn_add_sub #(.N(W1)) adder1_ys (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod1_ys[W1-1:0]), .b_i(prod1_ys[2*W1-1:W1]), .sum_o(prod1_mod_ys));
    
  mrsn_add_sub #(.N(W1)) adder1_x  (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(prod1_mod_xc),     .b_i(prod1_mod_ys),        .sum_o(prod1_mod_x));

  mrsn_add_sub #(.N(W1)) adder1_xs (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod1_xs[W1-1:0]), .b_i(prod1_xs[2*W1-1:W1]), .sum_o(prod1_mod_xs));

  mrsn_add_sub #(.N(W1)) adder1_yc (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod1_yc[W1-1:0]), .b_i(prod1_yc[2*W1-1:W1]), .sum_o(prod1_mod_yc));
    
  mrsn_add_sub #(.N(W1)) adder1_y  (.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(prod1_mod_xs),     .b_i(prod1_mod_yc),        .sum_o(prod1_mod_y));

   
  always_ff @(posedge clk_i or negedge rst_ni) begin
    if (~rst_ni) begin
       z_re <= '0;
       z_im <= '0;
    end else begin
       if(en_i) begin
	  z_re <= {{W_NOTUSED{1'b0}}, prod1_mod_x, prod0_mod_x};
	  z_im <= {{W_NOTUSED{1'b0}}, prod1_mod_y, prod0_mod_y};
       end
    end
  end

endmodule



