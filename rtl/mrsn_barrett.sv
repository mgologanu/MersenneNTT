/*
 

 The signed coefficient r obtained via CRT is further reduced to a
 signed remainder "rq" mod q, where q is an arbitrary prime in the
 range [2^9 .. 2^14-1]. The user must supply both q as a 32 bit
 integer and Rq = floor(2^31/q) again as a 32 bit integer, to be used
 in Barrett reduction.
 
 
 For the CRT, in order to have a unique remainder mod p0*p1 we use
 Newton's representation with signed coefficients:
 
 r = t0 + t1 * p0                                                       (1)
 
 with t0 in [-(p0-1)/2 ... (p0-1)/2], t1 in  [-(p1-1)/2 ... (p1-1)/2]
 
 so that r is automatically in [-(p0*p1-1)/2 ... (p0*p1-1)/2].
 
 Note that r from (1) is a 32 bit signed integer. The usual Barrett
 reduction for prime q requires -[(q-1)/2]^2 < r < [(q-1)/2]^2 and
 this is not true for q in the range [2^9 .. 2^14-1].
 
 One can use Barrett reduction with Rq = floor(2^46/q), but this
 requires multiplication of large numbers.
 
 Instead we use a two step Barrett reduction requiring 4
 multiplications with relatively small numbers (adapted to 18x25
 signed multiplier for Xilinx).
 
 Rq = round(2^33/q) - unsigned max 24 bits.
 
 1. r = t0 + t1 * 2^13 - t1 
 
 2. x = floor(r/2^14) is 18 bit signed integer (floor is equivalent to arithmetic right shift for negative)
 
 3. y = floor(x * Rq / 2^19);
 
 4. z = r - y * q;
 
 5. second barret reduction gives the final result without any check
 
 output =  z - round((z * Rq)/2^33) * q 
 
  */

// 7 cycles
 
`include "mrsn_ntt.svh"

module mrsn_barrett #(
    parameter WIDTH = 32	       

) (
   input logic			   clk_i,
   input logic			   rst_ni,
   input logic			   en_i,
 
   input logic [WIDTH-1:0]	   q,
   input logic [WIDTH-1:0]	   Rq,
   input logic signed [WIDTH-1:0]  r,
   output logic signed [WIDTH-1:0] zq

);
 
   logic signed [WIDTH-1:0]	r_q;
   
   
   logic signed [43-1:0]	y0, z0, z1, yy0, zz0, zz1, yy_halfup;
   
   logic signed [18-1:0]	x, z, z_q, q_q, yy, q_qq;

   logic signed [25-1:0]	y, Rq_q, Rq_qq;
   

   

						   
   
   pipe_reg #(.WIDTH(18), .DEPTH(1)) pipe_reg_q (.clk_i, .rst_ni, .en_i, .input_i({1'b0,q[17-1:0]}), .output_o(q_q));

   
   pipe_reg #(.WIDTH(WIDTH), .DEPTH(2)) pipe_reg_r (.clk_i, .rst_ni, .en_i, .input_i(r), .output_o(r_q));

   
   //x = floor(r/2^14);
      assign x[18-1:0] = r[18+14-1:14];
      assign Rq_q      = {1'b0,Rq[24-1:0]};
   
 
   //y = floor(x * Rq / 2^19);
   
   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 y0  <= '0;
      end else begin
	 if (en_i) begin
	    y0 <= x * Rq_q;
	 end
      end
   end

   pipe_reg #(.WIDTH(25), .DEPTH(3)) pipe_reg_Rq_q (.clk_i, .rst_ni, .en_i, .input_i(Rq_q), .output_o(Rq_qq));
   
   assign y[25-1:0] = {y0[43-1],  y0[24+19-1:19]};
   
   
   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 z0  <= '0;
	 z1  <= '0;
      end else begin
	 if (en_i) begin
	    z0 <=  y * q_q;
	    z1 <= {{(43-WIDTH){r_q[WIDTH-1]}},  r_q[WIDTH-1:0]} - z0;
	 end
      end
   end


   pipe_reg #(.WIDTH(18), .DEPTH(4)) pipe_reg_q_q (.clk_i, .rst_ni, .en_i, .input_i(q_q), .output_o(q_qq));

   
   assign z[18-1:0] = z1[18-1:0];
   

   //zz = z - round((z * Rq)/2^33) * q


   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 yy0  <= '0;
      end else begin
	 if (en_i) begin
	    yy0 <= z * Rq_qq;
	 end
      end
   end

   
   pipe_reg #(.WIDTH(18), .DEPTH(3)) pipe_reg_z (.clk_i, .rst_ni, .en_i, .input_i(z), .output_o(z_q));
   
   always_comb begin
      assign yy_halfup = yy0[(43-1):0] + { {(10){1'b0}}, 1'b1, {(43-10-1){1'b0}} };
   end

   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 yy  <= '0;
      end else begin
	 if (en_i) begin
	    yy <=  {{(18-10){yy_halfup[43-1]}}, yy_halfup[(43-1):(43-10)]};
	 end
      end
   end


   
   
   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 zz0  <= '0;
	 zz1  <= '0;
      end else begin
	 if (en_i) begin
	    zz0 <=  yy * q_qq;
	    zz1 <= {{(43-18){z_q[18-1]}},  z_q[18-1:0]} - zz0;
	 end
      end
   end

   always_comb begin
      assign zq = {{(WIDTH-18){zz1[18-1]}},zz1[18-1:0]};;
   end

   
endmodule
   
