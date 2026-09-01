/*
 
 Given positive remainders r0 mod p0 and r1 mod p1 with p0 and p1
 Mersenne primes (p0 = 2^13-1 and p1 = 2^19-1), representing some
 coefficient of the polynomial product in rings Z_p0[X]/(X^n+1) and
 Z_p1[X]/(X^n+1), use the Chinese Remainder Theorem to calculate a
 signed remainder r mod p0*p1, representing the same coefficient of
 the product in Z[X].
 
 
 In order to have a unique remainder mod p0*p1 we use Newton's
 representation with signed coefficients:
 
 r = t0 + t1 * p0                                                       (1)
 
 with t0 in [-(p0-1)/2 ... (p0-1)/2], t1 in  [-(p1-1)/2 ... (p1-1)/2]
 
 so that r is automatically in [-(p0*p1-1)/2 ... (p0*p1-1)/2].
 
  */

// 7 cycles
 
`include "mrsn_ntt.svh"

module mrsn_crt #(
    parameter WIDTH = 32	       

) (
   input logic			   clk_i,
   input logic			   rst_ni,
   input logic			   en_i,
   input logic [W0-1:0]		   r0,
   input logic [W1-1:0]		   r1,
   output logic signed [WIDTH-1:0] r

);
   
   localparam			W0 = 13;
   localparam			W1 = 19;
   
   
   localparam			W10 = W1-W0;
   
   
   logic [W0-1:0]		s0, t0_n, t0_p, t0, u0;
   
   
   logic [W1-1:0]		a1, s1, t1_n, t1_p, t1, u1;
   
   logic signed [WIDTH-1:0]	v0, v1, w0, w1, rr;
   
   
   

   
   /*
    First step: - 3 cycles
    
    (r0, r1) -> Newton form (s0, s1) with r = s0 + s1 * p0 
    
          with s0 in [0 ... p0-1], s1 in [0 ... p1-1]
    
    s0 = r0
    
    s1 = mod(mod(r1 - r0, p1) * q0, p1)
    
    where q0 * p0 = 1 mod(p1)
    
    For 
    p0 = 2^13 - 1
    p1 = 2^19 - 1
    
    q0 = 8321 
    
    */

   
   //a1 = r1 - r0 mod(p1)
   
   mrsn_add_sub #(.N(W1)) adder1 (.clk_i, .rst_ni, .en_i, .mode_i(`SUB), .a_i(r1), .b_i({{W10{1'b0}},r0}), .sum_o(a1)); 


   //s1 =  mod(a1 * q0, p1) - 2 cycles
   
   mrsn_mult_8321_p19 #(.W(W1)) mult1 (.clk_i, .rst_ni, .en_i, .x_i(a1), .p_o(s1));

   
   //s0 = r0
   pipe_reg #(.WIDTH(W0), .DEPTH(3)) pipe_reg_r0 (.clk_i, .rst_ni, .en_i, .input_i(r0), .output_o(s0));


   
   /* 
    Second step - 2 cycles
    
    Signed Newton form :  (t0, t1),  r = t0 + t1 * p0,  
     
              with t0 in [-(p0-1)/2 ... (p0-1)/2], t1 in  [-(p1-1)/2 ... (p1-1)/2]
   
    
     if (s0 > (p0-1)/2) 
          t0 = s0 - p0
          t1  = mod(s1 + 1, p1); 
     else 
          t0 = s0
          t1 = s1
   
    
    if (t1 > (p1-1)/2)
        t1 = t1 - p1
    else
        t1 = t1
   
    Note: s0 has 13 bits, p0 = 13'b1. 
    
    s0 > (p0-1)/2 =  is equivalent to s0[12] = 1 
    
    If we extend with 1 bit to have signed:
     
    s0 = {0 s0}   p0 = {0 111...1}
    
    s0 - p0 as 2 complement signed is s0 + ~p0 + 1 with ~p0 = {1000...000}
    
    Therefore: s0 - p0 is:
    
    0   1 s11 s10   .... s0 + 
    1   0   0   0   ....  0 + 
                          1
    -----------------------
    
    = s0 + 1 as unsigned 13 bit. The only case with carry is s0 = 011111..111
    
    011111..11
           + 1 
    ---------- 
    1000000000
    
    but in this case s0 = p0 (=0 mod(p0)) and then s0-p0 is zero as it
    should be.
        
       
    */

   mrsn_add_sub #(.N(W1)) adder2(.clk_i, .rst_ni, .en_i, .mode_i(`ADD), .a_i(s1), .b_i({{(W1-1){1'b0}},1'b1}), .sum_o(t1_p)); 

           
   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 t0_n <= '0;
	 t1_n <= '0;
	 t0_p <= '0;
      end else begin
	 if (en_i) begin
	    t0_n <= s0;
	    t1_n <= s1;
	    t0_p <= s0 + 1;
	 end
      end
   end // always_ff @ (posedge clk_i or negedge rst_ni)



   always_comb begin
      if (t0_n[W0-1] == 1'b1) begin
         assign  t0 = t0_p;
	 assign  t1 = t1_p;
      end else begin
         assign t0 = t0_n;
	 assign t1 = t1_n;
      end
   end
   

   pipe_reg #(.WIDTH(W0), .DEPTH(1)) pipe_reg_t0 (.clk_i, .rst_ni, .en_i, .input_i(t0), .output_o(u0));
   
   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 u1 <= '0;
      end else begin
	 if (en_i) begin
	    if (t1[W1-1] == 1'b1) begin
	       u1 <= t1 + 1;
	    end else begin
	       u1 <= t1;
	    end
	 end
      end
   end // always_ff @ (posedge clk_i or negedge rst_ni)
   



   
   //sign extend u0 and v0 to WIDTH bits
 
   always_comb begin
      assign  v0[WIDTH-1:0] = {{(WIDTH-W0){u0[W0-1]}},u0[W0-1:0]};
      assign  v1[WIDTH-1:0] = {{(WIDTH-W1){u1[W1-1]}},u1[W1-1:0]};
    end


   
   /*
    
    Third step: - 2 cycles
       
   
    r = t0 + t1 * 2^13 - t1
    
    */

   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 w0 <= '0;
	 w1 <= '0;
      end else begin
	 if (en_i) begin
	    w0 <= v0 - v1;
	    w1 <= v1 <<< 13;
	 end
      end
   end // always_ff @ (posedge clk_i or negedge rst_ni)

 
   
   always_ff @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 rr   <= '0;
      end else begin
	 if (en_i) begin
	    rr   <= w0 + w1;
	 end
      end
   end // always_ff @ (posedge clk_i or negedge rst_ni)



   always_comb begin
      assign r = rr;
   end

   
endmodule
   
