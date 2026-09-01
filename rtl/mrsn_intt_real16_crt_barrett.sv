`include "mrsn_ntt.svh"


// 23 cycles	      
module mrsn_intt_real16_crt_barrett #(
    parameter  WIDTH = 32,
    parameter  LEN   = 16
) (
    input logic		     clk_i,
    input logic		     rst_ni,
    input logic		     en_i,
    input logic [WIDTH-1:0]  q,
    input logic [WIDTH-1:0]  Rq,
    input logic [WIDTH-1:0]  c_i [LEN-1 : 0],
    output logic signed [WIDTH-1:0] zq  [LEN-1 : 0]
);


   
  localparam W0 = `W0;

  localparam W1 = `W1;


   logic [WIDTH-1:0] a [LEN-1 : 0];
   

   logic signed [WIDTH-1:0] b [LEN-1 : 0];
   

   logic [WIDTH-1:0]	    q_q, Rq_q;
   
   
   pipe_reg #(.WIDTH(WIDTH), .DEPTH(16)) pipe_reg_q (.clk_i, .rst_ni, .en_i, .input_i(q), .output_o(q_q));

   pipe_reg #(.WIDTH(WIDTH), .DEPTH(16)) pipe_reg_Rq (.clk_i, .rst_ni, .en_i, .input_i(Rq), .output_o(Rq_q));

   

   mrsn_intt_real16
    #(
      .WIDTH(WIDTH),
      .LEN(LEN)
      ) 
   intt_real16
     (
      .clk_i(clk_i),
      .rst_ni(rst_ni),
      .en_i(en_i),
      .c_i(c_i),
      .a_o(a)
      );

   


   
   generate 
      for (genvar i = 0; i < LEN; i++) begin


      mrsn_crt
      #(
	.WIDTH(WIDTH)
	)
      crt (
	   .clk_i(clk_i),
	   .rst_ni(rst_ni),
	   .en_i(en_i),
	   .r0(a[i][W0-1:0]),
	   .r1(a[i][W1+W0-1:W0]),
	   .r(b[i])    
	   );

      mrsn_barrett
	#(
	  .WIDTH(WIDTH)
	  )
      barrett (
	       .clk_i(clk_i),
	       .rst_ni(rst_ni),
	       .en_i(en_i),
	       .q(q_q),
	       .Rq(Rq_q),
	       .r(b[i]),
	       .zq(zq[i])
	       );
      end      
   endgenerate

	    
   
endmodule
