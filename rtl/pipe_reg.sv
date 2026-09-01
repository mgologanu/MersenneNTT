module pipe_reg
  #(
    parameter WIDTH = 1,
    parameter DEPTH = 1
   )
   (
    input logic		     clk_i,
    input logic		     rst_ni,
    input logic              en_i,
    input logic [WIDTH-1:0]  input_i,
    output logic [WIDTH-1:0] output_o
   );

   localparam		    LENGTH = WIDTH*DEPTH;
   
   logic [LENGTH-1:0]	    pipe;
 

   assign output_o = pipe[WIDTH-1:0];

   always @(posedge clk_i or negedge rst_ni) begin
      if (~rst_ni) begin
	 pipe <= {LENGTH{1'b0}};
      end else begin
	  if (en_i) begin
	     if (DEPTH > 1) begin
		pipe <= pipe >> WIDTH;
	     end
	     pipe[LENGTH-1:(DEPTH-1)*WIDTH] <= input_i;
	  end
      end
   end

endmodule // pipe_reg

