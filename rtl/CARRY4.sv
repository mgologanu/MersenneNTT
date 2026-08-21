/* -----\/----- EXCLUDED -----\/-----
module CARRY4
  (
   output [3:0]	CO, /-* Carry out for each stage *-/
   output [3:0]	O, /-* Carry chain data out  *-/
   input	CI, /-* Carry cascade input *-/
   input	CYINIT, /-* Carry init  *-/
   input [3:0]	DI, /-* Data input  *-/
   input [3:0]	S   /-*Select line  *-/
   );

   wire CIN = CI | CYINIT;

   muxcy muxcy0 (.O(CO[0]), .CI(CIN),   .DI(DI[0]), .S(S[0]));
   muxcy muxcy1 (.O(CO[1]), .CI(CO[0]), .DI(DI[1]), .S(S[1]));
   muxcy muxcy2 (.O(CO[2]), .CI(CO[1]), .DI(DI[2]), .S(S[2]));
   muxcy muxcy3 (.O(CO[3]), .CI(CO[2]), .DI(DI[3]), .S(S[3]));

   xorcy xorcy0 (.O(O[0]), .CI(CIN),   .LI(S[0]));
   xorcy xorcy1 (.O(O[1]), .CI(CO[0]), .LI(S[1]));
   xorcy xorcy2 (.O(O[2]), .CI(CO[1]), .LI(S[2]));
   xorcy xorcy3 (.O(O[3]), .CI(CO[2]), .LI(S[3]));

endmodule
 -----/\----- EXCLUDED -----/\----- */


module CARRY4
(
    // Carry cascade input
    input wire	      CI,
    // 
    input wire	      CYINIT,
 
 // Carry MUX data input = generate
    input wire [3:0]  DI,
   
 // Carry MUX select line = propagate
    input wire [3:0]  S,
 
 // Carry out of each stage of the chain
    output wire [3:0] CO,
    // Carry chain XOR general data out
    output wire [3:0] O
);
    wire _w_CO0 = S[0] ? CI | CYINIT : DI[0];
    wire _w_CO1 = S[1] ?      _w_CO0 : DI[1];
    wire _w_CO2 = S[2] ?      _w_CO1 : DI[2];
    wire _w_CO3 = S[3] ?      _w_CO2 : DI[3];

    assign CO   = { _w_CO3, _w_CO2, _w_CO1, _w_CO0 };

    assign O    =  S ^ { _w_CO2, _w_CO1, _w_CO0, CI | CYINIT };

endmodule
