/*

  NTT16 for 16 "complex" values 
 
  NTT is evaluated modulo Mersenne primes P0 = 2^13-1 and P1 = 2^19-1
  (or P1 = 2^17-1 - Work in progress)
 
  - register input    - 1 cycle
  - stage 0           - 1 cycle
  - stage 1           - 1 cycle
  - stage 2           - 5 cycles
  - stage 3           - 1 cycle
  - register output   - 1 cycle
 
  TOTAL:                10 cycles
 
  Input: 
 
  Output: 
 
 Each value packs in WIDTH bits the remainders modulo P0 and
 P1 as {zeros, rem1, rem0}
 
 */

`include "mrsn_ntt.svh"

module mrsn_ntt_complex16 #(
    parameter WIDTH = 32,
    parameter LEN   = 16
) (
    input  logic             clk_i,
    input  logic             rst_ni,
    input  logic             en_i,
    input  logic [WIDTH-1:0] a_re_i[LEN - 1 : 0],
    input  logic [WIDTH-1:0] a_im_i[LEN - 1 : 0],
    output logic [WIDTH-1:0] c_re_o[LEN - 1 : 0],
    output logic [WIDTH-1:0] c_im_o[LEN - 1 : 0]
);

  localparam W0 = `W0;

  localparam W1 = `W1;

  localparam W_NOTUSED = WIDTH - W0 - W1;

  localparam P0 = 2 ^ W0 - 1;
  localparam P1 = 2 ^ W1 - 1;

  localparam P0_om8 = (W0 - 1) / 2;
  localparam P1_om8 = (W1 - 1) / 2;

  logic [W0-1:0] c0_re_s0[LEN - 1 : 0];
  logic [W0-1:0] c0_im_s0[LEN - 1 : 0];
   
  logic [W0-1:0] c0_re_s4[LEN - 1 : 0];
  logic [W0-1:0] c0_im_s4[LEN - 1 : 0];

  logic [W1-1:0] c1_re_s0[LEN - 1 : 0];
  logic [W1-1:0] c1_im_s0[LEN - 1 : 0];

  logic [W1-1:0] c1_re_s4[LEN - 1 : 0];
  logic [W1-1:0] c1_im_s4[LEN - 1 : 0];

  genvar i;


  generate
    for (i = 0; i < LEN; i++) begin
      always_ff @(posedge clk_i or negedge rst_ni) begin
        if (~rst_ni) begin
          c0_re_s0[i] <= '0;
          c0_im_s0[i] <= '0;
          c1_re_s0[i] <= '0;
          c1_im_s0[i] <= '0;
        end else begin
          if (en_i) begin
            c0_re_s0[i] <= a_re_i[i][W0-1:0];
            c0_im_s0[i] <= a_im_i[i][W0-1:0];
            c1_re_s0[i] <= a_re_i[i][W1+W0-1:W0];
            c1_im_s0[i] <= a_im_i[i][W1+W0-1:W0];
          end
        end
      end
    end
  endgenerate

  mrsn_ntt_complex16_1 #(
      .W  (W0),
      .LEN(LEN)
  ) ntt0 (
      .clk_i,
      .rst_ni,
      .en_i,
      .a_re_i(c0_re_s0),
      .a_im_i(c0_im_s0),
      .c_re_o(c0_re_s4),
      .c_im_o(c0_im_s4)
  );

  mrsn_ntt_complex16_1 #(
      .W  (W1),
      .LEN(LEN)
  ) ntt1 (
      .clk_i,
      .rst_ni,
      .en_i,
      .a_re_i(c1_re_s0),
      .a_im_i(c1_im_s0),
      .c_re_o(c1_re_s4),
      .c_im_o(c1_im_s4)
  );



  // register output
  generate
    for (i = 0; i < LEN; i++) begin
      always_ff @(posedge clk_i or negedge rst_ni) begin
        if (~rst_ni) begin
          c_re_o[i] <= '0;
          c_im_o[i] <= '0;
        end else begin
          if (en_i) begin
            c_re_o[i] <= {{W_NOTUSED{1'b0}}, c1_re_s4[i], c0_re_s4[i]};
            c_im_o[i] <= {{W_NOTUSED{1'b0}}, c1_im_s4[i], c0_im_s4[i]};
          end
        end

      end
    end
  endgenerate

endmodule










