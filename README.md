# MersenneNTT - Work in progress!


Generic Number Theoretic Transform (NTT) using Mersenne primes for
polynomial multiplication.

It can be used for multiplication in arbitrary finite polynomial rings
Z_q[X]/(P(X)), both NTT friendly and NTT unfriendly choices, with
speeds comparable to the NTT friendly case.

Let n = deg(P(X)) and k = floor(log2(n)) + 2. 


For two polynomials a(x) and b(x) of degree n with coefficients in
Z_q, we first extend them by zeros up to degree 2^k, then do their
multiplication in Z_p[X]/(X^(2^k)-1) via MersenneNTT, then reduce
modulo both q and P(X).

The only condition is that:

((q-1)/2)^2*n < p, with p Mersenne prime.

Note that in many cases one of the polynomial has "small"
coefficients. In that case the condition is

max(abs(coef(a)) * max(abs(coef(b)) * n  < p.


Implementations:

- avx2 implementation for Mersenne prime 2^31-1, useful for typical
  polynomial multiplications as used in post-quantum cryptography (PQC).

- scalar implementation - Work in progress

- avx512 implementation - Work in progress

- Hardware implementations - Work in progress

- minimal implementation using 2^17-1 and/or 2^19-1.



