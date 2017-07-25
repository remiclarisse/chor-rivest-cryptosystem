# Pollard's Rho Algorithm

This folder contains an implementation `rho-pollard` of *Pollard's rho algorithm* for finding logarithms in the finite field Z/pZ, with p prime. Each integer is encoded with 64 bits, thus the largest prime that can be used is 4294967311 (approx 2^32, because of the squaring). Above that threshold, the program will almost certainly output the wrong answer.

The file `primes` contains a list of entries for `rho-pollard`. It is formed of 2748743 lines, on which the first number is prime (`-p`) and the second number is primitive (`-g`).

> Try `rho-pollard --help` for more info.

## What can be improved

1. use integers of size 128 ou 258 bits instead of 64 bits.
2. remove the condition testing the gcd: it is useless (see Pollard paper from '78).
3. improve (?) the search for the right d-th root of unity, or have (?) a smaller d.

---
## Bibliography

John Michael Pollard, **Monte Carlo Methods for Index Computation mod p**, in *Mathematics of Computation*, 1978
