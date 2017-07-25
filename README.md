# A study of the Chor-Rivest Cryptosystem

This work is/was done at the INRIA Saclay-Île de France by Rémi Clarisse (me writting) during an intership and supervised by Daniel Augot and Luca De Feo.

## The Discrete Logarithm Problem

Because Chor-Rivest cryptosystem uses discreet logarithms in finite fields, we have encountered severeal algorithms solving the discreet logarithm problem.

### Polhig-Hellman Method and Baby-Step-Giant-Step Algorithm

The `./sage/DLP/` folder contains an implementation of *Pohlig-Hellman* : Given a field Z/pZ and **knowing the factorisation of p-1** (strong assumption), solve the DLP for each maximal sub-r-group of (Z/pZ)\*, where r is a divisor of p-1, and reconstruct the original log. `pohlig-hellman.sage` uses an exhaustive seach in the cyclic sub-r-groups. There is also an implementation of the *baby-step giant-step algorithm* for solving the DLP.

#### What can be improved:

1. use Pohlig-Hellman with field extension of Z/pZ.
2. use BSGS in Pohlig-Hellman.

### Pollard's Rho Algorithm

The `./c/` folder contains an implementation `rho-pollard` of *Pollard's rho algorithm* for finding logarithms in the field Z/pZ, with p prime.

## The Chor-Rivest Cryptosystem

### Cryptosystem

The `./sage/chor-rivest-vaudenay/` folder contains an implementation `chor-rivest.sage` of the *Chor-Rivest cryptosystem*. Be thorough when choosing the parameters : p is a prime number and h is a smooth number. The **DLP must be reasonanbly solvable in GF(p^h)** : the generation of the public key depends on it. Regarding the usage, the names are eloquent enough : `CRGenerateKeys(p, h)` returns `[PubKey, PrivKey]`, `CREncrypt(m, PubKey)` encrypts the message m (of length p and Hamming's weight h) and `CRDecrypt(e, Keys)` gives back m, provided that e is the ciphertext of m.

### Attacks

The `./sage/chor-rivest-vaudenay/attacks-chor-rivest.sage` contains sevral implementations of proposed attacks against the Chor-Rivest cryptosystem.

## The logbook, a.k.a the tex file

The tex file `./tex/chor-rivest.tex` is kind of a logbook of the work and research done, and 'stuff' learned along the way.

> It is in french (:
