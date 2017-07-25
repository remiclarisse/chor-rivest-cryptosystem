# Chor-Rivest Cryptosystem with Joux Enhancement

## How it works

Launch `sage` in your favorite shell and attach or load the file `chor-rivest.sage`, like so:

    sage: load("chor-rivest.sage")

Define parameters for the cryptosystem: `q` a prime-power of 2 (3 may also work) and `k` a prime number, e.g.

    sage: q = 2 ** 7
    sage: k = 29

Generate keys with `CRGenerateKeys`:

    sage: %time [PubKey, PrivKey] = CRGenerateKeys (q, k)

Pick a message with `generateRandomMessage`:

    sage: m = generateRandomMessage (q, k)

Cipher it with `CREncrypt`:

    sage: e = CREncrypt (m, PubKey)

Now `e` is the ciphertext ! To uncipher it with `CRDecrypt`:

    sage: CRDecrypt (e, PubKey, PrivKey) == m

### Remark

On the construction of the finite field K.

---
## Bibliography

Benny Chor and Ronald L. Rivest, **A Knapsack-Type Public Key Cryptosystem Based on Arithmetic in Finite Fields**, in *IEEE Transactions on Information Theory*, 1988

Antoine Joux, **A new index calculus algorithm with complexity L(1/4+o(1)) in small characteristic**, 2013
