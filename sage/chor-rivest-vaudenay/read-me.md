# Chor-Rivest Cryptosystem and Vaudenay's Attack

## How it works
Launch `sage` in your favorite shell and attach or load the files `chor-rivest.sage` and `attack-chor-rivest.sage`, like so:

    sage: load("chor-rivest.sage")
    sage: load("attack-chor-rivest.sage")

Define parameters for the cryptosystem: `p` a prime number and `h` a smooth integer, e.g.

    sage: p = 197
    sage: h = 24

Generate keys with `CRGenerateKeys`:

    sage: %time [PubKey, PrivKey] = CRGenerateKeys (p, h)

Pick a message with `generateRandomMessage`:

    sage: m = generateRandomMessage (p, h)

Cipher it with `CREncrypt`:

    sage: e = CREncrypt (m, PubKey)

Now `e` is the ciphertext ! To uncipher it with `CRDecrypt`:

    sage: CRDecrypt (e, PubKey, PrivKey) == m

Use Vaudenay's attack, with `VaudenayAttack`, to make an equivalent key and try to uncipher `e` with it:

    sage: %time EquivPrivKey = VaudenayAttack (PubKey)
    sage: crackedMessage = CRDecrypt (e, PubKey, EquivPrivKey)

## What can be improved
Besides everything (python meh!), many computing can be done in parallel, and the computing of
logarithms have already been parallelized !

---
### Bibliography

Benny Chor and Ronald L. Rivest, **A Knapsack-Type Public Key Cryptosystem Based on Arithmetic in Finite Fields**, in *IEEE Transactions on Information Theory*, 1988

Serge Vaudenay, **Cryptanalysis of the Chor-Rivest Cryptosystem**, in *Journal of Cryptology*, 2000
