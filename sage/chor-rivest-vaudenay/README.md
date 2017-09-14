# Chor-Rivest Cryptosystem and Vaudenay's Attack

## How does it work?

Launch SageMath in your favorite shell and attach or load the files `chor-rivest-prime-number.sage` and `vaudenay-attack.sage`, like so:

    sage: load("chor-rivest-prime-number.sage")
    sage: load("vaudenay-attack.sage")

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


### Vaudenay's attack

Use Vaudenay's attack, with `VaudenayAttack`, to make an equivalent key then try to uncipher `e` with it:

    sage: %time EquivPrivKey = VaudenayAttack (PubKey)
    sage: crackedMessage = CRDecrypt (e, PubKey, EquivPrivKey)
    sage: crackedMessage == m

---
## Bibliography

Benny Chor and Ronald L. Rivest, **A Knapsack-Type Public Key Cryptosystem Based on Arithmetic in Finite Fields**, in *IEEE Transactions on Information Theory*, 1988

Serge Vaudenay, **Cryptanalysis of the Chor-Rivest Cryptosystem**, in *Journal of Cryptology*, 2000
