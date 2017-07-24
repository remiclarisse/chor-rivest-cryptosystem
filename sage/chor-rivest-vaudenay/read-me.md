# Chor-Rivest Cryptosystem and Vaudenay's Attack

## How it works
Launch `sage` in your favorite shell and attach or load the files `chor-rivest.sage` and `attack-chor-rivest.sage`, like so:

    sage: load("chor-rivest.sage")
    sage: load("attack-chor-rivest.sage")

Define parameters for the cryptosystem: p a prime number and h a smooth integer, e.g.

    sage: p = 197
    sage: h = 24

Generate keys with `CRGenerateKeys`:

    sage: %time [PubKey, PrivKey] = CRGenerateKeys (p, h)

Pick a message:

    sage: m = generateRandomMessage (p, h)

Cipher it:

    sage: e = CREncrypt (m, PubKey)

Now `e` is the ciphertext ! To uncipher it :

    sage: CRDecrypt (e, PubKey, PrivKey) == m

Use Vaudenay's attack to make an equivalent key and try to uncipher `e` with it:

    sage: %time EquivPrivKey = VaudenayAttack (PubKey)
    sage: crackedMessage = CRDecrypt (e, PubKey, EquivPrivKey)

## What can be improved
Besides everything (python meh!), many computing can be done in parallel, and the computing of
logarithms have already been parallelized !
