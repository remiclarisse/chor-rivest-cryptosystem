# Chor-Rivest Cryptosystem and Vaudenay's Attack

## How it works
Launch `sage` in your favorite shell and attach or load the files `chor-rivest.sage` and `attack-chor-rivest.sage`, like so :
    sage: load("chor-rivest.sage")
    sage: load("attack-chor-rivest.sage")
 Définir les paramètres du cryptosystème : p un nombre premier et h un entier friable, e.g.
    sage: p = 197
    sage: h = 24
Générer une paire de clé grâce à 'CRGenerateKeys' :
    sage: %time [PubKey, PrivKey] = CRGenerateKeys (p, h)
Générer le message à chiffer, puis le coder, et le transformer :
    sage: m = generateRandomMessage (p, h)
Chiffrer le message :
    sage: e = CREncrypt (m, PubKey)
Maintenant e est le chiffré du message codé ! Pour le déchiffrer :
    sage: CRDecrypt (e, PubKey, PrivKey) == m
Pour l'attaque de Vaudenay :
    sage: %time EquivPrivKey = VaudenayAttack (PubKey)
    sage: crackedMessage = CRDecrypt (e, PubKey, EquivPrivKey)
Le calcul des logarithmes est (grossièrement) parallelisé ...
