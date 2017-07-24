# Cloner le projet dans le dossier 'home'
# Une fois Sage lancé, attacher les deux fichiers 'chor-rivest.sage' et 'attack-chor-rivest' à Sage :
sage: attach("~/chor-rivest-cryptosystem/sage/chor-rivest-vaudenay/chor-rivest.sage")
sage: attach("~/chor-rivest-cryptosystem/sage/chor-rivest-vaudenay/attack-chor-rivest.sage")
# Définir les paramètres du cryptosystème : p un nombre premier et h un entier friable, e.g.
sage: p = 197
sage: h = 24
# Générer une paire de clé grâce à 'CRGenerateKeys' :
sage: %time [PubKey, PrivKey] = CRGenerateKeys (p, h)
# Générer le message à chiffer, puis le coder, et le transformer :
sage: m = generateRandomMessage (p, h)
# Chiffrer le message :
sage: e = CREncrypt (m, PubKey)
# Maintenant e est le chiffré du message codé ! Pour le déchiffrer :
sage: CRDecrypt (e, PubKey, PrivKey) == m
# Pour l'attaque de Vaudenay :
sage: %time EquivPrivKey = VaudenayAttack (PubKey)
sage: crackedMessage = CRDecrypt (e, PubKey, EquivPrivKey)
# Le calcul des logarithmes est (grossièrement) parallelisé ...
