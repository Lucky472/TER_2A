N.B. : Il faut penser à changer nb_max_sommets et cfl dans mod_precision selon le maillage
que l'on utilise (si l'on veut essayer de pas avoir trop de 0 dans la matrice de sommets S,
sinon on peut fixer nb_max_sommets = 6 et cela marche pour tous les maillages de la liste
ci-dessous) !

La condition cfl = 1 marche pour les maillages :
    - test.typ2
    - mesh1_3.typ2 ; mesh2_3.typ2
    - mesh4_1.typ2 ; mesh4_1.typ2
    - mesh6.typ2 ; mesh7.typ2

    - mesh3_[i].typ2 ; i = 1, 2, 4, 5
    - mesh5.typ2 ; mesh5_reg.typ2

La condition cfl = 0.95 marche pour le maillage :
    - mesh3_3.typ2
Note : Le schema V.F diverge losque cfl = 1 pour ce maillage