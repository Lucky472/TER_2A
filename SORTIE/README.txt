N.B. : Il faut penser à changer nb_max_sommets dans mod_precision selon le maillage
que l'on utilise !

Le code marche pour les maillages :
    - test.typ2
    - mesh1_3.typ2 ; mesh2_3.typ2
    - mesh4_1.typ2 ; mesh4_1.typ2
    - mesh6.typ2 ; mesh7.typ2

    - mesh3_[i].typ2 ; i = 1, 2, 4, 5
    - mesh5.typ2 ; mesh5_reg.typ2

Le calcul diverge pour le maillage :
    - mesh3_3.typ2 : pourquoi ? :(