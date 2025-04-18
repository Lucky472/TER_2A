program chaleur

    use mod_precision
    use mod_caracteristiques
    use mod_maillage
    use mod_resolution
    use mod_sortie

    implicit none

    character(len = 50)                             :: fichier

    integer                                         :: j, i, k, n, rep, nrep, nplot, e, nb_mailles, nb_aretes, euler, probleme
    integer, dimension(:), allocatable              :: sommets_maille, cl_arete_bord, row_CSR, col_CSR
    integer, dimension(:, :), allocatable           :: noeud_maille, ar, trig

    real(kind = pr)                                 :: t, tmax, dt, sommeDt, Fie, err, T1, T2, T4, p, tol
    real(kind = pr), dimension(:), allocatable      :: aire_maille, l_arete, d_arete, Tn, Tnp1, val_CSR, b
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud, milieu_arete, milieu_maille

! Lecture dans le fichier parameters.dat :
!       fichier : Fichier contenant le maillage
! Note : Le fichier doit etre contenu dans le dossier TYP2
    open(10, file = "parameters.dat")
        read(10, *) fichier
    close(10)

! ----------------------------------------------------------------------------------------------
! Choix du probleme
! ----------------------------------------------------------------------------------------------

    print *, "----------------------------------------------------------"
    print *, "Veuillez choisir le problème à résoudre :"
    print *, "1) Cas test 1 : Cas 1D - Neumann homogène"
    print *, "2) Cas test 2 : Flux entrant sur les bords Haut et Bas"
    print *, "3) Cas test 3 : Plaque chauffante circulaire - Température constante - CL de Dirichlet"
    print *, "4) Cas test 4 : Solution manufacturée"
    print *, "5) Cas test 5 : Conditions aux limites (Neumann et Dirichlet) périodiques"
    print *, "6) Plaque chauffante circulaire - Température altérée périodiquement - CL de Dirichlet"
    print *, "7) Pièce chauffée par un radiateur"
    print *, "8) Plaque chauffante circulaire - Température constante - Coefficient d'échange"
    print *, "9) Cas test 6 : Convergence en temps et en espace"
    print *, "----------------------------------------------------------"
    read *, probleme

! Lecture du maillage
    call maillage("TYP2/"//fichier, probleme, nb_mailles, nb_aretes, sommets_maille, noeud_maille, coord_noeud            &
    &             , aire_maille, l_arete, d_arete, milieu_arete, milieu_maille, ar, trig, cl_arete_bord)

! Allocation des tableaux de temperature
    allocate(Tn(1:nb_mailles)) ; allocate(Tnp1(1:nb_mailles))

! Initialisation du temps max
    tmax = 1._pr

! ----------------------------------------------------------------------------------------------
! Choix du schema temporel
! ----------------------------------------------------------------------------------------------

    print *, "----------------------------------------------------------"
    print *, "Veuillez choisir le schéma temporel :"
    print *, "1) Euler Explicite"
    print *, "2) Euler Implicite"
    print *, "----------------------------------------------------------"
    read *, euler

    select case (euler)

! ----------------------------------------------------------------------------------------------
! Euler explicite
! ----------------------------------------------------------------------------------------------
    case (1)

! Calcul du pas de temps
        dt = 1._pr
        do i = 1, nb_mailles
            sommeDt = 0._pr

! Calcul du denominateur de la condition CFL
            do e = 1, sommets_maille(i)
                sommeDt = sommeDt + (l_arete(ar(i, e)))*(1._pr/d_arete(ar(i, e)))
            end do

            if (1._pr/dt <= (sommeDt/aire_maille(i))) then
                dt =  1._pr/(sommeDt/aire_maille(i))
            end if

        end do

! Application du coefficient cfl sur dt
        dt = cfl*dt

! Implementation du schema
        if (probleme == 4 .OR. probleme == 6 .OR. probleme == 9) then
            nrep = 3
        else
            nrep = 1
        end if

        do rep = 1, nrep
! Initialisation de Tnp1 et Tn
            Tn = Tinit(probleme, milieu_maille)
            Tnp1 = Tn
! Remise a zero du temps
            t = 0._pr
            n = FLOOR(tmax/dt) + 1

            do j = 1, n

                do e = 1 , nb_aretes
                    i = trig(e, 1) ; k = trig(e, 2)

                    if (k == 0) then
                        if ((10 <= cl_arete_bord(e)) .AND. (cl_arete_bord(e) <= 19)) then
! On est sur une arete de bord avec une condition de Dirichlet
                            Fie = -D(probleme, milieu_arete(e, :))*                                         &
                            & (Dirichlet(probleme, e, cl_arete_bord, t, milieu_arete(e, :)) - Tn(i))/d_arete(e)
                            Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(e)*Fie
                        else if ((20 <= cl_arete_bord(e) .AND. cl_arete_bord(e) <= 29)) then
! On est sur une arete de bord avec une condition de Neumann
                            Fie = Neumann(probleme, e, cl_arete_bord, t, milieu_arete(e, :), Tn(i))
                            Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(e)*Fie
                        end if

                    else
                        Fie = -D(probleme, milieu_arete(e, :))*(Tn(k) - Tn(i))/d_arete(e)
                        Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(e)*Fie
                        Tnp1(k) = Tnp1(k) + (dt/aire_maille(k))*l_arete(e)*Fie

                    end if
                end do

! Ajout du terme source
                do i = 1, nb_mailles
                    Tnp1(i) = Tnp1(i) + dt*Terme_source(probleme, t, milieu_maille(i, :))
                end do

                Tn = Tnp1
                t = t + dt

                if (rep == 1) then
! Pour eviter d'appeler sortie trop de fois pour rien
                    nplot = FLOOR(REAL(n)/10)
                    if (MODULO(j, nplot) == 0) then
                        call sortie(j, Tn, sommets_maille, noeud_maille, coord_noeud)
                    end if
                end if

            end do

            if (probleme == 4 .OR. probleme == 6 .OR. probleme == 9) then
                err = 0._pr
                do i = 1, nb_mailles
                    err = err + aire_maille(i)*Tnp1(i)**2
                end do
                print *, dt, err

                if (rep == 1) then
                    T1 = err
                else if (rep == 2) then
                    T2 = err
                else if (rep == 3) then
                    T4 = err
                end if
                dt = dt/2
            end if
        end do

        if (probleme == 4 .OR. probleme == 6 .OR. probleme == 9) then
            p = LOG((T1 - T2)/(T2 - T4))/LOG(2._pr)
            print *, "p = ", p
        end if

! ----------------------------------------------------------------------------------------------
! Euler Implicite
! ----------------------------------------------------------------------------------------------
    case (2)

        dt = 0.05_pr
        tol = 1e-10

        call make_A_CSR(probleme, dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar,                    &
        &               trig, cl_arete_bord, row_CSR, col_CSR, val_CSR)

        if (probleme == 4 .OR. probleme == 6 .OR. probleme == 9) then
            nrep = 3
        else
            nrep = 1
        end if

        do rep = 1, nrep
! Initialisation de Tnp1 et Tn
            Tn = Tinit(probleme, milieu_maille)
            Tnp1 = Tn
! Remise a zero du temps
            t = 0._pr
! Pour obtenir le meme temps final que euler explicite dans le cas ou on a un tmax entier
! Par exemple, sans ca, lorsque tmax = 1 et dt = 0.1, on effectuait 11 iterations et on avait
! donc tfinal = 1.1 au lieu de 1
            if (ABS(tmax/dt - INT(tmax/dt)) < tol) then
                n = FLOOR(tmax/dt)
            else 
                n = FLOOR(tmax/dt) + 1
        
            end if
            do j = 1, n
                
                call make_b(probleme, dt, t+dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, milieu_maille,   &
                &           ar, trig, cl_arete_bord, Tn, b)

                call conjugate_gradient(row_CSR, col_CSR, val_CSR, Tn, b, Tnp1)

                Tn = Tnp1
                t = t + dt
                
                if (rep == 1) then
                    call sortie(j, Tn, sommets_maille, noeud_maille, coord_noeud)
                end if

            end do

            if (probleme == 4 .OR. probleme == 6 .OR. probleme == 9) then
                err = 0._pr
                do i = 1, nb_mailles
                    err = err + aire_maille(i)*Tnp1(i)**2
                end do

                if (rep == 1) then
                    T1 = err
                else if (rep == 2) then
                    T2 = err
                else if (rep == 3) then
                    T4 = err
                end if
                dt = dt/2
            end if
        end do

        if (probleme == 4 .OR. probleme == 6 .OR. probleme == 9) then
            p = LOG((T1 - T2)/(T2 - T4))/LOG(2._pr)
            print *, "p = ", p
        end if

        deallocate(row_CSR, col_CSR, val_CSR, b)

    case default
        print *, "Il n'y a pas de schéma temporel associé à ce nombre ! A vous d'en implémenter un :)"

    end select

    deallocate(sommets_maille, cl_arete_bord, aire_maille, l_arete, d_arete     &
    &          , noeud_maille, ar, trig, coord_noeud, milieu_arete, Tn, Tnp1)
    
end program chaleur
