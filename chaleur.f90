program chaleur

    use mod_precision
    use mod_caracteristiques
    use mod_maillage
    use mod_resolution
    use mod_sortie

    implicit none

    character(len = 50)                             :: fichier

    integer                                         :: j, i, k, n, nplot, e, nb_mailles, nb_aretes, euler
    integer, dimension(:), allocatable              :: sommets_maille, cl_arete_bord, row_CSR, col_CSR
    integer, dimension(:, :), allocatable           :: noeud_maille, ar, trig

    real(kind = pr)                                 :: t, tmax, dt, sommeDt, Fie, x1, err
    real(kind = pr), dimension(:), allocatable      :: aire_maille, l_arete, d_arete, Tn, Tnp1, val_CSR, b, T1, T2, T4
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud, milieu_arete, milieu_maille, A

! Lecture dans le fichier parameters.dat :
!       fichier : Fichier contenant le maillage
! Note : Le fichier doit etre contenu dans le dossier TYP2
    open(10, file = "parameters.dat")
        read(10, *) fichier
    close(10)

! Lecture du maillage
    call maillage("TYP2/"//fichier, nb_mailles, nb_aretes, sommets_maille, noeud_maille, coord_noeud            &
    &             , aire_maille, l_arete, d_arete, milieu_arete, milieu_maille, ar, trig, cl_arete_bord)


! Allocation des tableaux de temperature
    allocate(Tn(1:nb_mailles)) ; allocate(Tnp1(1:nb_mailles))
! Initialisation de Tnp1 et Tn
    Tn = Tinit(milieu_maille)
    Tnp1 = Tn

! Initialisation du temps
    t = 0._pr ; tmax = 15._pr

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
        n = FLOOR(tmax/dt) + 1

        do j = 1, n

            do e = 1 , nb_aretes
                i = trig(e, 1) ; k = trig(e, 2)

                if (k == 0) then
                    if ((10 <= cl_arete_bord(e)) .AND. (cl_arete_bord(e) <= 19)) then
! On est sur une arete de bord avec une condition de Dirichlet
                        Fie = -D(milieu_arete(e, :))*(Dirichlet(e, cl_arete_bord, t, milieu_arete) - Tn(i))/d_arete(e)
                        Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(i)*Fie
                    else if ((20 <= cl_arete_bord(e) .AND. cl_arete_bord(e) <= 29)) then
! On est sur une arete de bord avec une condition de Neumann
                        Fie = Neumann(e, cl_arete_bord, t, milieu_arete)
                        Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(i)*Fie
                    end if

                else
                    Fie = -D(milieu_arete(e, :))*(Tn(k) - Tn(i))/d_arete(e)
                    Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(e)*Fie
                    Tnp1(k) = Tnp1(k) + (dt/aire_maille(k))*l_arete(e)*Fie

                end if
            end do

! Ajout du terme source
            do i = 1, nb_mailles
                Tnp1(i) = Tnp1(i) + dt*Terme_source(t, milieu_maille(i, :))
            end do

            Tn = Tnp1
            t = t + dt

            nplot = FLOOR(REAL(n)/100)
            if (j == 1 .or. MODULO(j, nplot) == 0) then
                call sortie(j, Tn, sommets_maille, noeud_maille, coord_noeud)
            end if

        end do

! ----------------------------------------------------------------------------------------------
! Euler Implicite
! ----------------------------------------------------------------------------------------------
    case (2)

        dt = 0.5_pr
        n = FLOOR(tmax/dt) + 1

        ! call make_A_matrix(dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar, &
        ! &                  trig, cl_arete_bord, A)

        ! print *, "A = "
        ! do i = 1, 4
        !     print *, A(i, :)
        ! end do

        call make_A_CSR(dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar,                    &
        &               trig, cl_arete_bord, row_CSR, col_CSR, val_CSR)

        ! print *, "row_CSR", row_CSR
        ! print *, "col_CSR", col_CSR
        ! print *, "val_CSR", val_CSR

        do j = 1, n
            
            call make_b(dt, t+dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, milieu_maille,   &
            &           ar, trig, cl_arete_bord, Tn, b)

            call conjugate_gradient(row_CSR, col_CSR, val_CSR, Tn, b, Tnp1)

            Tn = Tnp1
            t = t + dt
            
            call sortie(j, Tn, sommets_maille, noeud_maille, coord_noeud)

        end do

        ! print *, "b = ", b

    case default
        print *, "Il n'y a pas de schéma temporel associé à ce nombre ! A vous d'en implémenter un :)"

    end select

    err = 0._pr
    do i = 1, nb_mailles
        err = err + aire_maille(i)*(Tnp1(i) - Sol_ex(milieu_maille(i, :)))**2
    end do

    print *, MAXVAL(l_arete), SQRT(err)

    deallocate(sommets_maille, cl_arete_bord, aire_maille, l_arete, d_arete &
    &          , noeud_maille, ar, trig, coord_noeud, milieu_arete, Tn, Tnp1)
    
end program chaleur
