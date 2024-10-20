program chaleur

    use mod_precision
    use mod_maillage
    use mod_sortie

    implicit none

    character(len = 50)                             :: fichier

    integer                                         :: j, i, k, n, nplot, e, nb_mailles, nb_aretes
    integer, dimension(:), allocatable              :: sommets_maille, cl_arete_bord
    integer, dimension(:, :), allocatable           :: noeud_maille, ar, trig
    
    real(kind = pr)                                 :: Tinit, Tg, Td, Phih, Phib, D                 &
    &                                                  , t, tmax, dt, sommeDt, Fie
    real(kind = pr), dimension(:), allocatable      :: aire_maille, l_arete, d_arete, Tn, Tnp1
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud, milieu_arete

! Lecture dans le fichier parameters.dat :
!       Tinit : Condition initiale Tinit
!       Tg : Condition aux limites sur le bord gauche
!       Td : Condition aux limites sur le bord droit
!       Phih : Condition aux limites sur le bord haut
!       Phib : Condition aux limites sur le bord bas
!       D : Coefficient de diffusion
!       fichier : Fichier contenant le maillage
! Note : Le fichier doit etre contenu dans le dossier TYP2
    open(10, file = "parameters.dat")

        read(10, *) Tinit
        read(10, *) Tg
        read(10, *) Td
        read(10, *) Phih
        read(10, *) Phib
        read(10, *) D
        read(10, *) fichier

    close(10)

! Lecture du maillage
    call maillage("TYP2/"//fichier, nb_mailles, nb_aretes, sommets_maille, noeud_maille, coord_noeud &
    &             , aire_maille, l_arete, d_arete, milieu_arete, ar, trig, cl_arete_bord)


! Allocation des tableaux de temperature
    allocate(Tn(1:nb_mailles)) ; allocate(Tnp1(1:nb_mailles))
! Initialisation de Tnp1 et Tn
    Tn = Tinit
    Tnp1 = Tn  

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

! Implementation du schema
    t = 0._pr ; tmax = 1._pr
    n = FLOOR(tmax/dt) + 1

    do j = 1, n

        do e = 1 , nb_aretes
            i = trig(e, 1) ; k = trig(e, 2)

            if (k == 0) then
                if (cl_arete_bord(e)==10) then
                        Fie = -D*(Tg - Tn(i))/d_arete(e)
                        Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(i)*Fie
                else if (cl_arete_bord(e)==11) then
                        Fie = -D*(Td - Tn(i))/d_arete(e)
                        Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(i)*Fie
                else if (cl_arete_bord(e)==20) then
                        Fie = 0._pr
                end if

            else
                Fie = -D*(Tn(k) - Tn(i))/d_arete(e)
                Tnp1(i) = Tnp1(i) - (dt/aire_maille(i))*l_arete(e)*Fie
                Tnp1(k) = Tnp1(k) + (dt/aire_maille(k))*l_arete(e)*Fie

            end if
            
        end do

        Tn = Tnp1

        nplot = FLOOR(REAL(n)/10)
        if (j == 1 .or. MODULO(j, nplot) == 0) then
            call sortie(j, Tn, sommets_maille, noeud_maille, coord_noeud)
        end if

    end do

    deallocate(sommets_maille, cl_arete_bord, aire_maille, l_arete, d_arete &
    &          , noeud_maille, ar, trig, coord_noeud, milieu_arete, Tn, Tnp1)
    
end program chaleur
