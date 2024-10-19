module mod_maillage

    use mod_precision
    use mod_tri_maillage

! ----------------------------------------------------------------------------------------------
! Module contenant les subroutines pour mod_maillage :
!       maillage (prochainement)
!       lecture_maillage
!       connectivite
!       aire
!       barycentre
!       cl_arete
! Module inspire du module mod_maillage.f90 cree par Luc Mieussens :
! https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK2/VF/TP/mod_maillage.f90
! ----------------------------------------------------------------------------------------------

    implicit none

    contains

        subroutine maillage(fichier, nb_mailles, nb_aretes, sommets_maille, S, P, aire_maille &
        &                   , l_arete, d_arete, milieu_arete, ar, trig, cl_arete_bord)

! Fichier d'entree
            character(len = *), intent(in)                              :: fichier

! Sorties de la subroutine :
!       nb_mailles : nombre de mailles du maillage
!       nb_aretes : nombre d'aretes du maillage
!       sommets_maille : Tableau tel que sommets_maille(i) = nombre de sommets de la maille i
!       S : Tableau tel que S(i, 1:5) soient les, au plus, 5 sommets de la maille i
! Note : S(i, 5) = 0 lorsque la maille i n'a que 4 sommets
!       P : Tableau tel que P(i, 1:2) soient les coordonnes du sommet i
!       aire_maille : Tableau tel que aire_maille(i) = aire de la maille i
!       l_arete : Tableau tel que l_arete(e) = longueur de l'arete e
!       d_arete : Tableau tel que d_arete(e) = distance d_e de Fie
!       milieu_arete : Tableau tel que milieu_arete(e, 1:2) soient les coordonnes du centre
! de l'arete e
!       ar : Tableau tel que ar(i, 1:5) soient les numeros des aretes formant la maille i
! Note : ar(i, 5) = 0 lorsque la maille i est formee par 4 aretes
!       trig : Tableau tel que trig(e, 1:2) soient le numero des mailles ayant pour arete commune
! l'arete e
! Note : trig(e, 2) = 0 lorsque que l'arete e n'a pas de maille adjacente
!       cl_arete_bord : Tableau contenant les valeurs entieres suivantes :
!           - 0 : l'arete e n'est pas une arete de bord
!           - 10 : l'arete e est soumise a la premiere condition aux limites de Dirichlet
! (ici Tg)
!           - 11 : l'arete e est soumise a la seconde condition aux limites de Dirichlet
! (ici Td)
!           - 20 : l'arete e est soumise a la condition aux limites de Neumann homogene
            integer, intent(out)                                        :: nb_mailles, nb_aretes
            integer, dimension(:), allocatable, intent(out)             :: sommets_maille, cl_arete_bord
            integer, dimension(:, :), allocatable, intent(out)          :: S, ar, trig
            real(kind = pr), dimension(:), allocatable, intent(out)     :: aire_maille, l_arete, d_arete
            real(kind = pr), dimension(:, :), allocatable, intent(out)  :: P, milieu_arete

! Variables locales
            integer                                                     :: i, nb_noeuds, tg, td
            integer, dimension(: ,:), allocatable                       :: e
            real(kind = pr), dimension(:, :), allocatable               :: P_centre
            real(kind = pr), dimension(2)                               :: A, B, Gi, Gk, Xe

! Lecture du maillage
            call lecture_maillage(fichier, sommets_maille, S, P)

            nb_mailles = size(S, 1) ; nb_noeuds = size(P, 1)

            print *, "Nombre de sommets :", nb_noeuds
            print *, "Nombre de mailles :", nb_mailles

! Connectivite
            call connectivite(sommets_maille, S, P, e, ar, trig)

            nb_aretes = size(e, 1)

! Aire des mailles
            call aire(sommets_maille, S, P, aire_maille)

! Barycentre des mailles
            call barycentre(sommets_maille, S, P, P_centre)

! Longueur d_e
            allocate(d_arete(1:nb_aretes))
            d_arete = 0._pr

            do i = 1, nb_aretes

                tg = trig(i, 1) ; td = trig(i, 2)
                Gi = P_centre(tg, :)

                if (td /= 0) then
                    Gk = P_centre(td, :)
                    d_arete(i) = SQRT(SUM((Gk - Gi)**2))
                else
                    A = P(e(i, 1), :)
                    B = P(e(i, 2), :)
                    Xe = (A + B)/2._pr
                    d_arete(i) = SQRT(SUM((Xe - Gi)**2))
                end if

            end do

! Longueur des aretes et coordonnes du centre de l'arete Xe
            allocate(l_arete(1:nb_aretes)) ; allocate(milieu_arete(1:nb_aretes, 2))

            do i = 1, nb_aretes

                A = P(e(i, 1), :)
                B = P(e(i, 2), :)

                l_arete(i) = SQRT(SUM((B - A)**2))
                milieu_arete(i, :) = (A + B)/2._pr

            end do

! Condtions aux limites sur les aretes de bord
            call cl_arete(P, e, trig, cl_arete_bord)

        end subroutine maillage


        subroutine lecture_maillage(fichier, sommets_maille, S, P)

! Fichier d'entree
            character(len = *), intent(in)                              :: fichier
        
! Sorties de la subroutine :
!       sommets_maille : Tableau tel que sommets_maille(i) = nombre de sommets de la maille i
!       S : Tableau tel que S(i, 1:5) sont les au plus 5 sommets de la maille i
! Note : S(i, 5) = 0 lorsque la maille i n'a que 4 sommets
!       P : Tableau tel que P(i, 1:2) sont les coordonees du sommet i
            integer, dimension(:), allocatable, intent(out)             :: sommets_maille
            integer, dimension(:, :), allocatable, intent(out)          :: S
            real(kind = pr), dimension(: ,:), allocatable, intent(out)  :: P

! Variables locales
            integer                                                     :: nb_noeuds, nb_mailles
            integer                                                     :: i

            open(unit = 10, file = trim(adjustl(fichier)))
! Lit la premiere ligne et le nombre de sommets
            read(10, *)
            read(10, *) nb_noeuds

            allocate(P(1:nb_noeuds, 1:2))

            do i = 1, nb_noeuds
                
                read(10, *) P(i, 1), P(i, 2)

            end do

! Lit "cells" puis lit le nombre de mailles
            read(10, *)
            read(10, *) nb_mailles

            allocate(sommets_maille(1:nb_mailles), S(1:nb_mailles, 5))

            do i = 1, nb_mailles
! Lit le nombre de sommets de la maille puis leur numero
                read(10, *) sommets_maille(i), S(i, 1:sommets_maille(i))

            end do

        end subroutine lecture_maillage


        subroutine connectivite(sommets_maille, S, P, e, ar, trig)

! Tableaux d'entree
            integer, dimension(:), intent(in)                           :: sommets_maille
            integer, dimension(:, :), intent(in)                        :: S
            real(kind = pr), dimension(:, :), intent(in)                :: P

! Sorties de la subroutine :
!       e : Tableau tel que e(i, 1) et e(i, 2) soient les sommets de l'arete e de
! la maille i
!       ar : Tableau contenant le numero des aretes de chaque maille
! Note : ar(i, 5) = 0 indique que la maille a 4 aretes
!       trig : numero des mailles ayant en commun l'arete i
            integer, dimension(:, :), allocatable, intent(out)          :: e, ar, trig

! Variables locales
            integer                                                     :: nb_cotes
            integer, dimension(:, :), allocatable                       :: flag

            call fill_flag(sommets_maille, S, P, nb_cotes, flag)

            call fill_e(sommets_maille, S, nb_cotes, flag, e)

            call fill_ar_trig(sommets_maille, S, nb_cotes, flag, ar, trig)

        end subroutine connectivite


        subroutine aire(sommets_maille, S, P, aire_maille)

! Tableaux d'entree
            integer, dimension(:), intent(in)                           :: sommets_maille
            integer, dimension(:, :), intent(in)                        :: S
            real(kind = pr), dimension(:, :), intent(in)                :: P

! Sortie de la subroutine :
!       aire_maille : Tableau tel que aire_maille(i) = aire de la maille i
            real(kind = pr), dimension(:), allocatable, intent(out)     :: aire_maille

! Varibles locales
            integer                                                     :: k, l, ni, nj, nb_trig
            real(kind = pr)                                             :: sommeAire
! Stocke les coordonees des sommets ni et nj            
            real(kind = pr), dimension(1:2)                             :: ai, aj

            nb_trig = size(S, 1)

            allocate(aire_maille(1:nb_trig))

            do k = 1, nb_trig
                sommeAire = 0._pr
                do l = 1, sommets_maille(k)

                    call get_Si1_back(k, l, sommets_maille, S, ni, nj)

                    ai = P(ni, :) ; aj = P(nj, :)

                    sommeAire = sommeAire + (1._pr/2)*((ai(1)*aj(2) - aj(1)*ai(2)))

                end do
                aire_maille(k) = sommeAire
            end do

        end subroutine aire


        subroutine barycentre(sommets_maille, S, P, barycentre_maille)

! Tableaux d'entree
            integer, dimension(:), intent(in)                           :: sommets_maille
            integer, dimension(:, :), intent(in)                        :: S
            real(kind = pr), dimension(:, :), intent(in)                :: P

! Sortie de la subroutine :
!       barycentre_maille : Tableau tel que barycentre_maille(i) = barycentre de
! la maille i
            real(kind = pr), dimension(:, :), allocatable, intent(out)  :: barycentre_maille

! Variables locales
            integer                                                     :: k, l, ni, nj, nb_trig
            real(kind = pr)                                             :: sommeXnum, sommeYnum, sommeDeno
            real(kind = pr), dimension(2)                               :: ai, aj

            nb_trig = size(S, 1)

            allocate(barycentre_maille(1:nb_trig, 2))

            do k = 1, nb_trig
                sommeXnum = 0._pr ; sommeYnum = 0._pr ; sommeDeno = 0._pr
                do l = 1, sommets_maille(k)

                    call get_Si1_back(k, l, sommets_maille, S, ni, nj)

                    ai = P(ni, :) ; aj = P(nj, :)

                    sommeXnum = sommeXnum + (ai(1) + aj(1))*((ai(1)*aj(2)) - (ai(2)*aj(1)))
                    sommeYnum = sommeYnum + (ai(2) + aj(2))*((ai(1)*aj(2)) - (ai(2)*aj(1)))
                    sommeDeno = sommeDeno + 3._pr*((ai(1)*aj(2)) - (ai(2)*aj(1)))

                end do
                barycentre_maille(k, 1) = sommeXnum/sommeDeno
                barycentre_maille(k, 2) = sommeYnum/sommeDeno
            end do

        end subroutine barycentre


        subroutine cl_arete(P, e, trig, cl_arete_bord)

! Entrees de la subroutine
            integer, dimension(:, :), intent(in)                        :: e, trig
            real(kind = pr), dimension(:, :), intent(in)                :: P

! Sortie de la subroutine :
!       cl_arete_bord : Tableau contenant les valeurs entieres suivantes :
!           - 0 : l'arete e n'est pas une arete de bord
!           - 10 : l'arete e est soumise a la premiere condition aux limites de Dirichlet
! (ici Tg)
!           - 11 : l'arete e est soumise a la seconde condition aux limites de Dirichlet
! (ici Td)
!           - 20 : l'arete e est soumise a la condition aux limites de Neumann homogene
            integer, dimension(:), allocatable, intent(out)             :: cl_arete_bord

! Variables locales
            integer                                                     :: i, nb_aretes, ni, nj
            real(kind = pr), dimension(2)                               :: ai, aj

            nb_aretes = size(e, 1)

            allocate(cl_arete_bord(1:nb_aretes))

            cl_arete_bord = 0

            do i = 1, nb_aretes

                if (trig(i, 2) == 0) then
                    ni = e(i, 1) ; nj = e(i, 2)
                    ai = P(ni, :) ; aj = P(nj, :)

! Teste les differents cotes de maille i
                    if (ai(1) == 0._pr .and. aj(1) == 0._pr) then           ! Bord gauche
                        cl_arete_bord(i) = 10
                    else if (ai(1) == 1._pr .and. aj(1) == 1._pr) then      ! Bord droit
                        cl_arete_bord(i) = 11
                    else if (ai(2) == 0._pr .and. aj(2) == 0._pr          & ! Bord bas
                    &       .or. ai(2) == 1._pr .and. aj(2) == 1._pr) then  ! Bord haut
                        cl_arete_bord(i) = 20
                    end if

                end if

            end do

        end subroutine cl_arete
    
end module mod_maillage