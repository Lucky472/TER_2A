module maillage

    use precision
    use tri_maillage

! ----------------------------------------------------------------------------------------------
! Module contenant les subroutines pour mod_maillage :
!       maillage (prochainement)
!       lecture_maillage
!       connectivite
!       aire
! Module cree a partir du module mod_maillage.f90 cree par Luc Mieussens :
! https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK2/VF/TP/mod_maillage.f90
! ----------------------------------------------------------------------------------------------

    implicit none

    contains

        subroutine lecture_maillage(fichier, sommets_maille, noeud_maille, coord_noeud)

! Fichier d'entree
            character(len = *), intent(in)                              :: fichier
        
! Sorties de la subroutine :
!       sommets_maille : Tableau tel que sommets_maille(i) = nombre de sommets de la maille i
!       noeud_maille : Tableau tel que noeud_maille(i, 1:5) sont les au plus 5 sommets de la
! maille i
! Note : noeud_maille(i, 5) = 0 lorsque la maille i n'a que 4 sommets
!       coord_noeud : Tableau tel que coord_noeud(i, 1:2) sont les coordonees du sommet i
            integer, dimension(:), allocatable, intent(out)             :: sommets_maille
            integer, dimension(:, :), allocatable, intent(out)          :: noeud_maille
            real(kind = pr), dimension(: ,:), allocatable, intent(out)  :: coord_noeud

! Variables locales
            integer                                                     :: nb_noeuds, nb_mailles
            integer                                                     :: i

            open(unit = 10, file = trim(adjustl(fichier)))
! Lit la premiere ligne et le nombre de sommets
            read(10, *)
            read(10, *) nb_noeuds

            allocate(coord_noeud(1:nb_noeuds, 1:2))

            do i = 1, nb_noeuds
                
                read(10, *) coord_noeud(i, 1), coord_noeud(i, 2)

            end do

! Lit "cells" puis lit le nombre de mailles
            read(10, *)
            read(10, *) nb_mailles

            allocate(sommets_maille(1:nb_mailles), noeud_maille(1:nb_mailles, 5))

            do i = 1, nb_mailles
! Lit le nombre de sommets de la maille puis leur numero
                read(10, *) sommets_maille(i), noeud_maille(i, 1:sommets_maille(i))

            end do

        end subroutine lecture_maillage


        subroutine connectivite(sommets_maille, S, P, e, ar, trig)

! Tableaux d'entree
            integer, dimension(:), intent(in)                       :: sommets_maille
            integer, dimension(:, :), intent(in)                    :: S
            real(kind = pr), dimension(:, :), intent(in)            :: P

! Sorties de la subroutine :
!       e : Tableau tel que e(i, 1) et e(i, 2) soient les sommets de l'arete e de
! la maille i
!       ar : Tableau contenant le numero des aretes de chaque maille
! Note : ar(i, 5) = 0 indique que la maille a 4 aretes
!       trig : numero des mailles ayant en commun l'arete i
            integer, dimension(:, :), allocatable, intent(out)      :: e, ar, trig

! Variables locales
            integer                                                 :: nb_cotes
            integer, dimension(:, :), allocatable                   :: flag

            call fill_flag(sommets_maille, S, P, nb_cotes, flag)

            call fill_e(sommets_maille, S, nb_cotes, flag, e)

            call fill_ar_trig(sommets_maille, S, nb_cotes, flag, ar, trig)

        end subroutine connectivite


        subroutine aire(sommets_maille, S, P, aire_maille)

! Tableaux d'entree
            integer, dimension(:), intent(in)                       :: sommets_maille
            integer, dimension(:, :), intent(in)                    :: S
            real(kind = pr), dimension(:, :), intent(in)            :: P

! Sortie de la subroutine :
!       aire_maille : Tableau tel que aire_maille(i) = aire de la maille i
            real(kind = pr), dimension(:), allocatable, intent(out) :: aire_maille

! Varibles locales
            integer                                                 :: k, l, ni, nj, nb_trig
            real(kind = pr)                                         :: sommeAire
! Stocke les coordonees des sommets ni et nj            
            real(kind = pr), dimension(1:2)                         :: ai, aj

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
    
end module maillage