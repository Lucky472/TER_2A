module maillage

    use precision

    implicit none

    contains

        subroutine lecture_maillage(fichier, sommets_maille, noeud_maille, coord_noeud)

! Fichier d'entree
            character(len = *), intent(in)                              :: fichier
        
! Sorties de la subroutine
            integer, dimension(:), allocatable, intent(out)             :: sommets_maille
            integer, dimension(:, :), allocatable, intent(out)          :: noeud_maille
            real(kind = pr), dimension(: ,:), allocatable, intent(out)  :: coord_noeud

! Variables locales
            integer                                                     :: nb_noeuds, nb_mailles
            integer                                                     :: i

            open(unit = 10, file = trim(adjustl(fichier)))
! Lit la première ligne et le nombre de sommets
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
! Lit le nombre de sommets de la maille puis leur numéro
                read(10, *) sommets_maille(i), noeud_maille(i, 1:sommets_maille(i))

            end do

        end subroutine lecture_maillage
    
end module maillage