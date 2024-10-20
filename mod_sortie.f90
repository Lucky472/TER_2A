module mod_sortie

    use mod_precision

! ----------------------------------------------------------------------------------------------
! Module contenant les subroutines de sortie :
!       sortie
!       compte_cellule
! Module inspire du module mod_sortie.f90 cree par Luc Mieussens :
! https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK2/VF/TP/mod_sortie.f90
! et de la documentation du format VTK :
! https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
! ----------------------------------------------------------------------------------------------

    implicit none

    contains

        subroutine sortie(iter, T, sommets_maille, S, P)

! Entrees de la subroutine
            integer, intent(in)                                 :: iter
            integer, dimension(:), intent(in)                   :: sommets_maille
            integer, dimension(:, :), intent(in)                :: S
            real(kind = pr), dimension(:), intent(in)           :: T
            real(kind = pr), dimension(:, :), intent(in)        :: P

! Variables locales
            character(len = 30)                                 :: it_ch
            integer                                             :: i, nb_noeuds, nb_mailles, nb_max_cellules

            nb_noeuds = size(P, 1) ; nb_mailles = size(S, 1)
            write(it_ch, *) iter

            open(unit = 1, file = 'SORTIE/sortie_'//trim(adjustl(it_ch))//'.vtk')
  
            write(1, '(1A26)') '# vtk DataFile Version 2.0'
            write(1, *) 'test maillage'
            write(1, *) 'ASCII'
            write(1, *) 'DATASET UNSTRUCTURED_GRID'

            write(1, *) 'POINTS', nb_noeuds, ' double'
            do i= 1, nb_noeuds
                write(1, *) P(i, 1), P(i, 2), 0.0_pr
            end do

            nb_max_cellules = compte_cellule(sommets_maille)
            write(1, *) 'CELLS ', nb_mailles, nb_max_cellules
            do i = 1, nb_mailles
                write(1, *) sommets_maille(i), S(i, 1:sommets_maille(i)) - 1
            end do
  
            write(1,*) 'CELL_TYPES ', nb_mailles
            do i = 1, nb_mailles

! Evalue la forme de la cellule
                if (sommets_maille(i) >= 5) then
! Dans ce cas, on a un polygone :
!       VTK_POLYGON = 7
                    write(1, *) 7
                else if (sommets_maille(i) == 4) then
! Dans ce cas, on a un quadrilataire :
!       VTK_QUAD = 9
                    write(1, *) 9
                else if (sommets_maille(i) == 3) then
! Dans ce cas, on a un triangle :
!       VTK_TRIANGLE = 5
                    write(1, *) 5
                end if

            end do

! Transmet les informations sur les donnees connues au logiciel de visualisation
            write(1, *) 'CELL_DATA ', nb_mailles
            write(1, *) "SCALARS T double"
            write(1, *) "LOOKUP_TABLE default"

            do i = 1, nb_mailles
                write(1, *) T(i)
            end do

            close(1)

        end subroutine sortie

        function compte_cellule(sommets_maille) result(nb_max_cellules)

! Tableau d'entree de la fonction
            integer, dimension(:), intent(in)                   :: sommets_maille

! Sortie de la fonction :
!       nb_max_cellules : Taille de la liste de cellule i.e. la somme du nombre de sommets
! de toutes les mailles auquel on ajoute la taille du tableau sommets_maille
            integer                                             :: nb_max_cellules

! Variables locales
            integer                                             :: i, nb_mailles

            nb_max_cellules = 0
            nb_mailles = size(sommets_maille, 1)

! Somme de tous les sommets sur toutes les mailles
            do i = 1, nb_mailles

                nb_max_cellules = nb_max_cellules + sommets_maille(i)

            end do

! Ajout de la taille du tableau sommets_maille
            nb_max_cellules = nb_max_cellules + nb_mailles

        end function compte_cellule

end module mod_sortie