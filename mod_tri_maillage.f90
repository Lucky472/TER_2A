module tri_maillage

    use precision

    implicit none

! ----------------------------------------------------------------------------------------------
! Module contenant les subroutines de tri pour mod_maillage :
!       fill_flag
!       fill_e
!       fill_ar_trig
! Module cree a partir du module mod_maillage.f90 cree par Luc Mieussens
! ----------------------------------------------------------------------------------------------
    
    contains

        subroutine get_Si1_back(k, l, sommets_maille, S, ni, nj)

! Entiers d'entree pour la subroutine
            integer, intent(in)                                 :: k, l

! Tableaux d'entree pour la subroutine
            integer, dimension(:), intent(in)                   :: sommets_maille
            integer, dimension(:, :), intent(in)                :: S

! Sorties de la subroutine : les sommets ni = S(k, l) et nj = S(k, l+1)
! Retourne nj = S(k, 1) si l = sommets_maille(k)
            integer, intent(out)                                :: ni, nj

            if (l /= sommets_maille(k)) then
                ni = S(k, l); nj = S(k, l+1)
            else
                ni = S(k, l); nj = S(k, 1)
            end if

        end subroutine get_Si1_back


        subroutine fill_flag(sommets_maille, S, P, nb_cotes, flag)

! Tableaux d'entree pour la subroutine (tableau des sommets, des coordonees et du nombre de sommets
! sur la maille)
            integer, dimension(:), intent(in)                   :: sommets_maille
            integer, dimension(:, :), intent(in)                :: S
            real(kind = pr), dimension(:, :), intent(in)        :: P

! Tableau de sortie (contient le numero de l'arete liant les deux sommets)
            integer, dimension(:, :), allocatable, intent(out)  :: flag
! Nombre d'aretes (necessaire pour initialiser e, ar et trig)
            integer, intent(out)                                :: nb_cotes

! Variables locales (ni, nj correspondent aux n1, n2, n3 dans le code de Luc)
            integer                                             :: k, l, ni, nj, nb_trig, nb_noeuds

! Nombre de mailles
            nb_trig = size(S, 1)
! Nombre de noeuds
            nb_noeuds = size(P, 1)

! Le coefficient flag(ni, nj) contient l'arete liant les sommets ni et nj en fin de boucle
            allocate(flag(1:nb_noeuds, 1:nb_noeuds))

! Calcul du nombre d'aretes nb_cotes
            nb_cotes = 0
            flag = 0

            do k = 1, nb_trig
                do l = 1, sommets_maille(k)

                    call get_Si1_back(k, l, sommets_maille, S, ni, nj)

! Marque l'arete entre les sommets ni et nj si elle n'est pas marquee
                    if(flag(ni, nj)==0 .or. flag(nj, ni)==0) then
                        nb_cotes = nb_cotes + 1
                        flag(ni, nj) = nb_cotes; flag(nj, ni) = nb_cotes
                    end if

                end do
            end do
        
        end subroutine fill_flag


        subroutine fill_e(sommets_maille, S, nb_cotes, flag, e)

! Entrees de la subroutine
            integer, dimension(:), intent(in)                   :: sommets_maille
            integer, dimension(:,:), intent(in)                 :: S, flag
            integer, intent(in)                                 :: nb_cotes

! Sortie de la subroutine : tableau e tel que e(i, 1) et e(i, 2) soient les
! noeuds a chaque extremite de l'arete
            integer, dimension(:, :), allocatable, intent(out)  :: e

! Varibles locales
            integer                                             :: k, l, ni, nj, nb_trig

            nb_trig = size(S, 1)

            allocate(e(1:nb_cotes, 2))

            do k = 1, nb_trig
                do l = 1, sommets_maille(k)

                    call get_Si1_back(k, l, sommets_maille, S, ni, nj)

                    e(flag(ni, nj), 1) = ni ; e(flag(ni, nj), 2) = nj 

                end do
            end do

        end subroutine fill_e


        subroutine fill_ar_trig(sommets_maille, S, nb_cotes, flag, ar, trig)

! Entrees de la subroutine
            integer, dimension(:), intent(in)                   :: sommets_maille
            integer, dimension(:, :), intent(in)                :: S, flag
            integer, intent(in)                                 :: nb_cotes

! Sorties de la subroutine :
!       ar : numero des aretes de chaque maille
!       trig : numero des mailles ayant en commun l'arete i
            integer, dimension(:, :), allocatable, intent(out)  :: ar, trig

! Variables locales
            integer                                             :: k, l, ni, nj, nb_trig

            nb_trig = size(S, 1)

            allocate(ar(1:nb_trig, 5), trig(1:nb_cotes, 2))

! Remplissage de ar et trig
            ar = 0
            trig = 0

            do k = 1, nb_trig
                do l = 1, sommets_maille(k)

                    call get_Si1_back(k, l, sommets_maille, S, ni, nj)

! Recupere le numero associe a l'arete liant les sommets ni et nj
                    ar(k, l) = flag(ni, nj)

                    if (trig(ar(k, l), 1) /= 0) then
                        trig(ar(k, l), 2) = k
                    else
                        trig(ar(k, l), 1) = k
                    end if

                end do
            end do

        end subroutine fill_ar_trig

end module tri_maillage