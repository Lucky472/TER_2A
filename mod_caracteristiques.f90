module mod_caracteristiques

    use mod_precision

! ----------------------------------------------------------------------------------------------
! Module definissant toutes les caracteristiques de la plaque :
!       D (Coefficient de diffusion)
!       Dirichlet (Condition aux limites de Dirichlet)
!       Neumann (Condition aux limites de Neumann)
!       Tinit (Temperature initiale dans la plaque)
!       Terme_source
!
! N.B : Le code est prevu pour pouvoir gerer 9 conditions de Dirichet et Neumann differentes
! Pour en ajouter ou changer ces dernieres, il faut aussi modifer la subroutine cl_arete
! du module mod_maille
! ----------------------------------------------------------------------------------------------

    implicit none

! Dimensions de la plaque :
!       L : Longueur
!       H : Hauteur/Largeur
    real(kind = pr), parameter                                      :: L = 1._pr
    real(kind = pr), parameter                                      :: H = 1._pr

    contains

        function D(milieu_arete) result (DXe)

! Entree de la fonction
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_arete

! Sortie de la fonction :
!       DXe : Coefficient de diffusion evalue au centre de l'arete e
            real(kind = pr)                                         :: DXe

! ----------------------------------------------------------------------------------------------
! Coefficient de diffusion uniforme
! ----------------------------------------------------------------------------------------------
            DXe = 1._pr

! ----------------------------------------------------------------------------------------------
! Coefficient de diffusion non-uniforme
! ----------------------------------------------------------------------------------------------
            ! DXe = 1._pr + 1._pr/2*SIN(2*pi*milieu_arete(1,1))*SIN(2*pi*milieu_arete(1,2))

        end function D


        function Tinit(milieu_maille) result (TinitGi)

! Entree de la fonction
            real(kind = pr), dimension(:, :), intent(in)            :: milieu_maille

! Sortie de la fonction :
!       TinitGi(1:nb_mailles) : Temperature initiale au centre de la maille i
            real(kind = pr), dimension(:), allocatable              :: TinitGi

! Variables locales
            integer                                                 :: n

            n = size(milieu_maille, 1)
            allocate(TinitGi(n))

! ----------------------------------------------------------------------------------------------
! Temperature initiale uniforme
! ----------------------------------------------------------------------------------------------
            TinitGi = 100._pr

        end function Tinit


        function Dirichlet(e, cl_arete_bord, t, milieu_arete) result (TbXe)

! Entrees de la fonction
            integer, intent(in)                                     :: e
            integer, dimension(:), intent(in)                       :: cl_arete_bord
            real(kind = pr), intent(in)                             :: t
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_arete 

! Sortie de la fonction :
!       TbXe : Temperature de bord evaluee au centre de l'arete e
            real(kind = pr)                                         :: TbXe

! ----------------------------------------------------------------------------------------------
! Cas 1D - Temperature uniforme
! ----------------------------------------------------------------------------------------------
            if (cl_arete_bord(e) == 10) then
! On est sur le bord gauche
                TbXe = 100._pr
            else if (cl_arete_bord(e) == 11) then
! On est sur le bord droit
                TbXe = 100._pr
            end if

        end function Dirichlet


        function Neumann(e, cl_arete_bord, t, milieu_arete) result (PhibXe)

! Entrees de la fonction
            integer, intent(in)                                     :: e
            integer, dimension(:), intent(in)                       :: cl_arete_bord
            real(kind = pr), intent(in)                             :: t
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_arete 

! Sortie de la fonction :
!       PhibXe : Flux de bord evaluee au centre de l'arete e
            real(kind = pr)                                         :: PhibXe

! ----------------------------------------------------------------------------------------------
! Cas 1D - Neumann homogene
! ----------------------------------------------------------------------------------------------
            if (cl_arete_bord(e) == 20) then
! On est sur le bord haut ou bas
                PhibXe = 0._pr
            end if

        end function Neumann


        function Terme_source(milieu_maille) result(SGi)

! Entree de la fonction
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_maille
! Sortie de la fonction :
!       SGi : Terme source evalue au centre de la maille i
            real(kind = pr)                                         :: SGi

! Variables locales
            real(kind = pr)                                         :: x1, x2, r
            real(kind = pr), dimension(1:2)                         :: C

! ----------------------------------------------------------------------------------------------
! Cas de la plaque chauffante circulaire
! ----------------------------------------------------------------------------------------------
! Rayon de la plaque chauffante
            r = 0.25_pr
! Coordonees du centre de la plaque chauffante
            C = 0.5_pr

            x1 = milieu_maille(1, 1) ; x2 = milieu_maille(1, 2)
            if (((x1 - C(1))**2 + (x2 - C(2))**2) <= r**2) then
                SGi = 1000._pr
            else
                SGi = 0._pr
            end if

! ----------------------------------------------------------------------------------------------
! Cas manufacture
! ----------------------------------------------------------------------------------------------
            ! SGi = 1000._pr

        end function Terme_source


        function Sol_ex(milieu_maille) result(Tex)
! Entree de la fonction
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_maille
! Sortie de la fonction :
!       SGi : Solution exacte evaluee au centre de la maille i
            real(kind = pr)                                         :: Tex

! Variables locales
            real(kind = pr)                                         :: x1, x2, S, D, Tg, Td

            S = 1000._pr ; D = 1._pr ; Tg = 100._pr ; Td = 100._pr

            x1 = milieu_maille(1, 1) ; x2 = milieu_maille(1, 2)
            Tex = -S/(2*D)*x1**2 + (Td - Tg + S/(2*D))*x1 + Tg
            
        end function Sol_ex

end module mod_caracteristiques