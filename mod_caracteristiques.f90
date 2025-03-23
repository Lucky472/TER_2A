module mod_caracteristiques

    use mod_precision

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
                TbXe = 300._pr
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

end module mod_caracteristiques