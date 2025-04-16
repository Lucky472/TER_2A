module mod_caracteristiques

    use mod_precision

! ----------------------------------------------------------------------------------------------
! Module definissant toutes les caracteristiques de la plaque :
!       D (Coefficient de diffusion)
!       Dirichlet (Condition aux limites de Dirichlet)
!       Neumann (Condition aux limites de Neumann)
!       Heat_transfer_neumann (Condition au limite implicite de Neumann pour le coefficient
! de transfert thermique h)
!       Tinit (Temperature initiale dans la plaque)
!       Terme_source
!
! N.B : Le code est prevu pour pouvoir gerer 9 conditions de Dirichet et Neumann differentes
! Pour en ajouter ou changer ces dernieres, il faut aussi modifer la subroutine cl_arete
! du module mod_maillage
! ----------------------------------------------------------------------------------------------

    implicit none

! Dimensions de la plaque :
!       L : Longueur
!       H : Hauteur/Largeur
    real(kind = pr), parameter                                      :: L = 1._pr
    real(kind = pr), parameter                                      :: H = 1._pr

! ----------------------------------------------------------------------------------------------
! Piece chauffee par un radiateur
! ----------------------------------------------------------------------------------------------
!       hb : Hauteur du bloc de beton
!       hv : Hauteur de la vitre en verre (1 - 2*hb)
!       eb : Epaisseur du bloc de beton
!       ev : Epaisseur d'un vitrage
!       ea : Epaisseur d'air entre chaque vitrage
!       dm : Distance entre le bord gauche et le bloc de beton
!       Tint : Temperature a l'interieur de la piece
!       Text : Temperature a l'exterieur de la piece
!       Lrad : Longueur du radiateur
!       Hrad : Hauteur du radiateur
!       Hsol : Distance entre le sol et le bas du radiateur
!       drad : Distance entre le mur et le radiateur
!       trad : Temps pendant lequel le radiateur emet une source de chaleur
    real(kind = pr), parameter                                      :: hb = 0.2_pr
    real(kind = pr), parameter                                      :: hv = 0.6_pr
    real(kind = pr), parameter                                      :: eb = 0.3_pr
    real(kind = pr), parameter                                      :: ev = 0.003_pr
    real(kind = pr), parameter                                      :: ea = 0.009_pr
    real(kind = pr), parameter                                      :: dm = 0.05_pr
    real(kind = pr), parameter                                      :: Tint = 290._pr
    real(kind = pr), parameter                                      :: Text = 278._pr
    real(kind = pr), parameter                                      :: Lrad = 0.3_pr
    real(kind = pr), parameter                                      :: Hrad = 0.2_pr
    real(kind = pr), parameter                                      :: Hsol = 0.05_pr
    real(kind = pr), parameter                                      :: drad = 0.05_pr
    real(kind = pr), parameter                                      :: trad = 5_pr

! ----------------------------------------------------------------------------------------------
! Neumann - Coefficient de transfert thermique
! ----------------------------------------------------------------------------------------------
!       ht : Coefficient de transfert thermique (suppose constant ici)
!       Tinf : Temperature loin de la plaque
    real(kind = pr), parameter                                      :: ht = 10._pr
    real(kind = pr), parameter                                      :: Tinf = 293._pr

    contains

        function D(probleme, milieu_arete) result (DXe)

! Entrees de la fonction
            integer, intent(in)                                     :: probleme
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_arete

! Sortie de la fonction :
!       DXe : Coefficient de diffusion evalue au centre de l'arete e
            real(kind = pr)                                         :: DXe

! Variables locales
            real(kind = pr)                                         :: e, x1, x2

            x1 = milieu_arete(1, 1) ; x2 = milieu_arete(1, 2)

! ----------------------------------------------------------------------------------------------
! Coefficient de diffusion uniforme :
!       Cas tests 1, 2, 3, 4, 5
!       Plaque chauffante periodique
! ----------------------------------------------------------------------------------------------
            if (1 <= probleme .AND. probleme <= 6) then
                DXe = 1._pr

! ----------------------------------------------------------------------------------------------
! Piece chauffee par un radiateur
! ----------------------------------------------------------------------------------------------
            else if (probleme == 7) then
! Epaisseur totale de la vitre
                e = ea + 2*ev

                if ((dm <= x1 .AND. x1 <= dm+eb) .AND. ((x2 <= hb) .OR. (1._pr - hb <= x2))) then
! On est sur du beton
                    DXe = 1.03_pr
                else if (((dm + eb - e <= x1 .AND. x1 <= dm + eb - e + ev) .OR. (dm + eb - ev <= x1  &
                &         .AND. x1 <= dm + eb)) .AND. ((hb <= x2) .AND. (x2 <= hb + hv))) then
! On est sur du verre
                    DXe = 0.93_pr
                else
! On est sur de l'air
                    DXe = 0.025_pr
                end if

! ----------------------------------------------------------------------------------------------
! Cas test 6
! ----------------------------------------------------------------------------------------------
            else if (probleme == 9) then
                Dxe = 1._pr + 0.5_pr*SIN(2*pi*x1)*COS(2*pi*x2)
            end if

        end function D


        function Tinit(probleme, milieu_maille) result (TinitGi)

! Entrees de la fonction
            integer, intent(in)                                     :: probleme
            real(kind = pr), dimension(:, :), intent(in)            :: milieu_maille

! Sortie de la fonction :
!       TinitGi(1:nb_mailles) : Temperature initiale au centre de la maille i
            real(kind = pr), dimension(:), allocatable              :: TinitGi

! Variables locales
            integer                                                 :: nb_mailles, i
            real(kind = pr)                                         :: e, x1, x2, Tm

            nb_mailles = size(milieu_maille, 1)
            allocate(TinitGi(nb_mailles))

! ----------------------------------------------------------------------------------------------
! Temperature initiale uniforme :
!       Cas tests 1, 2, 3, 4, 5
!       Plaque chauffante periodique
! ----------------------------------------------------------------------------------------------
            if (1 <= probleme .AND. probleme <= 6) then
                TinitGi = 100._pr

! ----------------------------------------------------------------------------------------------
! Piece chauffee par un radiateur
! ----------------------------------------------------------------------------------------------
            else if (probleme == 7) then
                e = ea + 2*ev

                do i = 1, nb_mailles
                    x1 = milieu_maille(i, 1) ; x2 = milieu_maille(i, 2)
                    if (x1 <= dm) then
! On est a l'exterieur
                        TinitGi(i) = Text
                    else if (dm < x1 .AND. x1 <= dm + eb) then
! On est au niveau du mur, on a quatre cas :
!       On est dans le beton, on approxime par une temperature affine l'evolution de la temperature
! dans le bloc
!       On est a la hauteur de la vitre mais en dehors du double vitrage, on suppose que l'on est a
! la temperature exterieure
!       On est dans un des vitrages, on approxime par une temperature affine
!       On est dans l'air entre les vitrages, on suppose que la temperature reste contante
                        if (x2 < hb .OR. hb + hv < x2) then
! On est dans le beton
                            TinitGi(i) = (Tint - Text)/eb*(x1 - dm) + Text
                        else if (hb <= x2 .AND. x2 <= hb + hv) then
! Temperature moyenne entre l'exterieur et l'interieur, supposee presente dans l'air entre les vitrages
                            Tm = (Tint + Text)/2
                            if (x1 <= dm + eb - e) then
! On est encore dans l'air libre
                                TinitGi(i) = Text
                            else if (dm + eb - e < x1 .AND. x1 < dm + eb - e + ev) then
! On est sur le vitrage exterieur
                                TinitGi(i) = (Tm - Text)/ev*(x1 - (dm + eb - e)) + Text
                            else if (dm + eb - ev < x1 .AND. x1 < dm + eb) then
! On est sur le vitrage interieur
                                TinitGi(i) = (Tint - Tm)/ev*(x1 - (dm + eb - ev)) + Tm
                            else if (dm + eb - e + ev <= x1 .AND. x1 <= dm + eb - ev) then
! On est dans l'air entre les vitrages
                                TinitGi(i) = Tm
                            end if
                        end if
                    else if (dm + eb < x1) then
! On est a l'interieur
                        TinitGi(i) = Tint
                    end if
                end do
            
! ----------------------------------------------------------------------------------------------
! Plaque chauffante - Coefficient d'echange thermique
! ----------------------------------------------------------------------------------------------
            else if (probleme == 8) then
                TinitGi = Tinf
            
! ----------------------------------------------------------------------------------------------
! Cas test 6
! ----------------------------------------------------------------------------------------------
            else if (probleme == 9) then
                do i = 1, nb_mailles
                    x1 = milieu_maille(i, 1) ; x2 = milieu_maille(i, 2)
                    TinitGi(i) = Tinf + 50._pr*SIN(pi/2*x1)*COS(pi/2*x2)
                end do
            end if

        end function Tinit


        function Dirichlet(probleme, e, cl_arete_bord, t, milieu_arete) result (TbXe)

! Entrees de la fonction
            integer, intent(in)                                     :: probleme, e
            integer, dimension(:), intent(in)                       :: cl_arete_bord
            real(kind = pr), intent(in)                             :: t
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_arete 

! Sortie de la fonction :
!       TbXe : Temperature de bord evaluee au centre de l'arete e
            real(kind = pr)                                         :: TbXe

! Variables locales
            real(kind = pr)                                         :: x1, x2

! ----------------------------------------------------------------------------------------------
! Temperature de bord Tg et Td :
!       Cas tests 1, 2
! ----------------------------------------------------------------------------------------------
            if (probleme == 1 .OR. probleme == 2) then
                if (cl_arete_bord(e) == 10) then
! On est sur le bord gauche
                    TbXe = 100._pr
                else if (cl_arete_bord(e) == 11) then
! On est sur le bord droit
                    TbXe = 300._pr
                end if

! ----------------------------------------------------------------------------------------------
! Plaque chauffante circurlaire
! Solution manufacturee
! ----------------------------------------------------------------------------------------------
            else if (probleme == 3 .OR. probleme == 4 .OR. probleme == 6) then
                if (cl_arete_bord(e) == 10) then
! On est sur le bord gauche
                    TbXe = 100._pr
                else if (cl_arete_bord(e) == 11) then
! On est sur le bord droit, haut et bas
                    TbXe = 100._pr
                end if

! ----------------------------------------------------------------------------------------------
! Dirichlet pondere periodiquement :
!       Cas test 5
! ----------------------------------------------------------------------------------------------
            else if (probleme == 5) then
                if (cl_arete_bord(e) == 10) then
! On est sur le bord gauche
                    TbXe = 100._pr*(1._pr + 1._pr/5*SIN(2*pi*t))
                else if (cl_arete_bord(e) == 11) then
! On est sur le bord droit
                    TbXe = 300._pr*(1._pr + 1._pr/5*SIN(2*pi*t))
                end if

! ----------------------------------------------------------------------------------------------
! Piece chauffee par un radiateur
! ----------------------------------------------------------------------------------------------
            else if (probleme == 7) then
                if (cl_arete_bord(e) == 10) then
! On est sur le bord gauche
                    TbXe = Text
                else if (cl_arete_bord(e) == 11) then
! On est sur le bord droit
                    TbXe = Tint
                end if

! ----------------------------------------------------------------------------------------------
! Cas test 6
! ----------------------------------------------------------------------------------------------
            else if (probleme == 9) then
                x1 = milieu_arete(1, 1) ; x2 = milieu_arete(1, 2)
                if (cl_arete_bord(e) == 10) then
! On est sur le bord gauche
                    TbXe = 250._pr + 20._pr*COS(pi*x1*x2)*SIN(2*pi*t)
                else if (cl_arete_bord(e) == 11) then
! On est sur le bord droit
                    TbXe = 350._pr + 20._pr*SIN(pi*x1*x2)*COS(2*pi*t)
                end if
            end if

        end function Dirichlet


        function Neumann(probleme, e, cl_arete_bord, t, milieu_arete, Ti) result (PhibXe)

! Entrees de la fonction
            integer, intent(in)                                     :: probleme, e
            integer, dimension(:), intent(in)                       :: cl_arete_bord
            real(kind = pr), intent(in)                             :: t, Ti
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_arete 

! Sortie de la fonction :
!       PhibXe : Flux de bord evaluee au centre de l'arete e
            real(kind = pr)                                         :: PhibXe

! Variables locales
            real(kind = pr)                                         :: x1, x2

! ----------------------------------------------------------------------------------------------
! Neumann homogene :
!       Cas tests 1, 4
!       Piece chauffee par un radiateur
! ----------------------------------------------------------------------------------------------
            if (probleme == 1 .OR. probleme == 4 .OR. probleme == 7) then
                if (cl_arete_bord(e) == 20) then
! On est sur le bord haut ou bas
                    PhibXe = 0._pr
                end if

! ----------------------------------------------------------------------------------------------
! Neumann non homogene constant :
!       Cas test 2
! ----------------------------------------------------------------------------------------------
            else if (probleme == 2) then
                if (cl_arete_bord(e) == 20) then
! On est sur le bord haut ou bas
                    PhibXe = -1000._pr
                end if

! ----------------------------------------------------------------------------------------------
! Neumann non homogene pondere periodiquement :
!       Cas test 5
! ----------------------------------------------------------------------------------------------
            else if (probleme == 5) then
                if (cl_arete_bord(e) == 20) then
! On est sur le bord haut ou bas
                    PhibXe = -1000._pr*(1._pr + 1._pr/5*SIN(2*pi*t))
                end if

! ----------------------------------------------------------------------------------------------
! Coefficient d'echange :
!       Plaque chauffante circulaire
! ----------------------------------------------------------------------------------------------
            else if (probleme == 8) then
                if (cl_arete_bord(e) == 20) then
                    PhibXe = ht*(Ti - Tinf)
                end if

! ----------------------------------------------------------------------------------------------
! Cas test 6
! ----------------------------------------------------------------------------------------------
            else if (probleme == 9) then
                x1 = milieu_arete(1, 1) ; x2 = milieu_arete(1, 2)
                if (cl_arete_bord(e) == 20) then
                    PhibXe = 200._pr*SIN(2*pi*x1)*COS(2*pi*x2)*EXP(-0.1*t)
                end if
            end if

        end function Neumann


        function Heat_transfer_neumann(probleme, e, cl_arete_bord) result (PhibXe)

! Entrees de la fonction
            integer, intent(in)                                     :: probleme, e
            integer, dimension(:), intent(in)                       :: cl_arete_bord

! Sortie de la fonction :
!       PhibXe : Equivalent du flux a ajouter au second membre dans le cas du schema d'Euler
! Implicite
            real(kind = pr)                                         :: PhibXe

            if (probleme == 8) then
                if (cl_arete_bord(e) == 20) then
                    PhibXe = -ht*Tinf
                end if
            end if

        end function Heat_transfer_neumann


        function Terme_source(probleme, t, milieu_maille) result(SGi)

! Entrees de la fonction
            integer, intent(in)                                     :: probleme
            real(kind = pr), intent(in)                             :: t
            real(kind = pr), dimension(1, 2), intent(in)            :: milieu_maille
! Sortie de la fonction :
!       SGi : Terme source evalue au centre de la maille i
            real(kind = pr)                                         :: SGi

! Variables locales
            real(kind = pr)                                         :: x1, x2, r
            real(kind = pr), dimension(1:2)                         :: C

            if (probleme == 1 .OR. probleme == 2 .OR. probleme == 5) then
                SGi = 0._pr

! ----------------------------------------------------------------------------------------------
! Cas de la plaque chauffante circulaire
! ----------------------------------------------------------------------------------------------
            else if (probleme == 3 .OR. probleme == 8) then
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
            else if (probleme == 4) then
                SGi = 1000._pr

! ----------------------------------------------------------------------------------------------
! Cas de la plaque chauffante circulaire periodique
! ----------------------------------------------------------------------------------------------
            else if (probleme == 6) then
! Rayon de la plaque chauffante
                r = 0.25_pr
! Coordonees du centre de la plaque chauffante
                C = 0.5_pr

                x1 = milieu_maille(1, 1) ; x2 = milieu_maille(1, 2)
                if (((x1 - C(1))**2 + (x2 - C(2))**2) <= r**2) then
                    SGi = 1000._pr*(1._pr + 1._pr/5*SIN(2*pi*t)  + 0.05_pr * SIN(100*pi*t))
                else
                    SGi = 0._pr
                end if

! ----------------------------------------------------------------------------------------------
! Piece chauffee par un radiateur
! ----------------------------------------------------------------------------------------------
            else if (probleme == 7) then
                x1 = milieu_maille(1, 1) ; x2 = milieu_maille(1, 2)

                if ((dm + eb + drad <= x1 .AND. x1 <= dm + eb + drad + Lrad) .AND. (Hsol <= x2      &
                &    .AND. x2 <= Hrad)) then
! On est bien dans la zone du radiateur
                    if (t<=trad) then
! Le radiateur est en train d'emettre une source de chaleur
                        SGi = 35._pr
                    else
                        SGi = 0._pr
                    end if
                else
                    SGi = 0._pr
                end if

! ----------------------------------------------------------------------------------------------
! Cas test 6
! ----------------------------------------------------------------------------------------------
            else if (probleme == 9) then
                x1 = milieu_maille(1, 1) ; x2 = milieu_maille(1, 2)
                SGi = 1000._pr*SIN(2*pi*x1)*COS(2*pi*x2)*SIN(0.1*pi*t)
            end if

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