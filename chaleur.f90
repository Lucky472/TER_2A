program chaleur

    use mod_precision
    use mod_tri_maillage
    use mod_maillage
    use mod_sortie

    implicit none

    integer                                         :: i, nb_mailles, nb_aretes

    integer, dimension(:), allocatable              :: sommets_maille, arete_bord
    real(kind = pr), dimension(:), allocatable      :: aire_maille, l_arete, d_arete, T
    integer, dimension(:, :), allocatable           :: noeud_maille, ar, trig
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud, milieu_arete


    call maillage('TYP2/mesh3_1.typ2', nb_mailles, nb_aretes, sommets_maille, noeud_maille, coord_noeud &
    &             , aire_maille, l_arete, d_arete, milieu_arete, ar, trig, arete_bord)

    ! print *, "aire_maille =" 
    ! do i = 1, size(aire_maille, 1)
    !     print *, aire_maille(i)
    ! end do

    ! print *, "l_arete =" 
    ! do i = 1, size(l_arete, 1)
    !     print *, l_arete(i)
    ! end do

    ! print *, "d_arete =" 
    ! do i = 1, size(d_arete, 1)
    !     print *, d_arete(i)
    ! end do

    ! print *, "milieu_arete =" 
    ! do i = 1, size(milieu_arete, 1)
    !     print *, milieu_arete(i, :)
    ! end do

    ! print *, "trig ="
    ! do i = 1, size(trig, 1)
    !     print *, trig(i, :)
    ! end do

    print *, "arete_bord ="
    do i = 1, size(arete_bord)
        print *, arete_bord(i)
    end do

    allocate(T(1:nb_mailles))
    T = 300._pr
    call sortie(0, T, sommets_maille, noeud_maille, coord_noeud)

    deallocate(sommets_maille, arete_bord, aire_maille, l_arete, d_arete &
    &          , noeud_maille, ar, trig, coord_noeud, milieu_arete, T)
    
end program chaleur
