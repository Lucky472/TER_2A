program chaleur

    use mod_precision
    use mod_tri_maillage
    use mod_maillage

    implicit none

    integer                                         :: i, nb_mailles, nb_aretes

    integer, dimension(:), allocatable              :: sommets_maille
    real(kind = pr), dimension(:), allocatable      :: aire_maille, l_arete, d_arete
    integer, dimension(:, :), allocatable           :: noeud_maille, ar, trig
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud, milieu_arete


    call maillage('TYP2/mesh3_1.typ2', nb_mailles, nb_aretes, sommets_maille, noeud_maille, coord_noeud &
    &             , aire_maille, l_arete, d_arete, milieu_arete, ar, trig)

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
    
end program chaleur
