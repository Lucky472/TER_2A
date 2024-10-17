program test

    use precision
    use maillage

    implicit none

    integer                                         :: i

    integer, dimension(:), allocatable              :: S
    integer, dimension(:, :), allocatable           :: noeud_maille
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud

    call lecture_maillage('TYP2/mesh3_3.typ2', S, noeud_maille, coord_noeud)

    do i=1, size(noeud_maille, 1)
        print *, noeud_maille(i, :)
    end do
    
end program test
