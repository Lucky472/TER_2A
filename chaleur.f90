program chaleur

    use precision
    use tri_maillage
    use maillage

    implicit none

    integer                                         :: i, nb_cotes

    integer, dimension(:), allocatable              :: sommets_maille
    integer, dimension(:, :), allocatable           :: noeud_maille, flag, e, ar, trig
    real(kind = pr), dimension(:, :), allocatable   :: coord_noeud

    call lecture_maillage('TYP2/test.typ2', sommets_maille, noeud_maille, coord_noeud)

    do i=1, size(noeud_maille, 1)
       print *, noeud_maille(i, :)
    end do

    call fill_flag(sommets_maille, noeud_maille, coord_noeud, nb_cotes, flag)

    print *, "Nombre d'aretes :", nb_cotes
    print *, "flag = "
    do i = 1, size(flag, 1)
        print *, flag(i, :)
    end do

    call fill_e(sommets_maille, noeud_maille, nb_cotes, flag, e)

    print *, "e = "
    do i = 1, size(e, 1)
        print *, e(i, :)
    end do

    call fill_ar_trig(sommets_maille, noeud_maille, nb_cotes, flag, ar, trig)

    print *, "ar = "
    do i = 1, size(ar, 1)
        print *, ar(i, :)
    end do
    
end program chaleur
