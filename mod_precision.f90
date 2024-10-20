module mod_precision

    integer, parameter          :: pr = 8
    integer, parameter          :: nb_max_sommets = 5

! Note : Il faut cfl <= 0.95 pour que le schema V.F converge sur mesh3_3
    real(kind = pr)             :: cfl = 0.95

end module mod_precision