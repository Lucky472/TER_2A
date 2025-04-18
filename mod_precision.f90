module mod_precision

    integer, parameter          :: pr = 8
    integer, parameter          :: nb_max_sommets = 8
    integer, parameter          :: Nmax = 1e3
    real(kind = pr), parameter  :: pi = 4._pr*ATAN(1._pr)
    real(kind = pr), parameter  :: epsilon = 1e-6

! Note : Il faut cfl <= 0.95 pour que le schema V.F explicite converge sur mesh3_3
    real(kind = pr)             :: cfl = 0.95

end module mod_precision