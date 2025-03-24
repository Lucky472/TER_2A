module mod_resolution

    use mod_precision
    use mod_caracteristiques
    use mod_maillage

! ----------------------------------------------------------------------------------------------
! Module contenant les subroutines pour mod_resolution :
!       make_A_matrix (temporaire, a des fins de verification)
!       make_A_COO
!       make_A_CSR
!       make_b
!       conjugate_gradient
! ----------------------------------------------------------------------------------------------

    implicit none

    contains

! Subroutine a des fins de verification : creee la matrice A sous forme pleine
        subroutine make_A_matrix(dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar,  &
        &                        trig, cl_arete_bord, A)

!
            integer, intent(in)                                         :: nb_mailles
            integer, dimension(:), intent(in)                           :: cl_arete_bord
            integer, dimension(:, :), intent(in)                        :: ar, trig
            real(kind = pr), intent(in)                                 :: dt
            real(kind = pr), dimension(:), intent(in)                   :: aire_maille, l_arete, d_arete
            real(kind = pr), dimension(:, :), intent(in)                :: milieu_arete

!
            real(kind = pr), dimension(:, :), allocatable, intent(out)  :: A

! Variables locales
            integer :: i, j, k
            integer, dimension(nb_max_sommets) :: e
            
            allocate(A(1:nb_mailles, 1:nb_mailles))
            
            A = 0._pr

            do i = 1, nb_mailles
                A(i, i) = aire_maille(i)
                e = ar(i, :)
                do j = 1, nb_max_sommets
                    if (e(j) /= 0) then
                        k = trig(e(j), 2)
                        if (k == 0 .AND. (cl_arete_bord(e(j)) == 10 .OR. cl_arete_bord(e(j)) == 11)) then
! On est sur une arete de bord avec une condition de Dirichlet
                            A(i, i) = A(i, i) + dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))
                        else if (k == 0 .AND. cl_arete_bord(e(j)) == 20) then
! On est sur une arete de bord avec une condition de Neumann
                            A(i, i) = A(i, i)
                        else if (k /= 0) then
! On est sur une arete interieure
                            A(i, i) = A(i, i) + dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))
                            if (k /= i) then
                                A(i, k) = A(i, k) - dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))
                                A(k, i) = A(i, k)
                            end if
                        end if
                    end if   
                end do
            end do
            
        end subroutine make_A_matrix
            

        
        subroutine make_A_COO(dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar,  &
        &                     trig, cl_arete_bord, A_row, A_col, A_val)

! Entrees de la subroutine
            integer, intent(in)                                         :: nb_mailles
            integer, dimension(:), intent(in)                           :: cl_arete_bord
            integer, dimension(:, :), intent(in)                        :: ar, trig
            real(kind = pr), intent(in)                                 :: dt
            real(kind = pr), dimension(:), intent(in)                   :: aire_maille, l_arete, d_arete
            real(kind = pr), dimension(:, :), intent(in)                :: milieu_arete

! Sorties de la subroutine :
!       A_row(1:nb_elements_non_nuls) : Tableau row_COO, ie l'indice i d'un element non nul de la matrice A
! (tries par colonne)
!       A_col(1:nb_elements_non_nuls) : Tableau col_COO, ie l'indice j d'un element non nul de la matrice A
! (tries par colonne)
!       A_val(1:nb_elements_non_nuls) : Tableau val_COO, ie la valeur d'un element non nul A(i,j) sachant que
! les termes extra-diagonaux sont d'abord ranges. Par exemple : (1,31) suivi de son symetrique (31,1) et la
! colonne est fermee par le terme diagonal (1,1)
            integer, dimension(:), allocatable, intent(out)             :: A_row, A_col
            real(kind = pr), dimension(:), allocatable, intent(out)     :: A_val

! Variables locales
            integer                                                     :: i, j, k
            integer, dimension(nb_max_sommets)                          :: e
            real(kind = pr)                                             :: Aii, Aik

            allocate(A_row(0), A_col(0), A_val(0))

            do i = 1, nb_mailles
                Aii = aire_maille(i)
                e = ar(i, :)
! N.B : On ne peut pas contruire la matrice A comme on le fait pour les flux, ie faire une boucle sur les
! aretes et ajouter (resp. soustraire) leur contribution aux deux mailles
! Ainsi, pour chaque maille, le tableau e recupere et stocke les aretes de celle-ci puis on ajoute la
! la contribution de chaque arete aux differents termes de A
                do j = 1, nb_max_sommets
                    Aik = 0._pr
                    if (e(j) /= 0) then
                        k = trig(e(j), 2)
                        if (k == 0 .AND. (cl_arete_bord(e(j)) == 10 .OR. cl_arete_bord(e(j)) == 11)) then
! On est sur une arete de bord avec une condition de Dirichlet
                            Aii = Aii + dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))
                        else if (k == 0 .AND. cl_arete_bord(e(j)) == 20) then
! On est sur une arete de bord avec une condition de Neumann
                            Aii = Aii
                        else if (k /= 0) then
! On est sur une arete interieure
                            Aii = Aii + dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))
                            if (k /= i) then
                                Aik = Aik - dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))
! Allocation pour A(i, k)
                                call push_back_int(A_row, i) ; call push_back_int(A_col, k)
                                call push_back_real(A_val, Aik)
! Allocation pour A(k, i)
                                call push_back_int(A_row, k) ; call push_back_int(A_col, i)
                                call push_back_real(A_val, Aik)
                            end if
                        end if
                    end if
                end do
! Allocation pour A(i, i)
                call push_back_int(A_row, i) ; call push_back_int(A_col, i)
                call push_back_real(A_val, Aii)
            end do

        end subroutine make_A_COO


        subroutine make_A_CSR(dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar,  &
        &                     trig, cl_arete_bord, row_CSR, col_CSR, val_CSR)

! Entrees de la subroutine
            integer, intent(in)                                         :: nb_mailles
            integer, dimension(:), intent(in)                           :: cl_arete_bord
            integer, dimension(:, :), intent(in)                        :: ar, trig
            real(kind = pr), intent(in)                                 :: dt
            real(kind = pr), dimension(:), intent(in)                   :: aire_maille, l_arete, d_arete
            real(kind = pr), dimension(:, :), intent(in)                :: milieu_arete
        
! Sorties de la subroutine :
!       row_CSR(1:nb_mailles+1) : Tableau row_CSR, ie row_CSR(i) = nombre d'elements non nuls de A sur
! les i precedentes colonnes
! N.B : row_CSR(1) = 0
!       col_CSR(1:nb_elements_non_nuls) : Tableau col_CSR, ie 1 <= col_CSR(j) = k <= nb_mailles l'indice
! de la ligne sur lequel se situe l'element non nul
! N.B : Pour chaque colonne i de A, on sait qu'il y a (row_CSR(i+1) - row_CSR(i)) elements non nuls
!       val_CSR(1:nb_elements_non_nuls) : Tableau val_CSR, ie le tableau des valeurs des elements non nuls
! de A
            integer, dimension(:), allocatable, intent(out)             :: row_CSR, col_CSR
            real(kind = pr), dimension(:), allocatable, intent(out)     :: val_CSR
        
! Variables locales
            integer                                                     :: i, j, n, unique_count
            integer, dimension(:), allocatable                          :: row_COO, col_COO, row_sorted, col_sorted
            real(kind = pr), dimension(:), allocatable                  :: val_COO, val_sorted

            call make_A_COO(dt, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, ar,    &
            &               trig, cl_arete_bord, row_COO, col_COO, val_COO)
        
            n = size(row_COO)
    
            allocate(row_sorted(1:n), col_sorted(1:n), val_sorted(1:n))
            row_sorted = row_COO ; col_sorted = col_COO ; val_sorted = val_COO

            deallocate(row_COO, col_COO, val_COO)
        
! Tri des indices (row, col)
            call sort_A_elements(n, row_sorted, col_sorted, val_sorted)
        
            allocate(row_CSR(1:nb_mailles + 1))
            allocate(col_CSR(1:n), val_CSR(1:n))

            row_CSR = 0 ; unique_count = 0
        
! Construction de col_CSR et val_CSR en évitant les doublons
            do i = 1, n
                if (i == 1 .OR. row_sorted(i) /= row_sorted(i-1) .OR. col_sorted(i) /= col_sorted(i-1)) then
                    unique_count = unique_count + 1
                    col_CSR(unique_count) = col_sorted(i)
                    val_CSR(unique_count) = val_sorted(i)
                else
! Somme des valeurs en cas de doublon
                    val_CSR(unique_count) = val_CSR(unique_count) + val_sorted(i)
                end if
            end do
        
! Construction du tableau row_CSR
            j = 1
            do i = 1, nb_mailles + 1
                do while (j <= unique_count .AND. row_sorted(j) == i-1)
                    j = j + 1
                end do
                row_CSR(i) = j - 1
            end do
        
            deallocate(row_sorted, col_sorted, val_sorted)
        
        end subroutine make_A_CSR


        subroutine make_b(dt, t, nb_mailles, aire_maille, l_arete, d_arete, milieu_arete, milieu_maille,    &
        &                 ar, trig, cl_arete_bord, Tn, b)

! Entrees de la subroutine
        integer, intent(in)                                         :: nb_mailles
        integer, dimension(:), intent(in)                           :: cl_arete_bord
        integer, dimension(:, :), intent(in)                        :: ar, trig
        real(kind = pr), intent(in)                                 :: dt, t
        real(kind = pr), dimension(:), intent(in)                   :: aire_maille, l_arete, d_arete, Tn
        real(kind = pr), dimension(:, :), intent(in)                :: milieu_arete, milieu_maille

! Sortie de la subroutine :
!       b(1:nb_mailles) : Second membre de la formulation matricielle du schema temporel
! d'Euler Implicite
        real(kind = pr), dimension(:), allocatable, intent(out)     :: b

! Variables locales
        integer                                                     :: i, j, k
        integer, dimension(nb_max_sommets)                          :: e
        real(kind = pr)                                             :: bi

        allocate(b(1:nb_mailles))

        do i = 1, nb_mailles
            bi = 0._pr
            e = ar(i, :)
! cf. subroutine make_A_COO pour la description du tableau e
            do j = 1, nb_max_sommets
                if (e(j) /= 0) then
                    k = trig(e(j), 2)
                    if (k == 0 .AND. (10 <= cl_arete_bord(e(j))) .AND. (cl_arete_bord(e(j)) <= 19)) then
! On est sur une arete de bord avec une condition de Dirichlet
                        bi = bi + dt*(l_arete(e(j))/d_arete(e(j)))*D(milieu_arete(e(j), :))*    &
                        &    Dirichlet(e(j), cl_arete_bord, t, milieu_arete)
                    else if (k == 0 .AND. (20 <= cl_arete_bord(e(j))) .AND. (cl_arete_bord(e(j)) <= 29)) then
! On est sur une arete de bord avec une condition de Neumann (ici Phih = Phib)
                        bi = bi - dt*l_arete(e(j))*Neumann(e(j), cl_arete_bord, t, milieu_arete)
                    end if
                end if
            end do
            b(i) = bi + aire_maille(i)*(Tn(i) + dt*Terme_source(milieu_maille(i, :)))
        end do

        end subroutine make_b


        subroutine conjugate_gradient(row_CSR, col_CSR, val_CSR, x0, bn, x)

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: row_CSR, col_CSR
            real(kind = pr), dimension(:), intent(in)                   :: val_CSR, x0, bn

! Sortie de la subroutine :
!       x : Resultat de l'algorithme
            real(kind = pr), dimension(1:size(x0)), intent(out)         :: x

! Variables locales
            integer                                                     :: k
            real(kind = pr)                                             :: rho0, rho, delta, alpha, gamma
            real(kind = pr), dimension(1:size(x0))                      :: r, p, q, result

            call CSR_dot_product(row_CSR, col_CSR, val_CSR, x0, result)
            x = x0 ; r = bn - result ; p = r
            k = 0
            rho0 = dot_product(r, r)

            do while (norm2(r)/norm2(bn) > epsilon .AND. k <= Nmax)
                call CSR_dot_product(row_CSR, col_CSR, val_CSR, p, q)
                delta = dot_product(p, q)

                if (delta < 1e-15) then
                    return
                end if
                alpha = rho0/delta

                x = x + alpha*p
                r = r - alpha*q
                rho = dot_product(r, r)

                gamma = rho/rho0
                p = gamma*p + r
                rho0 = rho
                k = k + 1
            end do

        end subroutine conjugate_gradient
  
! ----------------------------------------------------------------------------------------------
! Subroutines annexes implementant quelques methodes elementaires :
!       push_back_int
!       push_back_real
!       sort_A_elements
!       CSR_dot_product
! ----------------------------------------------------------------------------------------------

        subroutine push_back_int(vector, new_element)

! Entrees et sortie de la subroutine
            integer, intent(in)                                         :: new_element
            integer, dimension(:), allocatable, intent(inout)           :: vector

! Variables locales
            integer                                                     :: n
            integer, dimension(:), allocatable                          :: temp

            n = size(vector)

            allocate(temp(1:n+1))

            if (n > 0) then
                temp(1:n) = vector
            end if
            temp(n+1) = new_element

            deallocate(vector) ; allocate(vector(1:n+1))

            vector = temp
            deallocate(temp)
            
        end subroutine push_back_int


        subroutine push_back_real(vector, new_element)

! Entrees et sortie de la subroutine
            real(kind = pr), intent(in)                                 :: new_element
            real(kind = pr), dimension(:), allocatable, intent(inout)   :: vector

! Variables locales
            integer                                                     :: n
            real(kind = pr), dimension(:), allocatable                  :: temp

            n = size(vector)

            allocate(temp(1:n+1))

            if (n > 0) then
                temp(1:n) = vector
            end if
            temp(n+1) = new_element

            deallocate(vector) ; allocate(vector(1:n+1))

            vector = temp
            deallocate(temp)

        end subroutine push_back_real


        subroutine sort_A_elements(n, row_sorted, col_sorted, val_sorted)

            integer, intent(in)                                         :: n
            integer, dimension(:), intent(inout)                        :: row_sorted, col_sorted
            real(kind = pr), dimension(:), intent(inout)                :: val_sorted

! Variables locales
            integer                                                     :: i, j, temp_row, temp_col
            real(kind = pr)                                             :: temp_val
        
            do i = 1, n-1
                do j = i + 1, n
                    if (row_sorted(i) > row_sorted(j) .OR.                                          &
                    & (row_sorted(i) == row_sorted(j) .AND. col_sorted(i) > col_sorted(j))) then
! On prend avantage de la symétrie de la matrice et de la façon dont est stockée A
! au format COO (les éléments extra-diagonaux sont stockés directement en suivant)
                        temp_row = row_sorted(i) ; temp_col = col_sorted(i)
                        temp_val = val_sorted(i)

                        row_sorted(i) = row_sorted(j) ; col_sorted(i) = col_sorted(j)
                        val_sorted(i) = val_sorted(j)

                        row_sorted(j) = temp_row ;col_sorted(j) = temp_col
                        val_sorted(j) = temp_val

                    end if
                end do
            end do

        end subroutine sort_A_elements


        subroutine CSR_dot_product(row_CSR, col_CSR, val_CSR, x, result)

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: row_CSR, col_CSR
            real(kind = pr), dimension(:), intent(in)                   :: val_CSR, x

! Sortie de la subroutine :
!       result(1:size(x)) : Le resultat du produit matrice vecteur
            real(kind = pr), dimension(1:size(x)), intent(out)          :: result

! Variables locales
            integer                                                     :: i, j

            result = 0._pr

            if (size(row_CSR) - 1 /= size(x)) then
                print *, "Produit incorrect !"
                return
            else
                do i = 1, size(x)
                    do j = row_CSR(i), row_CSR(i+1)-1
                        result(i) = result(i) + x(col_CSR(j+1)) * val_CSR(j+1)
                    end do
                end do
            end if

        end subroutine CSR_dot_product

end module mod_resolution