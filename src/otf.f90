module otf
    use iso_fortran_env
    implicit none

    integer, parameter :: wp = real64

    real(kind=wp), public, parameter :: pi = 4.0_wp*datan(1.0_wp)
    real(kind=wp), public, parameter :: one = 1.0_wp

contains

! =============================================================================

    pure function ackley_mult(x, a, b, c) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a, b, c
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d
        d = size(x, 2)

        fx = -a*exp(-b*sqrt((1.d0/d)*sum(x**2, 2)))&
             - exp((1.d0/d)*sum(cos(c*x), 2))&
             + a + exp(1.d0)
    end function ackley_mult

    pure function ackley_single(x, a, b, c) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: a, b, c
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = ackley_mult(reshape(x, [1, size(x)]), a, b, c)
        fx = fx_res(1)
    end function ackley_single

! =============================================================================

    pure function rosenbrock_mult(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2)-1) :: x1, x2

        integer :: d
        d = size(x, 2)

        if (d < 2) &
            error stop "Rosenbrock function takes in two or more dimensions"

        x1 = x(:,:d-1)
        x2 = x(:,2:)

        fx = sum(100*(x2-x1**2)**2+(1-x1)**2, 2) + 1
    end function rosenbrock_mult

    pure function rosenbrock_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = rosenbrock_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function rosenbrock_single

! =============================================================================

    pure function eggholder_mult(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1)) :: x1, x2

        integer :: d
        d = size(x, 2)

        if (d /= 2) &
            error stop "Eggholder function takes exactly two dimensions"

        x1 = x(:,1)
        x2 = x(:,2)

        fx = -(x2+47)*sin(sqrt(abs(x2+(x1/2)+47))) &
                  -x1*sin(sqrt(abs(x1-(x2+47))))
    end function eggholder_mult

    pure function eggholder_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = eggholder_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function eggholder_single

! =============================================================================

    pure function cross_in_tray_mult(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1)) :: x1, x2

        integer :: d
        d = size(x, 2)

        if (d /= 2) &
            error stop "Cross-in-Tray function takes exactly two dimensions"

        x1 = x(:,1)
        x2 = x(:,2)

        fx = -0.0001*(abs( &
            sin(x1)*sin(x2)*exp(abs( &
                100 - (sqrt(x1**2 + x2**2)/pi) &
            )) &
        )+1)**0.1
    end function cross_in_tray_mult

    pure function cross_in_tray_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = cross_in_tray_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function cross_in_tray_single

! =============================================================================

    pure function griewank_mult(x, a) result(fx)
        !! https://www.sfu.ca/~ssurjano/griewank.html

        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2)) :: seq

        integer :: n, d, i, j

        n = size(x, 1)
        d = size(x, 2)

        seq = reshape([((sqrt(real(j, kind=wp)), i=1, n), j=1, d)], [n, d])

        fx = sum(x**2/a, 2) &
             - product(cos( &
                x/sqrt(seq) &
             ), 2) + 2
    end function griewank_mult

    pure function griewank_single(x, a) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = griewank_mult(reshape(x, [1, size(x)]), a)
        fx = fx_res(1)
    end function griewank_single

! =============================================================================

    pure function levy_n13_mult(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d

        d = size(x, 2)

        fx = sin(3*pi*x(:,1))**2 &
             + ((x(:,d)-1)**2)*(1+sin(2*pi*x(:,d))**2) &
             + sum(((x(:,:d-1)-1)**2)*(1+sin(3*pi*x(:,2:))**2), 2) &
             + 1
    end function levy_n13_mult

    pure function levy_n13_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = levy_n13_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function levy_n13_single

! =============================================================================

    pure function rastrigin_mult(x, a) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d

        d = size(x, 2)

        fx = a*d + sum(x**2 - a*cos(2*pi*x), 2) + 1
    end function rastrigin_mult

    pure function rastrigin_single(x, a) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = rastrigin_mult(reshape(x, [1, size(x)]), a)
        fx = fx_res(1)
    end function rastrigin_single

! =============================================================================

    pure function levy_mult(x) result(fx)
        !! https://www.sfu.ca/~ssurjano/levy.html

        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d

        d = size(x, 2)

        fx = (sin(pi*w(x(:,1)))**2) &
             + sum((w(x(:,:d-1))-1)**2 &
                   * (one + 10*(sin(pi*w(x(:,:d-1))+one)**2)), 2) &
             + (((w(x(:,d))-one)**2)*(1+sin(2*pi*w(x(:,d)))**2))
    contains

        pure elemental function w(x) result (wx)
            real(kind=wp), intent(in) :: x
            real(kind=wp) :: wx

            wx = 1 + ((x-1)/4.0_wp)
        end function w


    end function levy_mult

    pure function levy_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = levy_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function levy_single

! =============================================================================

    pure function schwefel_mult(x) result(fx)
        !! https://www.sfu.ca/~ssurjano/schwef.html
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d

        d = size(x, 2)

        fx = (418.9829_wp*real(d, kind=wp)) - sum(x*sin(sqrt(abs(x))), 2)
    end function schwefel_mult

    pure function schwefel_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = schwefel_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function schwefel_single

! =============================================================================

    pure function perm0_mult(x, beta) result(fx)
        !! https://www.sfu.ca/~ssurjano/perm0db.html

        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: beta
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2), size(x, 2)) :: &
            seq, inv, poly_x

        integer :: n, d, i, j, k

        n = size(x, 1)
        d = size(x, 2)

        seq = reshape([(((k, i=1, n), j=1, d), k=1, d)], [n, d, d])
        inv = reshape([(((one/(k**j), i=1, n), j=1, d), k=1, d)], [n, d, d])
        do concurrent (j=1:d)
            poly_x(:,j,:) = x**j
        end do

        fx = sum(sum((seq+beta)*(poly_x-inv), 3), 2)

    end function perm0_mult

    pure function perm0_single(x, beta) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: beta
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = perm0_mult(reshape(x, [1, size(x)]), beta)
        fx = fx_res(1)
    end function perm0_single

! =============================================================================

    pure function perm_mult(x, beta) result(fx)
        !! https://www.sfu.ca/~ssurjano/permdb.html
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: beta
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2), size(x, 2)) :: &
            outer, inner, x3d

        integer :: n, d, i, j, k

        n = size(x, 1)
        d = size(x, 2)

        outer = reshape([(((j, i=1, n), j=1, d), k=1, d)], [n, d, d])
        inner = reshape([(((k, i=1, n), j=1, d), k=1, d)], [n, d, d])
        do concurrent (j=1:d)
            ! this j is the same as in the definition of inv
            x3d(:,j,:) = x
        end do

        fx = sum(sum((inner**outer+beta)*(((x3d/inner)**outer)-1), 3)**2, 2)
    end function perm_mult

    pure function perm_single(x, beta) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: beta
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = perm_mult(reshape(x, [1, size(x)]), beta)
        fx = fx_res(1)
    end function perm_single

! =============================================================================

    pure function trid_mult(x) result(fx)
        !! https://www.sfu.ca/~ssurjano/trid.html
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d

        d = size(x, 2)

        fx = sum((x-1)**2, 2) - sum(x(:,2:)*x(:,:d-1), 2)
    end function trid_mult

    pure function trid_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = trid_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function trid_single

! =============================================================================

    pure function zakharov_mult(x) result(fx)
        !! https://www.sfu.ca/~ssurjano/zakharov.html
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2)) :: seq

        integer :: n, d, i, j

        n = size(x, 1)
        d = size(x, 2)

        seq = reshape([((j, i=1, n), j=1, d)], [n, d])

        fx = sum(x**2, 2) &
             + sum(0.5_wp*seq*x, 2)**2 &
             + sum(0.5_wp*seq*x, 2)**4
    end function zakharov_mult

    pure function zakharov_single(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = zakharov_mult(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function zakharov_single

! =============================================================================

end module otf
