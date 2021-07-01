module otf
    use iso_fortran_env
    implicit none

    integer, private, parameter :: wp = real64

    real(kind=wp), public, parameter :: pi = 4.d0*datan(1.d0)

contains

    pure function ackley_2d(x, a, b, c) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a, b, c
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d 
        d = size(x, 2)

        fx = -a*exp(-b*sqrt((1.d0/d)*sum(x**2, 2)))&
             - exp((1.d0/d)*sum(cos(c*x), 2))&
             + a + exp(1.d0)
    end function ackley_2d

    pure function rosenbrock_2d(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2)-1) :: x1, x2

        integer :: d 
        d = size(x, 2)

        if (d < 2) &
            error stop "Rosenbrock function takes in two or more dimensions"

        x1 = x(:,:d-1)
        x2 = x(:,2:)

        fx = sum(100*(x2-x1**2)**2+(x1-1)**2, 2) + 1
    end function rosenbrock_2d

    pure function eggholder_2d(x) result(fx)
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
    end function eggholder_2d

    pure function cross_in_tray_2d(x) result(fx)
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
    end function cross_in_tray_2d

    pure function griewank_2d(x, a) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2)) :: seq

        integer :: d, i

        d = size(x, 2)

        seq = 1.d0
        do i=2, d
            seq(:,i) = seq(:,i)*i
        end do

        fx = sum(x**2/a, 2) &
             - product(cos( &
                x/sqrt(seq) &
             ), 2) + 2
    end function griewank_2d

    pure function levy_n13_2d(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: n

        n = size(x, 2)

        fx = sin(3*pi*x(:,1))**2 &
             + ((x(:,n)-1)**2)*(1+sin(2*pi*x(:,n))**2) &
             + sum(((x(:,:n-1)-1)**2)*(1+sin(3*pi*x(:,2:))**2), 2) &
             + 1
    end function levy_n13_2d

    pure function rastrigin_2d(x, a) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: n

        n = size(x, 2)

        fx = a*n + sum(x**2 - a*cos(2*pi*x), 2) + 1
    end function rastrigin_2d

! =============================================================================

    pure function ackley_1d(x, a, b, c) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: a, b, c
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = ackley_2d(reshape(x, [1, size(x)]), a, b, c)
        fx = fx_res(1)
    end function ackley_1d

    pure function rosenbrock_1d(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = rosenbrock_2d(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function rosenbrock_1d

    pure function eggholder_1d(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = eggholder_2d(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function eggholder_1d

    pure function cross_in_tray_1d(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = cross_in_tray_2d(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function cross_in_tray_1d

    pure function griewank_1d(x, a) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = griewank_2d(reshape(x, [1, size(x)]), a)
        fx = fx_res(1)
    end function griewank_1d

    pure function levy_n13_1d(x) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = levy_n13_2d(reshape(x, [1, size(x)]))
        fx = fx_res(1)
    end function levy_n13_1d

    pure function rastrigin_1d(x, a) result(fx)
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp), intent(in) :: a
        real(kind=wp) :: fx

        real(kind=wp), dimension(1:1) :: fx_res

        fx_res = rastrigin_2d(reshape(x, [1, size(x)]), a)
        fx = fx_res(1)
    end function rastrigin_1d

end module otf
