module otf
    use iso_fortran_env
    implicit none

    integer, private, parameter :: wp = real64

    real(kind=wp), public, parameter :: pi = 4.d0*datan(1.d0)

contains

    pure function auckley(x, a, b, c) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), intent(in) :: a, b, c
        real(kind=wp), dimension(size(x, 1)) :: fx

        integer :: d 
        d = size(x, 2)

        fx = &
            -a*exp(-b*sqrt(&
                (1.d0/d)*sum(x**2, 2)))&
            -exp(&
                (1.d0/d)*sum(cos(c*x), 2))&
            +a+exp(1.d0)
    end function auckley

    pure function rosenbrock(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1), size(x, 2)-1) :: x1, x2

        integer :: d 
        d = size(x, 2)

        if (d < 2) &
            error stop "Rosenbrock function takes in two or more dimensions"

        x1 = x(:,:d-1)
        x2 = x(:,2:)

        fx = sum(100*(x2-x1**2)**2+(x1-1)**2, 2)
    end function rosenbrock

    pure function eggholder(x) result(fx)
        real(kind=wp), dimension(:,:), intent(in) :: x
        real(kind=wp), dimension(size(x, 1)) :: fx

        real(kind=wp), dimension(size(x, 1)) :: x1, x2

        integer :: d 
        d = size(x, 2)

        if (d /= 2) &
            error stop "Eggholder function takes in exactly two dimensions"

        fx = -(x2+47)*sin(sqrt(abs(x2+(x1/2)+47))) &
             -x1*sin(sqrt(abs(x1-(x2+47))))
    end function eggholder

end module otf
