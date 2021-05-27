program viz_tests

    use sobseq
    use ogpf
    use otf

    implicit none

    integer, parameter :: n_points = 1024

    type(multi_dim_sobol_state) :: sobol
    type(gpf) :: gp

    real(kind=wp), dimension(n_points,2) :: points, x
    real(kind=wp), dimension(n_points,1) :: fx 

    integer :: i = 0

    sobol = multi_dim_sobol_state(2)
    call sobol%md_populate(points)

    call gp%xlabel("Dimension 1") 
    call gp%ylabel("Dimension 2")
    call gp%options("set style data lines")

    mainloop: do

        x = points
        fx = 0

        print *
        print *, "Optimization Function Visualization"
        print *, "1: Auckley Function (a=20, b=0.2, 2pi)"
        print *
        print *, "Select function to visualize or 0 to exit"

        read (*,*) i

        select case(i)
        case (1)
            x = points*2*32.768 - 32.768
            fx(:,1) = auckley(x, 20.d0, 0.2d0, 2*pi)
        case (2)
            x = points*2*2.048 - 2.048
            fx(:,1) = rosenbrock(x)
        case (3)
            x = points*2*512 - 512
            fx(:,1) = eggholder(x)
        case (0)
            exit mainloop
        case default
            print *
            print *, "Invalid visualization number, try again"
            cycle
        end select

        print *, "Plotting Visualization ..."

        call gp%surf(points(:,1:1), points(:,2:2), fx, "with points", palette="jet")

        print *
        print *, "Press any button to continue..."
        read *
    end do mainloop

end program viz_tests
