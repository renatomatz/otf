program viz_tests

    use sobseq
    use ogpf
    use otf

    implicit none

    integer, parameter :: n_points = 2**12

    type(multi_dim_sobol_state) :: sobol
    type(gpf) :: gp

    real(kind=wp), dimension(n_points,2) :: points, x
    real(kind=wp), dimension(n_points,1) :: fx 

    integer :: i = 0

    sobol = multi_dim_sobol_state(2)
    call sobol%populate(points)

    call gp%xlabel("Dimension 1") 
    call gp%ylabel("Dimension 2")
    call gp%options("set style data points")

    mainloop: do

        x = points
        fx = 0

        print *
        print *, "Optimization Function Visualization"
        print *, "1: Auckley Function (a=20, b=0.2, 2pi)"
        print *, "2: Rosenbrock Function"
        print *, "3: Eggholder Function"
        print *, "4: Cross-in-Tray Function"
        print *, "5: Griewank Function"
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
        case (4)
            x = points*2*10 - 10
            fx(:,1) = cross_in_tray(x)
        case (5)
            x = points*2*6 - 6
            fx(:,1) = griewank(x)
        case (0)
            exit mainloop
        case default
            print *
            print *, "Invalid visualization number, try again"
            cycle
        end select

        print *, "Plotting Visualization ..."

        call gp%surf(x(:,1:1), x(:,2:2), fx, "with points ps 1.2")

        print *
        print *, "Press any button to continue..."
        read *
    end do mainloop

    call execute_command_line("rm ogpf_temp_script.gp")

end program viz_tests
