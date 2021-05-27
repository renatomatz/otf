program viz_tests

    use sobseq
    use ogpf
    use otf

    implicit none

    integer, parameter :: n_points = 2**12 ! power must be a multiple of 2
    integer, parameter :: n_points_sqrt = int(sqrt(real(n_points)))

    type(multi_dim_sobol_state) :: sobol
    type(gpf) :: gp

    real(kind=wp), allocatable, dimension(:,:) :: temp_x, temp_y, temp_fx

    real(kind=wp), dimension(n_points,2) :: points, x
    real(kind=wp), dimension(n_points,1) :: fx 

    integer :: i = 0, pgen

    call gp%xlabel("Dimension 1") 
    call gp%ylabel("Dimension 2")

    pointsloop: do
        print *
        print *, "1 Sobol:"
        print *, "2 Linear:"
        print *
        print *, "Select Point Generator or 0 to Exit"
        read (*,*) i

        pgen = i
        select case (pgen)
        case (1)
            sobol = multi_dim_sobol_state(2)
            call sobol%populate(points)
            call gp%options("set style data points")
        case (2)

            allocate(temp_x(n_points_sqrt, n_points_sqrt))
            allocate(temp_y(n_points_sqrt, n_points_sqrt))
            allocate(temp_fx(n_points_sqrt, n_points_sqrt))

            call meshgrid(temp_x, temp_y, &
                         [(real(i, kind=wp), i=1, n_points_sqrt)])

            points(:,1) = reshape(temp_x, [n_points])
            points(:,2) = reshape(temp_y, [n_points])
            points = points/n_points_sqrt

            call gp%options("set style data lines")
            call gp%options("set hidden3d")
        case (0)
            print *, "Selected 0, Terminating Program"
            stop
        case default
            print *, "Invalid selection, try again"
            cycle pointsloop
        end select

        exit 
    end do pointsloop

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

        select case (pgen)
        case (1)
            call gp%surf(x(:,1:1), x(:,2:2), fx)
        case (2)
            temp_x = reshape(x(:,1), [n_points_sqrt, n_points_sqrt])
            temp_y = reshape(x(:,2), [n_points_sqrt, n_points_sqrt])
            temp_fx = reshape(fx, [n_points_sqrt, n_points_sqrt])
            call gp%surf(temp_x, temp_y, temp_fx, palette="jet")
        end select
    end do mainloop

    call execute_command_line("rm ogpf_temp_script.gp")

end program viz_tests
