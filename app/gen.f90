program viz_tests

    use sobseq
    use otf

    implicit none

    integer, parameter :: n_points = 2**12 ! power must be a multiple of 2

    type(multi_dim_sobol_state) :: sobol

    real(kind=wp), allocatable, dimension(:,:) :: points, fx

    allocate (points(n_points,2), fx(n_points,3))

    sobol = multi_dim_sobol_state(2)
    call sobol%populate(points)

    open(unit=8, file="sobol.txt", status="new")
    write (8,'(2F19.15)') transpose(points)
    close(unit=8)

    fx = 0

    fx(:,1:2) = points*2*32.768 - 32.768
    fx(:,3) = ackley_mult(fx(:,1:2), 20.d0, 0.2d0, 2*pi)
    open(unit=8, file="ackley.txt", status="new")
    write (8,'(3F20.14)') transpose(fx)
    close(unit=8)

    fx(:,1:2) = points*2*2.048 - 2.048
    fx(:,3) = rosenbrock_mult(fx(:,1:2))
    open(unit=8, file="rosenbrock.txt", status="new")
    write (8,'(3F20.14)') transpose(fx)
    close(unit=8)

    fx(:,1:2) = points*2*512 - 512
    fx(:,3) = eggholder_mult(fx(:,1:2))
    open(unit=8, file="eggholder.txt", status="new")
    write (8,'(3F20.14)') transpose(fx)
    close(unit=8)

    fx(:,1:2) = points*2*10 - 10
    fx(:,3) = cross_in_tray_mult(fx(:,1:2))
    open(unit=8, file="cross_in_tray.txt", status="new")
    write (8,'(3F20.14)') transpose(fx)
    close(unit=8)

    fx(:,1:2) = points*2*6 - 6
    fx(:,3) = griewank_mult(fx(:,1:2), 200.0_wp)
    open(unit=8, file="griewank.txt", status="new")
    write (8,'(3F20.14)') transpose(fx)
    close(unit=8)

end program viz_tests
