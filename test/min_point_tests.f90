program min_point_tests

    use otf
    use cafut

    implicit none

    integer, parameter :: d = 10
    real(kind=wp), dimension(d) :: argmin_x

    type(TestSuite) :: ts
    real(kind=wp) :: fx, tgt
    integer :: i

    ts = TestSuite("Minimum Point Tests")

    argmin_x(1:2) = [512.0_wp, 404.2319_wp]
    tgt = -959.6407_wp
    fx = eggholder_single(argmin_x(1:2))
    call ts%add(TestRealVal("Eggholder", 0.0001_wp), tgt, fx)

    argmin_x = 0.0_wp
    tgt = 0.0_wp
    fx = ackley_single(argmin_x, 20.d0, 0.2d0, 2*pi)
    call ts%add(TestRealVal("Ackley"), fx, tgt)

    argmin_x = 1.0_wp
    tgt = 1.0_wp
    fx = rosenbrock_single(argmin_x)
    call ts%add(TestRealVal("Rosenbrock"), fx, tgt)

    argmin_x = 0.0_wp
    tgt = 1.0_wp
    fx = griewank_single(argmin_x, 200.0_wp)
    call ts%add(TestRealVal("Griewank"), fx, tgt)

    argmin_x = 1.0_wp
    tgt = 0.0_wp
    fx = levy_single(argmin_x)
    call ts%add(TestRealVal("Levy"), fx, tgt)

    argmin_x = 1.0_wp
    tgt = 1.0_wp
    fx = levy_n13_single(argmin_x)
    call ts%add(TestRealVal("Levy No. 13"), fx, tgt)

    argmin_x = 0.0_wp
    tgt = 1.0_wp
    fx = rastrigin_single(argmin_x, 10.0_wp)
    call ts%add(TestRealVal("Rastrigin"), fx, tgt)

    argmin_x = 420.9687_wp
    tgt = 0.0_wp
    fx = schwefel_single(argmin_x)
    call ts%add(TestRealVal("Schwefel", 0.001_wp), fx, tgt)

    argmin_x = [(real(one/i, kind=wp), i=1, d)]
    tgt = 0.0_wp
    fx = perm0_single(argmin_x, 2.0_wp)
    call ts%add(TestRealVal("Perm 0"), fx, tgt)

    argmin_x = [(real(i, kind=wp), i=1, d)]
    tgt = 0.0_wp
    fx = perm_single(argmin_x, 2.0_wp)
    call ts%add(TestRealVal("Perm"), fx, tgt)

    argmin_x = [(real(i*(d+one-i), kind=wp), i=1, d)]
    tgt = -(d*(d+4.0_wp)*(d-one))/6.0_wp
    fx = trid_single(argmin_x)
    call ts%add(TestRealVal("Trid"), fx, tgt)

    argmin_x = 0.0_wp
    tgt = 0.0_wp
    fx = zakharov_single(argmin_x)
    call ts%add(TestRealVal("Zakharov"), fx, tgt)

    call ts%runTests()

end program min_point_tests
