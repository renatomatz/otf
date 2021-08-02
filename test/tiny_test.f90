program tiny_test

    use otf

    implicit none

    integer, parameter :: d = 2
    real(kind=wp), dimension(d) :: x
    real(kind=wp) :: fx

    !integer :: i
    !x = [(real(one/i, kind=wp), i=1, d)]

    x = -0.5_wp
    fx = perm0_single(x, 10.0_wp)

    print *, fx

end program tiny_test
