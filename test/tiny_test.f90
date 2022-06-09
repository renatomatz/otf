program tiny_test

    use otf

    implicit none

    integer, parameter :: d = 4
    real(kind=wp), dimension(d) :: x
    real(kind=wp) :: fx

    !integer :: i
    !x = [(real(one/i, kind=wp), i=1, d)]

    x = -2.0_wp
    fx = perm_single(x, 0.5_wp)

    !print *, fx

end program tiny_test
