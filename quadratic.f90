program quadratic
    use, intrinsic :: iso_fortran_env, only : dp => real64

    implicit none

    real(dp) :: a, b, c
    real(dp) :: x1, x2
    real(dp) :: D

    real(dp), external :: discriminant

    print *, "Enter a, b, c"
    read *, a, b, c

    if( a .eq. 0.0_dp) stop "Fuck you"

    D = discriminant(a, b, c)
    print *, "Discriminant: D = ", D

    if (D .gt. 0.0) then
        call roots(a, b, c, x1, x2)
        print *, "Roots are ", x1, x2
    else if (D .eq. 0.0_dp) then
        call roots(a, b, c, x1, x2)
        print *, x1

    else
        print*, "No real roots"
    end if


end program quadratic

pure function discriminant(a, b, c) result(D)
    ! Computes the discriminant of a quadratic
    ! ax^2 + bx + c

    use, intrinsic :: iso_fortran_env, only : dp => real64

    implicit none

    real(dp), intent(in) :: a, b, c     ! coefficients
    real(dp) :: D                       ! discriminant

    D = b**2 - 4.0_dp * a * c

end function discriminant

subroutine roots(a, b, c, x1, x2)
    ! Computes the roots x1, x2 of a quadratic
    ! ax^2 + bx + c

    use, intrinsic :: iso_fortran_env, only : dp => real64

    implicit none

    real(dp), intent(in) :: a, b, c     ! coefficients
    real(dp), intent(out) :: x1, x2     ! roots
    real(dp) :: D                       ! discriminant

    real(dp), external :: discriminant            ! discriminant function

    ! degenerate quadratic
    if (a .eq. 0.0_dp) then
        print *, "Not a quadratic."
        return
    end if
    
    D = discriminant(a, b, c)

    if (D < 0.0_dp) then
        print "(A, F16.8)", "No real roots since D = ", D
        return
    end if
    
    D = sqrt(D)
    x1 = (-b + D) / (2.0_dp * a)
    x2 = (-b - D) / (2.0_dp * a)

end subroutine roots
