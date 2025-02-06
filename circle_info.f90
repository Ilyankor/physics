program circle_info
    implicit none

    integer, parameter :: dp = selected_real_kind(15)
    integer, parameter :: n = 10    ! number of circles

    integer :: i    ! iterator 
    integer :: io   ! file io

    real(dp), dimension(n) :: radius
    real(dp) :: area, perimeter

    ! ask user for radius of n circles
    do i = 1, n
        write(*, "(A)", advance="no") "Enter radius of circle: "
        read *, radius(i)
        print "(A, I2, A, F16.8)", "i = ", i, ": radius(i) = ", radius(i)
    end do

    open(newunit=io, file="log.txt")

    ! compute circle info and write to file
    do i = 1, n
        call area_of_circle(radius(i), perimeter, area)
        write(io, "(I2, A, F16.8, A, F16.8, A, F16.8)") &
            i, ": radius = ", radius(i), ", area = ", area, ", perimeter = ", perimeter
    end do

    close(io)

end program circle_info

subroutine area_of_circle(radius, perimeter, area)
    implicit none

    integer, parameter :: dp = selected_real_kind(15)
    real(dp), parameter :: pi = 3.141592653589793_dp

    real(dp), intent(in) :: radius
    real(dp), intent(inout) :: perimeter, area

    perimeter = 2.0_dp * pi * radius
    area = pi * radius**2

end subroutine area_of_circle
