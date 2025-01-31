program circle_area
    implicit none

    real :: pi
    real :: r
    
    pi = 3.141593
    r = 4.0

    print *, "Perimeter = ", 2.0 * pi * r
    print *, "Area =      ", pi * r**2
    
end program circle_area
