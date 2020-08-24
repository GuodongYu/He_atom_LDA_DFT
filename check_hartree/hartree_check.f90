program check
    use ode
    implicit none
    real(kind=8), parameter :: r00 = 0.0, rend = 10.0, U0 = 0.0, dU0 = 1.0
    integer(kind=4), parameter :: n = 5000
    real(kind=8), external :: f, u
    real(kind=8), dimension(n) :: UH, R
    integer(kind=4) :: i    

    call numerov3a(n, f, u, r00, rend, U0, dU0, R, UH)
    UH = UH + ( 1.0 - UH(n) ) / R(n) * R
    open(unit=15, file='UH_hydrogen.dat')
    do i = 1, n
        write(15,*), R(i), UH(i)
    end do
end program check

function f(x)
    implicit none
    real(kind=8) :: x, f
    f = 0.0
end function f

function u(x)
    implicit none
    real(kind=8) :: x, u
    u = -4.0 * x * exp( -2.0 * x )
end function u
