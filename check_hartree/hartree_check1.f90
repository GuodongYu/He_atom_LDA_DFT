program check
    use ode
    implicit none
    real(kind=8), parameter :: r00 = 0.0, rend = 10.0, U0 = 0.0, dU0 = 1.0
    integer(kind=4), parameter :: n = 5000
    real(kind=8), external :: f, u
    integer(kind=4) :: i    
    real(kind=8), allocatable :: fl(:), ul(:), rl(:), UHl(:)
    
    allocate(fl(n), ul(n), rl(n), uHl(n))

    do i = 1, n
        rl(i) = r00 + (i-1) * ( rend - r00 ) / ( n - 1 )
        fl(i) = 0.0
        ul(i) = u(rl(i))
    end do

    call numerov2b(n, rl, fl, ul, U0, dU0, UHl)
    UHl = UHl + ( 1.0 - UHl(n) ) / rl(n) * rl
    open(unit=15, file='UH_hydrogen1.dat')
    do i = 1, n
        write(15,*), Rl(i), UHl(i)
    end do
    deallocate(fl, ul, rl, uhl)
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
