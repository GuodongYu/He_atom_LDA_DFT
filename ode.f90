MODULE ode
implicit none
CONTAINS
    SUBROUTINE numerov1a(f, u, h, t, t0, x0, x1, xt)
        !
        !Description:
        !   to integrate the ordinary differential equation of d^2x/dt^2 = f(t)*x(t)+u(t) 
        !           with x(t0)=x0, x(t0+h)=x1. output x(t).
        !   
        !   Error: order 6
        ! 
        !input/ourput parameters:
        !   f    :   the external function of f(t) (real*8, external, in)
        !   u    :   the external function of u(t) (real*8, external, in)
        !   h    :   the step length (real*8, in)
        !   t    :   the solution at t time will be returned (real*8, in)
        !   t0   :   the beginning point for integration (real*8, in)
        !   x0   :   the initial boundary condition x(t0) (real*8, in)
        !   x1   :   the initial boundary condition x(t0+h) (real*8, in)
        !   xt   :   the solution x at x(t) (real*8, out)
        !
        IMPLICIT NONE 
        
        ! arguements
        REAL(KIND=8), INTENT(IN)  :: h     
        REAL(KIND=8), INTENT(IN)  :: t0    
        REAL(KIND=8), INTENT(IN)  :: x0    
        REAL(KIND=8), INTENT(IN)  :: x1    
        REAL(KIND=8), EXTERNAL    :: f, u  
        REAL(KIND=8), INTENT(IN)  :: t    
        REAL(KIND=8), INTENT(OUT) :: xt    
        ! local variables
        REAL(KIND=8) :: tt, xtt_h, wtt_2h, wtt_h, wtt 
        INTEGER(KIND=4) :: nstep 
        REAL(KIND=8) :: c 
        c = h**2 / 12.0
        nstep = int(t/h) 
        wtt_2h = ( 1.0 - c * f(t0) ) * x0  - c * u(t0)
        wtt_h = ( 1.0 - c * f(t0+h) ) * x1 - c * u(t0+h) 
        tt = t0 + h + h
        do while (tt <= nstep*h )
            xtt_h = ( wtt_h + c * u(tt-h) ) / ( 1 - c * f(tt-h) ) 
            wtt =  2.0 * wtt_h - wtt_2h + ( h**2 * ( f(tt-h) * xtt_h + u(tt-h) ) )
            wtt_2h = wtt_h
            wtt_h = wtt
            tt  = tt + h
        end do 
        xt = ( wtt + c * u(tt) ) / ( 1 - c * f(tt) )
    END SUBROUTINE numerov1a

    SUBROUTINE numerov1aa(f, u, h, t, t0, x0, x1, xt)
        !
        !Description:
        !   Same as subroutine numerov1a but saving the variables as array          
        !   Error: order 6
        ! 
        !input/output parameters:
        !   f    :   the external function of f(t) (real*8, external, in)
        !   u    :   the external function of u(t) (real*8, external, in)
        !   h    :   the step length (real*8, in)
        !   t    :   the solution at t time will be returned (real*8, in)
        !   t0   :   the beginning point of integration (real*8, in)
        !   x0   :   the initial condition x(t0) (real*8, in)
        !   x1   :   the initial condition x(t0+h) (real*8, in)
        !   xt   :   the solution x at x(t) (real*8, out)
        !
        IMPLICIT NONE 
        
        !Data dictornary
        REAL(KIND=8), INTENT(IN)  :: h     ! The step length 
        REAL(KIND=8), INTENT(IN)  :: t0    ! The beginning point
        REAL(KIND=8), INTENT(IN)  :: x0    ! The value of x(t0) 
        REAL(KIND=8), INTENT(IN)  :: x1    ! The value of x(t0+h)
        REAL(KIND=8), EXTERNAL    :: f, u  ! The known functions of f(t) and u(t)
        REAL(KIND=8), INTENT(IN)  :: t     ! The variable t 
        REAL(KIND=8), INTENT(OUT) :: xt    ! the x value at t
        REAL(KIND=8) :: c 
        integer(kind=4) :: n, i
        real(kind=8), allocatable :: tl(:), xl(:), wl(:) 

        n = abs( (t-t0) / h )
        allocate(tl(n))
        allocate(xl(n))
        allocate(wl(n))
        c = h**2 / 12.0
        tl(1) = t0
        do i=2, n
            tl(i) = tl(i-1) + h
        end do
        xl(1) = x0
        xl(2) = x1
        wl(1) = xl(1) * ( 1 - c * f(tl(1)) ) - c * u(tl(1))
        wl(2) = xl(2) * ( 1 - c * f(tl(2)) ) - c * u(tl(2))
        do i = 3, n
            wl(i) = 2*wl(i-1) - wl(i-2) + h**2 * ( f(tl(i-1))*xl(i-1) + u(tl(i-1)) )
            xl(i) = (wl(i) + c*u(tl(i))) / ( 1 - c * f(tl(i)) )
        end do
        xt = xl(n)
        deallocate(tl)
        deallocate(xl)
        deallocate(wl)
    END SUBROUTINE numerov1aa

    SUBROUTINE numerov1b(f, u, h, t, t0, x0, dx0, xt)
        !
        !Description:
        !   to integrate the ordinary differential equation of d^2x/dt^2 = f(t)*x(t)+u(t) 
        !           with boundary conditions x(t0) = x0, dx/dt = dx0 at t0. Return x(t).
        !   Error: order 6
        !
        !input/output parameters 
        !   f    :   the external function of f(t) (real*8, external, in)
        !   u    :   the external function of u(t) (real*8, external, in)
        !   h    :   the step length (real*8, in)
        !   t    :   the solution at t time will be returned (real*8, in)
        !   t0   :   the beginning point of integration (real*8, in)
        !   x0   :   the initial condition x(t0) (real*8, in)
        !   dx0  :   the initial condition dx/dt at t0 (real*8, in)
        !   xt   :   the solution x at x(t) (real*8, out)
        !
        IMPLICIT NONE
        
        !Data dictornary
        REAL(KIND=8), INTENT(IN)  :: h     
        REAL(KIND=8), INTENT(IN)  :: t0    
        REAL(KIND=8), INTENT(IN)  :: x0    
        REAL(KIND=8), INTENT(IN)  :: dx0   
        REAL(KIND=8), EXTERNAL    :: f, u  
        REAL(KIND=8), INTENT(IN)  :: t    
        REAL(KIND=8), INTENT(OUT) :: xt   
        REAL(KIND=8) :: x1                 

        x1 = x0 + h * dx0
        CALL numerov1a(f, u, h, t, t0, x0, x1, xt)
    END SUBROUTINE numerov1b

    SUBROUTINE numerov2a(n, tl, fl, ul, x0, x1, xl)
        !
        !Description:
        !   to integrate the ordinary differential equation of d^2x/dt^2 = f(t)*x(t)+u(t) 
        !   with boundary condition x(tl(1)) = x0, and x(tl(2)) = x1 
        !   Return solution x(tl) in xl
        !   Error: order 6
        !
        !input/output parameters
        !   n  :  the size of the integration points (ingeter*4, in)
        !   tl :  the t sampling array (real*8, dimension(n), in)
        !   fl :  the value array of f(tl) (real*8, dimension(n), in)
        !   ul :  the value array of u(tl) (real*8, dimension(n), in)
        !   x0 :  the boundary condition x(t0) value (real*8, in) 
        !   x1 :  the boundary condition x(t0+h) value (real*8, in) 
        !   xl :  the solution x as array (real*8, dimension(n), out)
        !
        IMPLICIT NONE
        
        ! arguments
        INTEGER :: n                            
        REAL(KIND=8), DIMENSION(n), INTENT(IN) :: tl    
        REAL(KIND=8), DIMENSION(n), INTENT(IN) :: fl    
        REAL(KIND=8), DIMENSION(n), INTENT(IN) :: ul    
        REAL(KIND=8), INTENT(IN)  :: x0                 
        REAL(KIND=8), INTENT(IN)  :: x1                 
        REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: xl   
        ! local variables
        REAL(KIND=8), allocatable :: wl(:)                
        REAL(KIND=8) :: h                               
        REAL(KIND=8) :: c 
        INTEGER(KIND=4) :: i  
        
        allocate(wl(n))
        h = tl(2) - tl(1)
        c = h**2 / 12.0
        xl(1) = x0
        xl(2) = x1
        wl(1) = xl(1) * ( 1.0 - c * fl(1) )  - c * ul(1)
        wl(2) = xl(2) * ( 1.0 - c * fl(2) )  - c * ul(2)
        do i = 3, n
            wl(i) = 2.0 * wl(i-1) - wl(i-2) + h**2 * ( fl(i-1) * xl(i-1) + ul(i-1) ) 
            xl(i) = ( wl(i) + c * ul(i)) / ( 1.0 - c * fl(i) )
        end do 
        deallocate(wl)
    END SUBROUTINE numerov2a

    SUBROUTINE numerov2b(n, tl, fl, ul, x0, dx0, xl)
        !
        !Description:
        !   to integrate the ordinary differential equation of d^2x/dt^2 = f(t)*x(t)+u(t) 
        !   with boundary condition x(t0)=x0 and x(t0+h)=x1. 
        !   Return a solution array xl saving x(tl).
        !   Error: order 6
        !
        !input/output parameters
        !   n  :  the size of the integration points (ingeter*4, in)
        !   tl :  the t sampling array (real*8, dimension(n), in)
        !   fl :  the value array of f(tl) (real*8, dimension(n), in)
        !   ul :  the value array of u(tl) (real*8, dimension(n), in)
        !   x0 :  the boundary condition x(t0) value (real*8, in) 
        !   x1 :  the boundary condition x(t0+h) value (real*8, in) 
        !   xl :  the solution x as array (real*8, dimension(n), out)
        !
        IMPLICIT NONE
        
        !Data dictornary
        INTEGER :: n                            
        REAL(KIND=8), DIMENSION(n), INTENT(IN) :: tl    
        REAL(KIND=8), DIMENSION(n), INTENT(IN) :: fl    
        REAL(KIND=8), DIMENSION(n), INTENT(IN) :: ul    
        REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: xl   
        REAL(KIND=8), INTENT(IN)  :: x0                 
        REAL(KIND=8), INTENT(IN)  :: dx0                
        REAL(KIND=8), DIMENSION(n) :: wl               
        REAL(KIND=8)  :: x1                            
        REAL(KIND=8) :: h                              
 
        h = tl(2) - tl(1)
        x1 = x0 + h * dx0
        CALL numerov2a(n, tl, fl, ul, x0, x1, xl)
    END SUBROUTINE numerov2b

    SUBROUTINE numerov3a(n, f, u, t0, tend, x0, dx0, t, x)
        !
        !Purpose:
        !   to integrate the ordinary differential equation of d^2x/dt^2 = f(t)*x(t)+u(t) 
        !           if the x values at the first two points (x0, x1) are known.
        !           Here, the function of f(t) and u(t) are already known.
        !   Error: order 6
        !
        IMPLICIT NONE
        
        !Data dictornary
        INTEGER(kind=4), intent(in) :: n                            ! The size of all arrays
        REAL(KIND=8), external :: f, u    ! The array saving discrete f values
        REAL(KIND=8), INTENT(IN)  :: x0, dx0, t0, tend                 ! The value of x(t0) 
        REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: t, x ! the x values as an array
        real(kind=8), allocatable :: w(:)
        REAL(KIND=8) :: h                               ! The step length 
        REAL(KIND=8) :: c 
        INTEGER(KIND=4) :: i  

        
        allocate(w(n))
        do i = 1, n
            t(i) = t0 + (i - 1) * (tend - t0) / (n - 1)
        end do

        h = t(2) - t(1)
        c = h**2 / 12.0
        x(1) = x0
        x(2) = x0 + dx0 * h

        w(1) = x(1) * ( 1.0 - c * f(t(1)) )  - c * u(t(1))
        w(2) = x(2) * ( 1.0 - c * f(t(2)) )  - c * u(t(2))
        do i = 3, n
            w(i) = 2.0 * w(i-1) - w(i-2) + h**2 * ( f(t(i-1)) * x(i-1) + u(t(i-1)) ) 
            x(i) = ( w(i) + c * u(t(i))) / ( 1.0 - c * f(t(i)) )
        end do 
        deallocate(w)
    END SUBROUTINE numerov3a
END MODULE ode

