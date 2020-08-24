!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Description:                                                       
!!     The routine for solving the He atom with LDA-DFT               
!!     Solving the radial equation by diaagnoliz the hamiltian matrix 
!!     Hamiltian matrix is from finite differential method             
!!
!!  History:
!!     V0 2018/01/03  Guodong Yu                                                                    
!!                          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE kindinit
    IMPLICIT NONE 
    ! Real value precision 
    INTEGER, PARAMETER :: SGL = 4 
    INTEGER, PARAMETER :: DBL = 8 
    ! Integer value precision
    INTEGER, PARAMETER :: SHORT = 4
    INTEGER, PARAMETER :: LONG = 8 
END MODULE kindinit

module init
    use kindinit
    implicit none
    REAL(KIND=DBL), PARAMETER :: PI = 3.14159265358979323846
    REAL(KIND=DBL), PARAMETER :: r0 = 0.0          ! the origin point
    REAL(KIND=DBL), PARAMETER :: RMAX = 100.0      ! the maximum radius
    INTEGER(KIND=SHORT), PARAMETER :: NUME = 2     ! the electron number 
    INTEGER(KIND=SHORT), PARAMETER :: Nlevl = 1    ! the calculated energy level number
    INTEGER(KIND=SHORT), PARAMETER :: Z = 2        ! The neuclei charge
    INTEGER(KIND=SHORT), PARAMETER :: NN = 8000   ! The radial points number
end module init

SUBROUTINE get_rl(rl)
    !
    ! Purpose:
    !   to sample the radial distance between r0 and rmax
    !
    ! input/output parameters
    !   rl  :  the sampling points [r0, rmax] (out, real(8), dimenstion(NN))
    !
    use init
    IMPLICIT NONE
    
    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: rl    

    ! local variables
    INTEGER(KIND=SHORT) :: i                          
    REAL(KIND=DBL) :: h
    h = ( RMAX - r0 ) / (NN - 1)        
    DO i = 1, NN  
        rl(i) = r0 + ( i - 1 ) * h
    END DO
END SUBROUTINE get_rl

SUBROUTINE kinetic(rl, T)
    !
    ! input/output parameters
    !    rl  :  the sampling points array (in, real*8, dimension(NN))
    !    T   :  the kinetic matrix (out, real*8, dimenstion(NN,NN))
    !
    ! Description:
    !   to get the differential matrix of kinetic operator
    !
    use init
    IMPLICIT NONE
    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl      
    REAL(KIND=DBL), DIMENSION(NN,NN), INTENT(OUT) :: T  

    ! local variables
    INTEGER(KIND=SHORT) :: i
    real(kind=dbl) :: h

    h = rl(2) - rl(1)            
    T = 0
    DO i = 2, NN - 1
        T(i, i-1) = -1.0 / ( 2.0 * h ** 2 )
        T(i, i+1) = -1.0 / ( 2.0 * h ** 2 )
        T(i, i) = 1.0 / h ** 2
    END DO
    T(1, 1) = 1.0 / h ** 2
    T(NN, NN) = 1.0 / h ** 2
    T(1, 2) = -1.0 / ( 2.0 * h ** 2 )
    T(NN, NN-1) = -1.0 / ( 2.0 * h ** 2 ) 
END SUBROUTINE kinetic
        
SUBROUTINE pot_nuclei(rl, vn)
    ! 
    ! input/output parameters
    !   rl : the sampling points array (in, real*8, dimension(NN))
    !   vn : the nuclei potential array (out, real*8, dimenstion(NN))
    ! 
    ! Purpose:
    !   To calculate the nuclei potential -Z/r
    !
    use init
    IMPLICIT NONE
    !arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl  
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: vn
    ! local variables
    INTEGER(KIND=SHORT) :: i                        
    DO i = 1, NN
        vn(i) = -Z / rl(i)
        !vn(i) = 0.5 * (rl(i) - 0.5*RMAX)**2 ! harmonic potential
    END DO
END SUBROUTINE pot_nuclei

SUBROUTINE pot_hartree(nl, rl, vH, eH)
    !
    ! input/output parameters
    !   nl  :   the charge density array (in, real*8, dimension(NN))   
    !   rl  :   the sampling points array (in, real*8, dimension(NN))
    !   vH  :   the hartree potential array (out, real*8, dimenstion(NN))
    !   eH  :   the hartree energy density array (out, real*8, dimenstion(NN))
    !           eH = 0.5 * vH
    ! Purpose:
    !   To get the Hartree potential by solving the Poisson equation
    !   namely, Laplacian vH = 4 * PI * n(r) 
    !
    ! After replce vH by U = r * vH
    ! U'' = -4 * pi * r * U , where U = r*vH, vH is Hartree potential
    ! Solve it using Numerov's method 
    ! Fit the equation to the form of x''(t) = f(t)x(t) + g(t)
    !
    use init
    USE ode
    IMPLICIT NONE
    
    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: nl     
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl     
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: vH    
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: eH    
    ! local variables
    REAL(KIND=DBL), allocatable :: fl(:), gl(:), ul(:)  
    REAL(KIND=DBL) :: U0 = 0.0                          
    REAL(KIND=DBL) :: dU0 = 1.0                        
    integer(kind=short) :: i
    
    allocate(fl(NN))
    allocate(gl(NN))
    allocate(ul(NN))
    fl = 0.0
    do i=1, NN
        gl(i) = -4.0 * pi * rl(i) * nl(i)
    end do
    CALL numerov2b(NN, rl, fl, gl, U0, dU0, Ul)
    Ul = Ul + ( NUME - Ul(NN)) / RMAX * rl
    do i = 1, NN
        vH(i) = Ul(i)/rl(i)
    end do
    eH = 0.5 * vH
    deallocate(fl)
    deallocate(gl)
    deallocate(ul)
END SUBROUTINE pot_hartree

SUBROUTINE pot_hartreeRL(nl, rl, vH, eH)
    !
    !Purpose:
    !   To get the Hartree potential by solving the Poisson equation
    !   namely, Laplacian vH = 4 * PI * n(r) 
    !
    ! After replce vH by U = r * vH
    ! U'' = -4 * pi * r * U , where U = r*vH, vH is Hartree potential
    ! Solve it using Numerov's method 
    ! Fit the equation to the form of x''(t) = f(t)x(t) + g(t)
    !
    use init
    USE ode
    IMPLICIT NONE
    
    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: nl    
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl   
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: vH    
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: eH    

    ! local variables
    REAL(KIND=DBL), allocatable :: fl(:), gl(:), ul(:), rl_tmp(:), ul_tmp(:)            
    REAL(KIND=DBL) :: U0 = NUME                          
    REAL(KIND=DBL) :: dU0 = -1.0                    
    integer(kind=short) :: i
   
    allocate(fl(NN), gl(NN), ul(NN), rl_tmp(NN), ul_tmp(NN))
    fl = 0.0
    do i=NN, 1, -1
        gl(i) = -4.0 * pi * rl(NN-i+1) * nl(NN-i+1)
        rl_tmp(i) = rl(NN-i+1)
    end do
    CALL numerov2b(NN, rl_tmp, fl, gl, U0, dU0, Ul)
    Ul = Ul - Ul(NN) + Ul(NN) / Rmax * rl_tmp
    do i = 1, NN
        Ul_tmp(i) = Ul(NN-i+1)
        vH(i) = Ul_tmp(i) / rl(i)
    end do
    eH = 0.5 * vH
    deallocate(fl,gl,ul, rl_tmp, ul_tmp)
END SUBROUTINE pot_hartreeRL

SUBROUTINE pot_xlda(nl, vx, ex)
    !    
    ! input/output parameters
    !   nl  :  the charge density array (in, real*8, dimension(NN))
    !   vx  :  the exchange potential under lda (out, real*8, dimension(NN))
    !   ex  :  the exchange energy density under lda (out, real*8, dimenstion(NN))
    !
    ! Purpose:
    !   to get the LDA exchange potential from charge density n 
    !
    use init
    IMPLICIT NONE

    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: nl    
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: vx   
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: ex  
    vx = -(3.0 * nl / PI) ** (1.0/3.0)
    ex = 3.0 / 4.0 * vx
END SUBROUTINE pot_xlda

SUBROUTINE pot_clda(nl, vc, ec)
    !
    ! input/output parameters
    !   nl  :  the charge density array (in, real*8, dimension(NN))
    !   vc  :  the correlation potential array (out, real*8, dimension(NN))
    !   ec  :  the correlation energy density array (out, real*8, dimension(NN))
    !
    ! Description:
    !   to get the LDA correlation potential from known charge density n 
    !
    use init
    IMPLICIT NONE

    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: nl    
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: vc   
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: ec 

    ! local variables
    REAL(KIND=DBL), allocatable ::  revx(:), atn(:), log1(:), log2(:)
    REAL(KIND=DBL) :: A, b, c, x0, Q

    allocate(revx(NN),atn(NN),log1(NN),log2(NN))
    A = 0.0621814
    b = 3.72744
    c = 12.9352
    x0 = -0.10498
    Q = SQRT(4*c - b**2)
    revx = ( 4.0 * PI * nl / 3.0 ) ** (1.0/6.0)

    atn = atan(Q*revx/(2.0+b*revx))
    log1 = -LOG(1.0 + b*revx + c* revx**2)
    log2 = log(1.0 + 2.0 * x0 * revx + x0**2 * revx**2) - log(1.0 + b*revx + c * revx**2)
    ec = A/2.0 * ( log1  + 2.0*b/Q * atn - b*x0/xx(x0) * ( log2 + 2.0*(b+2*x0)/Q * atn ) )

    vc = ec - (1.0/6.0) * A * ( (c-b*x0)*revx**2 - c*x0*revx**3 ) / &
                (1.0 + (b-x0)*revx + (c-b*x0)*revx**2 -c*x0*revx**3)

    deallocate(revx,atn,log1,log2) 
    contains 
        function xx(r)
            implicit none
            real(kind=dbl) :: r, xx
            xx = r**2 + b*r + c
        end function xx
END SUBROUTINE pot_clda

SUBROUTINE n_init( rl, n0l )
    !
    ! input/output parameters
    !   rl  :  the sampling points along radial direction (in, real*8, dimenstion(NN))
    !   n0l :  the  charge density (out, real*8, dimension(NN))
    !
    ! Description:
    !   to initial the normalized charge density function along radial direction
    !   according to n0(r) = c* exp(-a* (r-r0)**2 ) with a = 6.0, r0 = 1.0   
    !
    use init
    IMPLICIT NONE

    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: n0l        
    ! local variables
    REAL(KIND=SGL) :: a = 6.0, cent = 1.0    
    REAL(KIND=SGL) :: c, tmp, sum_tmp
    INTEGER(KIND=SHORT) :: i
    real(kind=dbl) :: h

    h = rl(2) - rl(1)
    sum_tmp = 0.0
    DO i = 1, NN
        tmp = 4.0 * PI * EXP( -a * ( rl(i) - cent )**2 ) * h * rl(i) ** 2
        sum_tmp = sum_tmp + tmp
    END DO
    c = NUME / sum_tmp
    n0l = c* EXP( -a * ( rl - cent )**2 )
END SUBROUTINE n_init

SUBROUTINE n_init1( rl, n0l )
    !  
    ! input/output parameters
    !   rl   :   the sampling points along r  (in, real*8, dimension(NN))
    !   n0l  :   the constant normalized charge density (out, real*8, dimension(NN))
    !
    !Purpose:
    !   to initial the normalized charge density function along radial direction
    !   n0(r) = c
    !
    use init
    IMPLICIT NONE

    ! arguements
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: n0l 
    ! local variables
    real(kind=dbl) :: com
    n0l = 1.0
    call integrate(n0l, com)
    n0l = 1.0 / com
END SUBROUTINE n_init1

SUBROUTINE hamiltonian(T, vH, vx, vc, vn, H)
    !
    ! input/output parameters
    !   T  :   the kinetic energy operator matrix  (in, real*8, dimension(NN,NN))
    !   vH :   the Hartree potential operator array (in, real*8, dimension(NN))
    !   vx :   the exchange potential operator array (in, real*8, dimension(NN))
    !   vc :   the correlation potential operator array (in, real*8, dimension(NN))
    !   vn :   the neuclei potential operator array (in, real*8, dimension(NN))
    !   H  :   the Hamiltnian matrix (out, real*8, dimension(NN,NN))
    !
    ! description
    !   To set up the hamiltonian matrix
    !
    use init
    IMPLICIT NONE

    ! arguements
    REAL(KIND=DBL), DIMENSION(NN, NN), INTENT(IN) :: T  
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: vH     
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: vx     
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: vc     
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: vn     
    REAL(KIND=DBL), DIMENSION(NN, NN), INTENT(OUT) :: H 
    ! local variables
    INTEGER(KIND=SHORT) :: i, j                         
    H = 0.0
    DO i = 1, NN
        Do j = 1, NN
            IF ( i == j ) THEN
                H(i, j) = T(i, j) + vn(i) + vH(i) + vx(i) + vc(i) 
            ELSE
                H(i, j) = T(i, j)
            END IF 
        END DO
    END DO
END SUBROUTINE hamiltonian

SUBROUTINE integrate(fl, com)
    ! 
    ! input/output parameters
    !   fl  :  the array along radial direction, such as charge density (in, real*8, dimension(NN))
    !   com :  the integration of fl (out, real*8)
    !
    ! Description:
    !   to integrate the fl function within the sphere with radius of rmax
    !   according to fl(r)*4*pi*r^2*sin(theta)dr*dtheta*dphi
    !
    USE INIT
    IMPLICIT NONE
    ! arguements
    real(kind = dbl), dimension(nn), intent(in) :: fl
    real(kind = dbl), intent(out) :: com
    ! local variables
    real(kind = dbl), allocatable :: rl(:)
    real(kind = dbl) :: tmp_sum, tmp
    integer(kind = short) :: i
    real(kind=dbl) :: h
    
    allocate(rl(NN))
    CALL get_rl(rl)
    h = rl(2) -rl(1)
    tmp_sum = 0.0
    DO i = 2, NN    
        tmp = 4.0 * PI * h * fl(i) * rl(i) **2
        tmp_sum = tmp_sum + tmp
    END DO
    com = tmp_sum
    deallocate(rl)
END SUBROUTINE integrate

SUBROUTINE diag_subdiag_parse(H, d, e, e2)
    !
    ! input/output parameters
    !   H  :  the hamiltonian matrix  (in, real*8, dimension(NN,NN))
    !   d  :  the diagnal elements array (out, real*8, dimension(NN))
    !   e  :  the subdiagnal elements array with e(1)=0.0 (out, real*8, dimension(NN))
    !   e2 :  the square of e (out, real*8, dimension(NN))
    !
    ! Description:
    !   to get the diagnal and subdiagnal elements preparing for solving the eigen problems
    !   using the subroutines imtqlv and tinvit in eispack library
    use init
    IMPLICIT NONE
    ! arguements
    REAL(KIND=DBL), DIMENSION(NN,NN), INTENT(IN) :: H   
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: d     
    REAL(KIND=DBL), DIMENSION(NN),INTENT(OUT) :: e  
    REAL(KIND=DBL), DIMENSION(NN),INTENT(OUT) :: e2 
    ! local variables                                      
    INTEGER(KIND=SHORT) :: i                          
    e = 0.0
    e2 = 0.0
    DO i = 1, NN
        d(i) = H(i,i)
        IF (i > 1) THEN
            e(i) = H(i-1, i) 
            e2(i) = e(i) ** 2
        END IF
    END DO
END SUBROUTINE diag_subdiag_parse

SUBROUTINE eigensolver(Hamiltn, eigs, vecs, ierrval, ierrvec)
    !
    ! input/output parameters
    !   Hamiltn  :   the hamiltion matrix (in, real*8, dimension(NN,NN))
    !   eigs     :   the calculated eigenvalues array (out, real*8, dimension(NN))
    !   vecs     :   the calculated eigenvectors array (out, real*8, dimension(NN-1,nelevl))
    !   ierrval  :   the error state when getting eigen values (out, integer*4)
    !   ierrvec  :   the error state when getting eigen vectors (out, integer*4)
    !
    ! Description:
    !   To do the eigen problems of Hamiltian. Because of the divergence of the neuclei potential
    !   at origin, so the first matrix element of Hamilton will be removed when solving the eigen
    !   problems. And the vecs have the dimension of (NN-1, nelevl)
    !
    use init
    IMPLICIT NONE

    ! Arguements
    REAL(KIND=DBL), DIMENSION(NN,NN), INTENT(IN) :: Hamiltn        
    REAL(KIND=DBL), DIMENSION(NN-1,nlevl), INTENT(OUT) :: vecs 
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: eigs        
    INTEGER(KIND=SHORT), INTENT(OUT) :: ierrval               
    INTEGER(KIND=SHORT), INTENT(OUT) :: ierrvec               
    ! local variables
    REAL(KIND=DBL), allocatable :: d(:), e(:), e2(:), rl(:),ind(:)        
    REAL(KIND=DBL) :: h

    allocate(d(NN),e(NN),e2(NN),rl(NN),ind(NN-1))
    CALL diag_subdiag_parse(Hamiltn, d, e, e2) 
    CALL imtqlv( NN-1, d(2:NN), e(2:NN), e2(2:NN), eigs, ind, ierrval )
    CALL tinvit( NN-1, d(2:NN), e(2:NN), e2(2:NN), Nlevl, eigs(1:Nlevl), ind(1:Nlevl), vecs, ierrvec)
    deallocate(d,e,e2,rl,ind)
END SUBROUTINE eigensolver

SUBROUTINE density_get(rl, vecs, n)
    !
    ! input/output parameters
    !   rl   :  the sampling points during [r0, rmax] (in, real*8, dimension(NN))
    !   vecs :  the occupied eigenvectors of hamilton (in, real*8, dimension(NN, nelevl))
    !   n    :  the normalized charge density array (out, real*8, dimension(NN))
    !  
    ! Description:
    !   To get the charge density distribution from the occupied eigen vectors 
    !
    use init
    IMPLICIT NONE

    ! arguements
    REAL(KIND=DBL), DIMENSION(NN-1, nlevl), INTENT(IN) :: vecs  
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: n            
    REAL(KIND=DBL), DIMENSION(NN), INTENT(IN) :: rl 
    ! local variables
    REAL(KIND=DBL) :: tmp, sum_tmp, com
    INTEGER(KIND=SHORT) :: i, j
    real(kind=dbl) :: h
    n(1) = 0.0_dbl
    h = rl(2) - rl(1)
    DO i = 2, NN
        sum_tmp = 0.0
        DO j = 1, nlevl
            tmp =  ( vecs(i-1, j) /  rl(i) ) ** 2 / (4.0 * PI * h)
            sum_tmp = sum_tmp +  tmp * 2.0
        END DO
        n(i) = sum_tmp
    END DO
END SUBROUTINE density_get

subroutine total_energy_from_n(rl, n, etot)
    !
    ! input/output parameters
    !   rl   :   the sampling points along radial direction  (in, real*8, dimension(NN))
    !   n    :   the charge density array (in, real*8, dimension(NN))
    !   etot :   the total energy (out, real*8)
    !
    ! Description:
    !   to calculate the total energy 
    !
    use init
    implicit none
    ! arguements
    real(kind=dbl), dimension(NN), intent(in) :: rl, n
    real(kind=dbl), intent(out) :: etot
    ! local variable
    real(kind=dbl), allocatable :: vn(:), T(:,:), vx(:), ex(:), vH(:), eH(:), vc(:), ec(:)
    real(kind=dbl), allocatable :: eigs(:), vecs(:,:), Hamiltn(:,:)
    real(kind=dbl):: sum_tmp, tmp, h
    integer(kind=short) :: i, ierrval, ierrvec

    allocate(vn(NN),T(NN,NN),Hamiltn(NN,NN),vx(NN),ex(NN),vH(NN),eH(NN),vc(NN),ec(NN))
    allocate(eigs(NN),vecs(NN-1,1))
    CALL pot_nuclei(rl, vn)  ! nuclei potential
    CALL kinetic(rl, T)      ! kinetic energy
    CALL pot_xlda(n, vx, ex)        ! exchange potential
    CALL pot_hartree(n, rl, vH, eH) ! Hartree potential
    CALL pot_clda(n, vc, ec)        ! the correlation potential
    CALL hamiltonian(T, vH, vx, vc, vn, Hamiltn) ! The hamiltonian
    CALL eigensolver(Hamiltn, eigs, vecs, ierrval, ierrvec)
    h = rl(2) - rl(1)
    sum_tmp = 0.0
    DO i = 2, NN
        tmp = 4.0 * PI * rl(i) ** 2 * n(i) * ( eH(i) + ex(i) + ec(i) - vH(i) - vx(i) - vc(i))
        sum_tmp = sum_tmp + tmp
    END DO            
    etot = 2.0 * eigs(1) + sum_tmp * h
    deallocate(vn,T,Hamiltn,vx,ex,vh,eh,vc,ec)
    deallocate(eigs)
end subroutine total_energy_from_n


SUBROUTINE SCFsolver(echeck, mix, tol, nmax, mix_start, Hout, nfinal)
    !    
    ! Input/output parameters 
    !   echeck    :   whether the total energy difference is taken for scf convergence check (logical, in)
    !   mix       :   the mixture percent of charge density of the last step (real*8, in)
    !   tol       :   the tolarance for convergence check (real*8, in)
    !   nmax      :   the maximum steps for scf, routine will stop after nmax steps (integer*4, in)
    !   mix_start :   after mix_start step, mixture method will take over 
    !   Hout      :   the converged hamiltonian (real*8, dimension(NN,NN), out)
    !   nfinal    :   the converged charge density (real*8, dimension(NN), out)
    !
    ! Description:
    !   To solve the LDA-DFT ks single-electron equation self-consistantly
    !
    use init
    IMPLICIT NONE
    !DATA dictionary
    LOGICAL, INTENT(IN) :: echeck           
    REAL(KIND=DBL), INTENT(IN) :: mix                 
    REAL(KIND=DBL), INTENT(IN) :: tol                 
    INTEGER(KIND=SHORT), INTENT(IN) :: nmax         
    INTEGER(KIND=SHORT), INTENT(IN) :: mix_start      
    REAL(KIND=DBL), DIMENSION(NN), INTENT(OUT) :: nfinal  
    REAL(KIND=DBL), DIMENSION(NN, NN), INTENT(OUT):: Hout 
    
    REAL(KIND=DBL), allocatable:: vH(:), eH(:), vn(:),vx(:),ex(:),vc(:), ec(:), vecs(:,:)   
    REAL(KIND=DBL), allocatable:: T(:,:), Hamiltn(:,:) 
    REAL(KIND=DBL), allocatable :: nin(:), nout(:)  
    REAL(KIND=DBL), allocatable :: rl(:)  
    INTEGER(KIND=SHORT) :: istep          
    REAL(KIND=DBL), allocatable :: eigs(:)        
    INTEGER(KIND=SHORT) :: ierrval, ierrvec           
    REAL(KIND=DBL) :: diff                 
    REAL(KIND=DBL) :: etot_in, etot_out ,etot         
    
    allocate(vH(NN),eH(NN),vn(NN),vx(NN),ex(NN),vc(NN),ec(NN))
    allocate(T(NN,NN),Hamiltn(NN,NN),vecs(NN-1,Nlevl))
    allocate(nin(NN),nout(NN), rl(NN), eigs(NN))
        

    CALL get_rl(rl)                
    CALL n_init(rl, nin)

    CALL pot_nuclei(rl, vn)  
    CALL kinetic(rl, T)      
    CALL pot_xlda(nin, vx, ex)        
    CALL pot_hartree(nin, rl, vH, eH) 
    CALL pot_clda(nin, vc, ec)        
    CALL hamiltonian(T, vH, vx, vc, vn, Hamiltn) 
    CALL eigensolver(Hamiltn, eigs, vecs, ierrval, ierrvec)
    CALL density_get(rl, vecs, nout)
    istep = 0
    IF ( echeck ) THEN
        CALL total_energy_from_n(rl, nout, etot_in)
        print *, etot_in
    END IF
    diff = 1.0E+20
    deallocate(vH,eH,vx,ex,vc,ec,vecs,eigs)
    DO WHILE ( diff > tol )
        allocate(vH(NN),eH(NN),vx(NN),ex(NN),vc(NN),ec(NN),vecs(NN,Nlevl),eigs(NN))
        istep = istep + 1
        IF (istep >= mix_start) THEN 
            nin = mix * nin + ( 1 - mix ) * nout  
        ELSE
            nin = nout
        END IF
        CALL pot_xlda(nin, vx, ex)
        CALL pot_hartree(nin, rl, vH, eH)
        CALL pot_clda(nin, vc, ec)
        CALL hamiltonian(T, vH, vx, vc, vn, Hamiltn)
        CALL eigensolver(Hamiltn, eigs, vecs, ierrval, ierrvec)
        CALL density_get(rl, vecs, nout)
        IF ( echeck ) THEN
            CALL total_energy_from_n(rl, nout, etot_out)
            print *, etot_out
            diff = ABS(etot_out - etot_in)
            etot_in = etot_out
        ELSE
            diff = MAXVAL(ABS(nout - nin))
        END IF

        IF ( istep >= nmax ) THEN
            EXIT
        END IF
        PRINT *, "step: ", istep, "diff: ", diff 
        !write(*,*) eigs
        deallocate(vH,eH,vx,ex,vc,ec,vecs,eigs)
    END DO
    PRINT *, "converged"
    nfinal = nout
    Hout = Hamiltn
    print *, 'Total energy: ', etot_out, 'Hartree' 
    deallocate(T,Hamiltn)
    deallocate(nin,nout, rl)
END SUBROUTINE SCFsolver

PROGRAM  main
    USE INIT
    IMPLICIT NONE
    REAL(KIND=DBL) :: mix=0.91, tol=0.000001
    REAL(KIND=DBL) :: etot
    INTEGER(KIND=SHORT) :: nmax=100, mix_start=300
    REAL(KIND=DBL), allocatable :: Hamiltn(:,:), n(:), rl(:)   ! converged Hamiltonian
    
    allocate(Hamiltn(NN,NN),n(NN), rl(NN))
    CALL get_rl(rl)
    CALL SCFsolver(.TRUE., mix, tol, nmax, mix_start, Hamiltn, n)
    CALL total_energy_from_n(rl, n, etot)
    deallocate(Hamiltn, n, rl)
END PROGRAM main


