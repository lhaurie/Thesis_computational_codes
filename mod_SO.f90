module parameters_so
    implicit none
    !number of solutions and number of momentum points for the average
    integer, parameter :: nsol=4 , Nav=300
    !Constants and temperature
    double precision,parameter ::Pi=4.d0*atan(1.d0),temp=0.001d0, beta=1/temp
    !electron density n, hopping parameters t, inter-orbital hopping parameter lambda, Coulomb interaction U
    double precision,parameter ::n=1.85d0, tx=1.d0, ty=1.d0,lambda=0.59d0, U=9.8d0
    !x0 are the initial guess, xx contains the self-consistent parameters. kx and ky runs over momentum.
    double precision :: xx(9),kx,ky, x01,x02,x03,x04,x05,x06,x07,x08, x09
end module parameters_so
    
    !***********************************************************************
    
module matrices_so
    use parameters_so
    implicit none
    real(8) :: m11x,m12x,m22x, m13, m14,m23,m24,m11y, m12y, m22y      !Elements of the M matrix
    double precision ::  al,alpha1_x, beta1_x, alpha1_y, beta1_y, nx, ny    !Elements of I^{-1} matrix, electron density
    real(8), dimension(int(nsol),int(nsol)) :: matM,matIinv,matI,matE  !Matrices definition
    integer :: i
    contains
    
    
    subroutine init_matrices()

      matM=0.d0
      matI=0.d0
      matIinv=0.d0
      matE=0.d0
      nx=(n+xx(2))/2.d0   !xx(2) defines nx and ny from delta n
      ny=(n-xx(2))/2.d0
      al=(cos(kx)+cos(ky))/2.d0   

      !Define the elements of I^{-1}
      alpha1_x=2.d0/(2.d0-nx)       
      beta1_x=2.d0/nx
      alpha1_y=2.d0/(2.d0-ny)
      beta1_y=2.d0/ny

      !Define the elements of the M matrix
      m11x=-xx(1)*(1.d0-nx/2.d0)- 4.d0*tx*xx(5)-4.d0*tx*al*(1.d0-nx+xx(3))-lambda*xx(7)
      m12x=4.d0*tx*xx(5) - 4.d0*tx*al*(nx/2.d0-xx(3))+lambda*xx(7)
      m22x=-(xx(1)-U)*nx/2.d0-4.d0*tx*xx(5)-4.d0*tx*al*xx(3)-lambda*xx(7)

      m11y=-xx(1)*(1.d0-ny/2.d0)- 4.d0*ty*xx(6)-4.d0*ty*al*(1.d0-ny+xx(4))-lambda*xx(8)
      m12y=4.d0*ty*xx(6) - 4.d0*ty*al*(ny/2.d0-xx(4))+lambda*xx(8)
      m22y=-(xx(1)-U)*ny/2.d0-4.d0*ty*xx(6)-4.d0*ty*al*xx(4)-lambda*xx(8)

      m13=-lambda*(1-nx/2.d0-ny/2.d0+xx(9))
      m14=-lambda*(ny/2.d0-xx(9))
      m23=-lambda*(nx/2.d0-xx(9))
      m24=-lambda*xx(9)

      !Define the matrices
      matI(1,1)=1.d0-nx/2.d0
      matI(2,2)=nx/2.0d0
      matI(3,3)=1.d0-ny/2.d0
      matI(4,4)=ny/2.d0

      matIinv(1,1)=alpha1_x
      matIinv(2,2)=beta1_x
      matIinv(3,3)=alpha1_y
      matIinv(4,4)=beta1_y

      matM(1,1)=m11x
      matM(1,2)=m12x
      matM(2,1)=m12x
      matM(2,2)=m22x

      matM(3,3)=m11y
      matM(3,4)=m12y
      matM(4,3)=m12y
      matM(4,4)=m22y

      matM(1,3)=m13
      matM(1,4)=m14
      matM(2,3)=m23
      matM(2,4)=m24
      ! M is symmetrical !
      matM(3,1)=m13
      matM(3,2)=m23
      matM(4,1)=m14
      matM(4,2)=m24

      matE=matmul(matM,matIinv)

    end subroutine init_matrices
end module matrices_so

!***************************************************************

module Cij_so
    use parameters_so
    use matrices_so
    implicit none
    real(8) :: C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44
    real(8) :: C011,C012,C013,C014,C021,C022,C023,C024,C031,C032,C033,C034,C041,C042,C043,C044
    contains

    subroutine coefficients(matCk1, matCk0)  
    !Returns the momentum dependent correlation functions

    real(8), dimension(nsol,nsol),intent(out) :: matCk1,matCk0
    integer :: i,k,j,lwork,lda,ldvl,ldvr,info
    double precision,dimension(nsol) :: wr, wi
    double precision, dimension(nsol,nsol) :: vl,vr
    complex(8) :: eigsk(nsol)
    double precision, allocatable :: work(:)
    double precision :: vr_inv(nsol,nsol)

    !Diagonalization of E, eigsk: eigenvalues, vr: right eigenvectors
    call init_matrices()
    lda=nsol
    ldvl=nsol
    ldvr=nsol
    lwork=10*nsol+10
    allocate(work(lwork))
      call dgeev('N','V',nsol,matE,nsol,wr, wi ,vl,ldvl,vr,ldvr,work,lwork,info)
      deallocate(work)
      do i=1,nsol
        eigsk(i)=wr(i)+complex(0,1)*wi(i)
      enddo

    !inverting the right eigenvector
      call inverse(vr,vr_inv,nsol)   
    !multiplying by I. vr(i,k) vr_inv(k,j) is the spectral weight (cf Anurag's proof)
    vr_inv=matmul(vr_inv,matI)
    matCk0=0.d0
    do i=1,nsol
      do j=1,nsol
        do k=1,nsol
          matCk0(i,j)=matCk0(i,j)+vr(i,k)*vr_inv(k,j)*1.0d0/2.d0*((1.d0+dtanh(beta*real(eigsk(k))/2.d0)))
        enddo
      enddo
    enddo
    !NN correlations
    matCk1=(al*matCk0)

    end subroutine coefficients




    subroutine averages(nstep)
    implicit none
    integer ::nstep,i,j
    real(8), dimension(nsol,nsol)::matCk1 ,matCk0
    real(8), dimension(nsol,nsol)::  outCk0, outCk1

    outCk0=0.d0
    outCk1=0.d0

    !Average of the correlation functions between 0 (included) and 2pi (excluded) for kx and ky
    do i=1,nstep
        do j=1,nstep
          kx=-pi+((2.0d0*pi*(i-1))/nstep)
          ky=-pi+((2.0d0*pi*(j-1))/nstep)
          call init_matrices()
          call coefficients(matCk1,matCk0)
          outCk1=outCk1 + matCk1  
          outCk0=outCk0 + matCk0
        enddo
    enddo
    outCk1=outCk1/(nstep**2)
    outCk0=outCk0/(nstep**2)

    !Associating explicitly the correlation functions
    C11=outCk1(1,1) 
    C12=outCk1(1,2)
    C13=outCk1(1,3)
    C14=outCk1(1,4)
    C21=outCk1(2,1)
    C22=outCk1(2,2)
    C23=outCk1(2,3)
    C24=outCk1(2,4)
    C31=outCk1(3,1)
    C32=outCk1(3,2)
    C33=outCk1(3,3)
    C34=outCk1(3,4)
    C41=outCk1(4,1)
    C42=outCk1(4,2)
    C43=outCk1(4,3)
    C44=outCk1(4,4)

    C011=outCk0(1,1)
    C012=outCk0(1,2)
    C013=outCk0(1,3)
    C014=outCk0(1,4)
    C021=outCk0(2,1)
    C022=outCk0(2,2)
    C023=outCk0(2,3)
    C024=outCk0(2,4)
    C031=outCk0(3,1)
    C032=outCk0(3,2)
    C033=outCk0(3,3)
    C034=outCk0(3,4)
    C041=outCk0(4,1)
    C042=outCk0(4,2)
    C043=outCk0(4,3)
    C044=outCk0(4,4)

    end subroutine averages
      
end module Cij_so

!----------------------------------------LLAPACK FUNCTIONS------------------------------------------------------------!


! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
subroutine zinverse(A,Ainv,n)
  integer :: n, info
  complex(8), dimension(n,n) :: A
  complex(8), dimension(n,n) :: Ainv
  complex(8), dimension(n) :: work  ! work array for LAPACK
  integer, dimension(n) :: ipiv   ! pivot indices


  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
    stop 'Matrix is numerically singular!'
  endif

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call ZGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
    stop 'Matrix inversion failed!'
  endif
end subroutine zinverse




! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
subroutine inverse(A,Ainv,n)
  integer :: n, info
  real(8), dimension(n,n) :: A
  real(8), dimension(n,n) :: Ainv
  real(8), dimension(n) :: work  ! work array for LAPACK
  integer, dimension(n) :: ipiv   ! pivot indices

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
    stop 'Matrix is numerically singular!'
  endif

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
    stop 'Matrix inversion failed!'
  endif
end subroutine inverse