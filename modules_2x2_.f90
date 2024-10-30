!**********************************GENERAL PARAMETERS*******************************
module parameters2x2
!General fixed parameters are here
   implicit none
   !nsol is the dimension of matrices (here 2, for the uniform case). Nav is the number of momentum points for averages.
   integer, parameter :: nsol=2,Nav=500
   !Constants and temperature
   real(8), parameter ::pi=4.d0*datan(1.d0), temp=0.001000d0, beta=1.0d0/temp
   ! Define hoppings parameters here
   real(8), parameter ::t=1.d0 
   ! Interaction strength and electron density (between 0 and 2, n=1 corresponds to half-filling).
   real(8),parameter ::U=8.000000d0, n=0.8

   ! All important matrices   
   real(8) :: M_matk(nsol,nsol), MatCk0_num(nsol,nsol), MatCk1_num(nsol,nsol), C0ij_num(nsol,nsol), C1ij_num(nsol,nsol)
   ! vector of the solutions xx: xx(1)=mu xx(2)=p xx(3)=e
   real(8), allocatable :: xx(:) 
   !Momentum variables & initial guess
   real(8) :: kx,ky,x01,x02,x03

end module parameters2x2

!************************MATRICES DEFINITIONS*********************************

module matrices2x2
    use parameters2x2
    implicit none
    contains

      subroutine init_matrices_num(eigsk)
          !Variables for llapack diagonalizer
          integer :: i,k,j,lwork,lda,ldvl,ldvr,info
          real(8),dimension(nsol) :: wr,wi
          real(8), dimension(nsol,nsol) :: vl,vr
          real(8), allocatable :: work(:)
          real(8) :: vr_inv(nsol,nsol)
          !eigsk contains the eigenvalues of the E matrix
          real(8) :: eigsk(nsol)
          !The inverse of I matrix and the E matrices. alphak contains the symmetry of the considered lattice (here a square lattice)
          real(8) :: I_matk(nsol,nsol),E_matk(nsol,nsol),alphak

          lda=nsol
          ldvl=nsol
          ldvr=nsol
          lwork=10*nsol+10
          I_matk=0.0d0
          allocate(work(lwork))
            alphak=(2.0d0*t*dcos(kx))+(2.0d0*t*dcos(ky))  !Fourier factor for a square lattice

            !Expression of the M matrix, derived analytically.
            M_matk(1,1)=-xx(1)*(1.0d0-n/2.0d0)-alphak*(1.0d0-n+xx(2))-4.0d0*t*xx(3)
            M_matk(1,2)=4.0d0*t*xx(3)-alphak*(n/2.0d0-xx(2))
            M_matk(2,1)=4.0d0*t*xx(3)-alphak*(n/2.0d0-xx(2))
            M_matk(2,2)=-(xx(1)-U)*n/2.0d0-4.0d0*t*xx(3)-alphak*xx(2)

            !Expression of the inverse of the I matrix.
            I_matk(1,1)=1.0d0/(1.0d0-(n/2.0d0))
            I_matk(2,2)=2.0d0/n

            !E matrix, obtained from the composite operators approximation.
            E_matk=matmul(M_matk,I_matk)
            !Diagonalization
            call dgeev('V','V',nsol,E_matk,nsol,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
            do i=1,nsol
              eigsk(i)=wr(i)
            enddo
            !Inverting the eigenvectors for the computation of spectral weights.
            call inverse(vr,vr_inv,nsol)
          deallocate(work)
          !Expression of the I matrix.
          I_matk(1,1)=1.0d0/I_matk(1,1)
          I_matk(2,2)=1.0d0/I_matk(2,2)
          vr_inv=matmul(vr_inv,I_matk)

          !Expression of the on-site momentum dependent correlation function.
          MatCk0_num=0.0d0
          do i=1,nsol
            do j=1,nsol
              do k=1,nsol
                MatCk0_num(i,j)=MatCk0_num(i,j)+(vr(i,k)*vr_inv(k,j)*1.0d0/2.d0*((1.d0+dtanh(beta*eigsk(k)/2.d0))))
              enddo
            enddo
          enddo
          !expression of the nearest-neighbors momentum dependent correlation function.
          MatCk1_num=(alphak*MatCk0_num)/(4.0d0*t)  
      end subroutine init_matrices_num
end module matrices2x2


!***************************************FUNCTIONS TO CALL MATCK AND AVERAGE IT******************************************
module Cij2x2
    use parameters2x2
    use matrices2x2
    implicit none
    !used for self-consistent equations:
    real(8)  :: C011,C012,C022   !Correlation functions on the same site
    real(8) :: C11_1, C12_1,C22_1 !Correlation function at different site
    contains

      subroutine averages_num
          implicit none
          integer :: i,j
          real(8), allocatable :: eigsk_num(:)
          allocate(eigsk_num(nsol))

          !We perform the average over momentum space of the correlation functions.
          C0ij_num=0.0d0
          C1ij_num=0.0d0
          do i=1,Nav
              do j=1,Nav
                kx=-pi+((2.0d0*pi*(i-1))/Nav)
                ky=-pi+((2.0d0*pi*(j-1))/Nav)
                call init_matrices_num(eigsk_num)  
                C0ij_num=C0ij_num+MatCk0_num
                C1ij_num=C1ij_num+MatCk1_num
              enddo
          enddo
          C0ij_num=C0ij_num/dfloat(Nav*Nav)
          C1ij_num=C1ij_num/dfloat(Nav*Nav)

          C11_1=C1ij_num(1,1)  
          C12_1=C1ij_num(1,2)
          C22_1=C1ij_num(2,2)

          C011=C0ij_num(1,1)
          C012=C0ij_num(1,2)
          C022=C0ij_num(2,2)

          deallocate(eigsk_num)

      end subroutine averages_num

end module Cij2x2

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


