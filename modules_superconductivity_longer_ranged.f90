
!*****************************General parameters******************************************
module parameters_as
    implicit none

    !General fixed parameters 
    integer, parameter :: nsol=4,Nav=501
    !Constants and temperature
    double precision,parameter ::Pi=4.d0*atan(1.d0),temp=1.e-5 ,beta=1/temp

    ! tight binding parameters.

    !double precision, parameter :: t=1, t2=-0.1636/(0.5951), t3=0.0519/(0.5951)&         ! tb1 from M. Norman article
    !&, t4=0.1117/2/(0.5951), t5=-0.0510/(0.5951)   
    !double precision, parameter ::t=1 ,t2=-0.2368/(0.6798), t3=0.0794/(0.6798), &    ! tb2 from M. Norman article
    !&t4=-0.0343/2/(0.6798), t5=-0.0011/(0.6798)  
    double precision, parameter :: t=1, t2=-0.1461/(0.5881), t3=-0.0095/(0.5881), &      ! tb3 from M. Norman article
    &t4=0.1298/2/(0.5881), t5=-0.0069/(0.5881)   
    !double precision, parameter :: t=1, t2=-0.0740/(0.7823), t3=0.0587/(0.7823)&         ! tb4 from M. Norman article
    !&, t4=0.1398/2/(0.7823), t5=0.0174/(0.7823)  

    !double precision,parameter  :: t=1, t2=0, t3=0, t4=0, t5=0   !Nearest neighbors case
    !Coulomb interaction and density per site
    double precision,parameter ::U=abs(8*t),n=0.6000d0
    !initial guess
    double precision :: xx(12),kx,ky, x01,x02,x03,x04, x05, x06, x07,x08,x09,x010,x011,x012    
end module parameters_as

    !*********************************Matrices definition**************************************

module matrices_as
    use parameters_as
    implicit none

    !Fourier factors for a square lattice for nearest, next nearest, next-next nearest [...] neighbors. ga imposes d-wave symmetry.
    double precision :: al,al2,al3,al4,al5,ga
    !Elements of the M matrix. 
    double precision :: m11,m12,m22,m13,eps 
    !Elements of the I^{-1} matrix
    double precision ::  alpha1, beta1
    !Eigenvalues of the E matrix
    double precision, dimension (int(nsol)) ::e          
    
    contains
      subroutine init_matrices()

        !Fourier factors for a square lattice
        al=(cos(kx)+cos(ky))/2.d0  
        al2=cos(kx)*cos(ky)
        al3=(cos(2*kx)+cos(2*ky))/2
        al4=(cos(2*kx)*cos(ky)+cos(kx)*cos(2*ky))/2
        al5=(cos(2*kx)*cos(2*ky))
        !d-wave Fourier factor for superconductivity.
        ga=(cos(kx)-cos(ky))/2.d0

        !Elements of I^{-1} matrix
        alpha1=2.d0/(2.d0-n)
        beta1=2.d0/n

        !tight-binding energy
        eps=4.d0*t*al+4.d0*t2*al2+4.d0*t3*al3+8*t4*al4+4*t5*al5   

        !Elements of the M matrix
        m11= -xx(1)*(1.d0-n/2.d0)-eps*(1.d0-n)-4.d0*t*al*xx(2)-4*t2*al2*xx(3)-4*t3*al3*xx(4)-8*t4*al4*xx(5)-4*t5*al5*xx(6)&
        & - 4*t*xx(7)-4*t2*xx(8)-4*t3*xx(9)-8*t4*xx(10)-4*t5*xx(11) 

        m12=-eps*(n/2.d0)+4.d0*t*al*xx(2)+4*t2*al2*xx(3)+4*t3*al3*xx(4)+8*t4*al4*xx(5)+4*t5*al5*xx(6) +4.d0*t*xx(7) +4.d0*t2*xx(8)&
        & + 4*t3*xx(9)+8*t4*xx(10)+4*t5*xx(11)                           

        m22=-(xx(1)-U)*n/2.d0 -4.d0*t*al*xx(2)-4.d0*t2*al2*xx(3) -4*t3*al3*xx(4)-8*t4*al4*xx(5)-4*t5*al5*xx(6) &
        &-4*t*xx(7)-4*t2*xx(8)-4*t3*xx(9)-8*t4*xx(10)-4*t5*xx(11)

        m13=-4.d0*t*ga*xx(12)   

        !Eigenvalues of the E matrix: an analytical form have been directly obtained with mathematica (speed up the code).

        e(1)=(-(Sqrt(2.d0)*Sqrt(4*m13**2 + m22**2*(-2 + n)**2 + n*(-2*m12**2*(-2 + n) + m11**2*n) - Sqrt(16*m13**4 +& 
            &4*m13**2*(8*m22**2 + 8*(m12 - m22)*(m12 + m22)*n - 2*(m11**2 + 8*m11*m12 + 10*m12**2 + 4*m11*m22 + &
            &8*m12*m22 + m22**2)*n**2 + 4*(m11 + 2*m12 + m22)**2*n**3 - (m11 + 2*m12 + m22)**2*n**4) + (2*m22 +&
            & m11*n - m22*n)**2*(8*m12**2*n + (-2*m22 + (m11 - 2*m12 + m22)*n)*(-2*m22 + (m11 + 2*m12 + m22)*n)))))/((-2.d0 + n)*n))
            

        e(2)=(Sqrt(2.d0)*Sqrt(4*m13**2 + m22**2*(-2 + n)**2 + n*(-2*m12**2*(-2 + n) + m11**2*n) - Sqrt(16*m13**4 + &
            &4*m13**2*(8*m22**2 + 8*(m12 - m22)*(m12 + m22)*n - 2*(m11**2 + 8*m11*m12 + 10*m12**2 + 4*m11*m22 +&
            & 8*m12*m22 + m22**2)*n**2 + 4*(m11 + 2*m12 + m22)**2*n**3 - (m11 + 2*m12 + m22)**2*n**4) + (2*m22 +&
            & m11*n - m22*n)**2*(8*m12**2*n + (-2*m22 + (m11 - 2*m12 + m22)*n)*(-2*m22 + (m11 + 2*m12 + m22)*n)))))/((-2 + n)*n)


        e(3)=(-(Sqrt(2.d0)*Sqrt(4*m13**2 + m22**2*(-2 + n)**2 + n*(-2*m12**2*(-2 + n) + m11**2*n) + Sqrt(16*m13**4 + &
            &4*m13**2*(8*m22**2 + 8*(m12 - m22)*(m12 + m22)*n - 2*(m11**2 + 8*m11*m12 + 10*m12**2 + 4*m11*m22 +&
            & 8*m12*m22 + m22**2)*n**2 + 4*(m11 + 2*m12 + m22)**2*n**3 - (m11 + 2*m12 + m22)**2*n**4) + (2*m22 +&
            & m11*n - m22*n)**2*(8*m12**2*n + (-2*m22 + (m11 - 2*m12 + m22)*n)*(-2*m22 + (m11 + 2*m12 + m22)*n)))))/((-2 + n)*n))

        e(4)=(Sqrt(2.d0)*Sqrt(4*m13**2 + m22**2*(-2 + n)**2 + n*(-2*m12**2*(-2 + n) + m11**2*n) + Sqrt(16*m13**4 + &
            &4*m13**2*(8*m22**2 + 8*(m12 - m22)*(m12 + m22)*n - 2*(m11**2 + 8*m11*m12 + 10*m12**2 + 4*m11*m22 +&
            & 8*m12*m22 + m22**2)*n**2 + 4*(m11 + 2*m12 + m22)**2*n**3 - (m11 + 2*m12 + m22)**2*n**4) + (2*m22 + &
            &m11*n - m22*n)**2*(8*m12**2*n + (-2*m22 + (m11 - 2*m12 + m22)*n)*(-2*m22 + (m11 + 2*m12 + m22)*n)))))/((-2 + n)*n)


    end subroutine init_matrices

end module matrices_as

!***************************************************************

module Cij_as
    use parameters_as
    use matrices_as
    implicit none
    !Correlation functions (on site)
    double precision :: C011,C012,C013,C014,C022,C024, C033, C034, C044, C023
    !Correlation functions (NN)
    double precision  :: C11_1,C12_1,C13_1,C14_1,C22_1,C24_1, C33_1, C34_1, C44_1, C23_1
    !Correlation functions (NNN) 
    double precision ::  C11_2,C12_2,C13_2,C14_2,C22_2,C24_2 , C33_2, C34_2, C44_2, C23_2
    double precision :: C11_3, C12_3, C13_3, C14_3, C22_3, C24_3, C33_3, C34_3, C44_3, C23_3
    double precision :: C11_4, C12_4, C13_4, C14_4, C22_4, C24_4, C33_4, C34_4, C44_4, C23_4
    double precision :: C11_5, C12_5, C13_5, C14_5, C22_5, C24_5, C33_5, C34_5, C44_5, C23_5

    contains
    subroutine coefficients(matCk5, matCk4, matCk3, matCk2, matCk1, matCk0)  
        integer :: i,j,k
        !Elements involved in spectral weights definitions, spectral weights
        double precision, dimension (nsol,nsol,nsol) :: lambda, kappa
        !Correlation functions matrices
        double precision, dimension(nsol,nsol),intent(out) :: matCk5, matCk4, matCk3, matCk2,matCk1, matCk0
        
        !definitions and initialization of spectral weights
        call init_matrices()
        do i=1,4
            do j=1,4
                do k=1,4
                lambda(i,j,k)=0.d0
                kappa(i,j,k)=0.d0
                enddo
            enddo
        enddo


        !The expressions below are direct obtained from mathematica to greatly speed up the code.

        do i=1,4
        lambda(i,1,1)=-0.5*(8*(e(i) + m11 + 2*m12)*m13**2 - 8*(m12 - m13)*(m12 + m13)*m22 + 8*(e(i) + m11)*m22**2 + &
        &4*e(i)*(m12 - m22)*(m12 + m22)*n - 2*e(i)**2*(e(i) + m11)*n**2 + e(i)**3*n**3)/n**2

        lambda(i,1,2)=(4*m12**3 + 8*m12*m13**2 + 4*m13**2*(e(i) + m11 + m22) + m12*(-2*(e(i) + m11) + e(i)*n)*(2*m22 &
        &+ e(i)*n))/((-2 + n)*n)

        lambda(i,1,3)=m13*(e(i)**2 - (4*(m12 + m22)**2)/n**2)

        lambda(i,1,4)=-((m13*(-2*(e(i) + m11 + m12) + e(i)*n)*(-2*(m12 + m22) + e(i)*n))/((-2 + n)*n))

        lambda(i,2,3)=-((m13*(2*(m11 + m12) + e(i)*(-2 + n))*(2*(m12 + m22) + e(i)*n))/((-2 + n)*n))

        lambda(i,2,2)=e(i)**2*m22 - (4*(m11*(-m12**2 + m13**2) + m11**2*m22 + m13**2*(2*m12 + m22)))/(-2 + n)**2 &
        &+ (e(i)**3*n)/2. + (2*e(i)*(-2*m13**2 + m12**2*(-2 + n) - m11**2*n))/(-2 + n)**2

        lambda(i,2,4)=(m13*(-4*(m11 + m12)**2 + e(i)**2*(-2 + n)**2))/(-2 + n)**2


        lambda(i,3,3)=-(e(i)**2*m11) - (e(i)**3*(-2 + n))/2. + (4*((m11 + 2*m12)*m13**2 + (-m12**2 + m13**2)*m22 &
        &+ m11*m22**2))/n**2 - (2*e(i)*(2*m13**2 - m22**2*(-2 + n) + m12**2*n))/n**2

        lambda(i,3,4)=(-4*m12**3 + 4*m13**2*(e(i) - m11 - m22) - m12*(8*m13**2 + (2*m11 + e(i)*(-2 + n))*(-2*m22 &
        &+ e(i)*n)))/((-2 + n)*n)

        lambda(i,4,4)=-(e(i)**2*m22) + (4*(m11*(-m12**2 + m13**2) + m11**2*m22 + m13**2*(2*m12 + m22)))/(-2 + n)**2 &
        &+ (e(i)**3*n)/2. + (2*e(i)*(-2*m13**2 + m12**2*(-2 + n) - m11**2*n))/(-2 + n)**2

        enddo

        !Spectral weights
        kappa(1,:,:)=lambda(1,:,:)/((e(1)-e(2))*(e(1)-e(3))*(e(1)-e(4)))     
        kappa(2,:,:)=lambda(2,:,:)/((e(2)-e(1))*(e(2)-e(3))*(e(2)-e(4)))
        kappa(3,:,:)=lambda(3,:,:)/((e(3)-e(1))*(e(3)-e(2))*(e(3)-e(4)))   
        kappa(4,:,:)=lambda(4,:,:)/((e(4)-e(1))*(e(4)-e(2))*(e(4)-e(3)))


        !On-site correlation function
        matCk0=((1.d0 + tanh(beta*e(1)/2.d0))*kappa(1,:,:)+(1.d0 + tanh(beta*e(2)/2.d0))*kappa(2,:,:) &
        &+(1.d0+tanh(beta*e(3)/2.d0))*kappa(3,:,:)+(1.d0+tanh(beta*e(4)/2.d0))*kappa(4,:,:)  )/2.d0 

        !Longer ranged correlation functions
        matCk1=cos(kx)*matCk0  !symmetric in x/y (d-wave symmetry)
        matCk2=al2*matCk0  
        matCk3=al3*matCk0
        matCk4=al4*matCk0
        matCk5=al5*matCk0

    end subroutine coefficients 

    
    
    subroutine averages(nstep)
      ! Average over Brillouin zone
      implicit none
      integer ::nstep,i,j
      double precision :: dx,eps
      double precision, dimension(nsol,nsol)::matCk1, matCk2, matCk5, matCk4, matCk3  ,matCk0
      double precision, dimension(nsol,nsol)::  outCk0, outCk1, outCk2, outCk3, outCk4, outCk5

      eps=-Pi/2
      dx= 2.d0*Pi/(nstep-1)
      outCk0=0.d0
      outCk1=0.d0
      outCk2=0.d0
      outCk3=0.d0
      outCk4=0.d0
      outCk5=0.d0
 
      
      do i=0,nstep-1
         kx=eps+dx*float(i)

         do j=0,nstep-1
            ky=eps+dx*float(j)

            call init_matrices()
            !Computation of the correlation function for each momentum
            call coefficients(matCk5, matCk4, matCk3, matCk2, matCk1,matCk0)  

            outCk1=outCk1 + matCk1  
            outCk2=outCk2 + matCk2
            outCk3=outCk3 + matCk3
            outCk4=outCk4 + matCk4
            outCk5=outCk5 + matCk5
            outCk0=outCk0 + matCk0
         enddo
      enddo
      !Averaged correlation functions over Brillouin zone.
      outCk1=outCk1/(nstep**2)
      outCk2=outCk2/(nstep**2)
      outCk3=outCk3/(nstep**2)
      outCk4=outCk4/(nstep**2)
      outCk5=outCk5/(nstep**2)
      outCk0=outCk0/(nstep**2)


      !We only return the coefficients that interest us for the self consistent equations
      C11_1=outCk1(1,1)  
      C12_1=outCk1(1,2)
      C22_1=outCk1(2,2)
      C13_1=outCk1(1,3)
      C14_1=outCk1(1,4)
      C24_1=outCk1(2,4)
      C33_1=outCk1(3,3)
      C34_1=outCk1(3,4)
      C44_1=outCk1(4,4)
      C23_1=outCk1(2,3)


      C11_2=outCk2(1,1)
      C12_2=outCk2(1,2)
      C22_2=outCk2(2,2)
      C13_2=outCk2(1,3)
      C14_2=outCk2(1,4)
      C24_2=outCk2(2,4)
      C33_2=outCk2(3,3)
      C34_2=outCk2(3,4)
      C44_2=outCk2(4,4)
      C23_2=outCk2(2,3)


      C11_3=outCk3(1,1)
      C12_3=outCk3(1,2)
      C22_3=outCk3(2,2)
      C13_3=outCk3(1,3)
      C14_3=outCk3(1,4)
      C24_3=outCk3(2,4)
      C33_3=outCk3(3,3)
      C34_3=outCk3(3,4)
      C44_3=outCk3(4,4)
      C23_3=outCk3(2,3)


      C11_4=outCk4(1,1)
      C12_4=outCk4(1,2)
      C22_4=outCk4(2,2)
      C13_4=outCk4(1,3)
      C14_4=outCk4(1,4)
      C24_4=outCk4(2,4)
      C33_4=outCk4(3,3)
      C34_4=outCk4(3,4)
      C44_4=outCk4(4,4)
      C23_4=outCk4(2,3)


      C11_5=outCk5(1,1)
      C12_5=outCk5(1,2)
      C22_5=outCk5(2,2)
      C13_5=outCk5(1,3)
      C14_5=outCk5(1,4)
      C24_5=outCk5(2,4)
      C33_5=outCk5(3,3)
      C34_5=outCk5(3,4)
      C44_5=outCk5(4,4)
      C23_5=outCk5(2,3)

      C011=outCk0(1,1)
      C012=outCk0(1,2)
      C022=outCk0(2,2)
      C013=outCk0(1,3)
      C014=outCk0(1,4)
      C024=outCK0(2,4)
      C023=outCk0(2,3)
      C033=outCk0(3,3)
      C034=outCk0(3,4)
      C044=outCk0(4,4)
      C023=outCk0(2,3)

      end subroutine averages
       
    end module Cij_as
    