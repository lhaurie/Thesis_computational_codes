Program MinSquares
    use parameters_as
    use matrices_as
    use Cij_as
      implicit None
      integer :: nn2=12       !Number of unknown parameters
      real :: start, finish  ! For measuring time
      character(len=40) ::tsvFormat  !For saving the datas in a text file properly if necessary
      double precision :: func(12)    !functions to minimize
      LOGICAL check          ! For checking if minsquares converged

      call cpu_time(start)
    

    !_____________These are the  initial guess for the solving the equations-----------------!
!x1=mu,x2=p1, x3=p2, x4=p3, x5=p4, x6=p5, x7=e1, x8=e2, x9=e3, x10=e4, x11=e5, x12=theta 


!How to initialize : if doping in hole (n<1) : a good thumb rule is to take mu close to zero. If doping in electron, take mu close to U.
! all the p are usually close to zero. Theta is usually taken around 10**-2. e are usually 10**-3.

!Example of initial guesses for tb3 at n=0.6
x01=0.51846773584650611E-001

x02=0.36972099091210284E-002
x03=0.82740538087644325E-001
x04=0.85889393300319375E-001
x05=0.86498971844751762E-001
x06=0.85185997956890691E-001

x07=-0.11966307419769243
x08=-0.30454982364105605E-001
x09=-0.25072182474478216E-001
x010=-0.28576152953000315E-001
x011=0.30474071002574609E-001

x012=0.12461840280814500E-001


!_________________________________________________!

! Initialization
xx(1)=x01
xx(2)=x02
xx(3)=x03
xx(4)=x04
xx(5)=x05
xx(6)=x06
xx(7)=x07
xx(8)=x08
xx(9)=x09
xx(10)=x010
xx(11)=x011
xx(12)=x012

!Minimization using broydn method
call broydn(xx,nn2,check) 
call funcv(nn2,func)

!display an error message if we haven't converge at the end.
if (check) then   
    write(*,*) '!!!!!!!!!!!!!!Convergence problems!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) 'Reinitilize the program with the outputs from this run as new input guess'
  endif
call cpu_time(finish)


!output when the code is done (converged or not)
write(*,*) 'Nav',Nav
write(*,*) '----------------------------Set parameters-------------------------------------------'
write(*,*) 'n',n,'U',U,'T',temp
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Init guess            mu               ','       p1           ','            p2          ','             e1&
&     ', '             e2                     '
write(*,*)   x01,x02, x03, x04, x05, x06, x07,x08,x09,x010,x011,x012
write(*,*)
write(*,*)
write(*,*) 'Final output:            mu               ','       p1           ','            p2          ','             e1&
&     ', '             e2                     '
write(*,*)   xx(1),",",xx(2),",", xx(3),",", xx(4),",", xx(5), ",", xx(6),",", xx(7),",",xx(8),",",xx(9),",", xx(10)&
& , ",", xx(11), "," , xx(12)
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Function value after the minimization'
write(*,*) '        f1                        f2                       f3                       f4                f5  '
write(*,*) func(1), func(2), func(3), func(4), func(5), func(6), func(7), func(8), func(9), func(10), func(11), func(12)
write(*,*)
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Delta:', C13_1+2.d0*C14_1+C24_1
write(*,*) 'C011',C011,'C012',C012,'C022',C022
write(*,*) 'C013',C013,'C014',C014,'C024',C024
write(*,*) 'C023',C023,'C033',C033,'C044',C044
write(*,*) 'C034',C034
write(*,*) 'C11_1',C11_1,'C12_1',C12_1,'C22_1',C22_1
write(*,*) 'C13_1', C13_1, 'C14_1', C14_1, 'C24_1', C24_1
write(*,*) 'C23_1', C23_1, 'C33_1', C33_1, 'C44_1', C44_1
write(*,*) 'C34_1', C34_1
write(*,*) 'C11_2',C11_2,'C12_2',C12_2,'C22_2',C22_2
write(*,*) 'C13_2', C13_2, 'C14_2', C14_2, 'C24_2', C24_2
write(*,*) 'C23_2',C23_2, 'C33_2', C33_2, 'C44_2', C44_2
write(*,*) 'C34_2', C24_2
write(*,*) 
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Time to run the program'
write(*,*) finish-start,'seconds',(finish-start)/60,'minutes'

!We can save the datas in a text file by uncommenting the line below.

!open(unit=1,file="tb4Allt1_U8t_T=0_Roth_as.txt",position="append", action="write")     
!tsvFormat = '(*(G0.17,:,"'//achar(9)//'"))'
!write(1,tsvFormat) n,  xx(1),xx(2),xx(3),xx(4),xx(5),xx(6), xx(7), xx(8), xx(9), xx(10), xx(11), xx(12), func(1),&
!& func(2), func(3), func(4), func(5), func(6),func(7), func(8), func(9), func(10), func(11), func(12), &
!&C011,C012,C013, C014,C022,C023,C024,C033, C034,C044, C11_1, C12_1, C13_1, C14_1, C22_1, C23_1, C24_1, C33_1,C34_1,C44_1
!close(1)
end Program MinSquares
    
    !--------------------------------------------------------------------------------------------------------------------!
    
    
    SUBROUTINE funcv(nn3,func_vec)
      !Here are defined the functions that we want to minimize for the self consistency.
    use matrices_as
    use Cij_as
    implicit none


    integer::nn3,i
    !Parameters involved in Roth self-consistent equations for p parameter.
    real(8) ::  Dd, phi, xi_1, xia
    double precision :: rho1_1, rho1_2, rho1_3, rho1_4, rho1_5
    double precision :: rho2_1, rho2_2, rho2_3, rho2_4, rho2_5 
    double precision :: rho3_1, rho3_2, rho3_3, rho3_4, rho3_5  
    double precision :: p1, p2, p3, p4, p5 

    !Self-consistent equations |xx_new - xx_old|
    real(8) :: func_vec(nn3)
    do i=1,nn3
      func_vec(i)=0.0d0
    enddo

    !Returns longer ranged correlation functions averaged on momentum.
    call averages(Nav)
    
    !Parameters in the self-consistent equations for p

    Dd=n/2.d0-C022 - C012       
    phi=(n**2-4.d0*Dd)/(n*(2.d0-n))   
    !phi=2/n*(C012+C022)- 2.d0/(2.d0-n)*(C012+C011)  !equivalent expressions for phi
    !phi=2.d0/(2.d0-n)*(C033+C034)-2.d0/n*(C034+C044)  !This comes from the symmetries relations C033=-C011, C044=-C022 and C034=-C012).
    rho1_1=2.d0/(2.d0-n)*(C11_1+C12_1)**2 + 2.d0/n*(C12_1+C22_1)**2  
    rho1_2=2.d0/(2.d0-n)*(C11_2+C12_2)**2 + 2.d0/n*(C12_2+C22_2)**2
    rho1_3=2.d0/(2.d0-n)*(C11_3+C12_3)**2 + 2.d0/n*(C12_3+C22_3)**2
    rho1_4=2.d0/(2.d0-n)*(C11_4+C12_4)**2 + 2.d0/n*(C12_4+C22_4)**2
    rho1_5=2.d0/(2.d0-n)*(C11_5+C12_5)**2 + 2.d0/n*(C12_5+C22_5)**2

    rho2_1=2.d0/(2.d0-n)*(C13_1+C14_1)*(C13_1+C23_1)+ 2.d0/n*(C23_1+C24_1)*(C14_1+C24_1) 
    rho2_2=2.d0/(2.d0-n)*(C13_2+C14_2)*(C13_2+C23_2)+ 2.d0/n*(C23_2+C24_2)*(C14_2+C24_2)  
    rho2_3=2.d0/(2.d0-n)*(C13_3+C14_3)*(C13_3+C23_3)+ 2.d0/n*(C23_3+C24_3)*(C14_3+C24_3)  
    rho2_4=2.d0/(2.d0-n)*(C13_4+C14_4)*(C13_4+C23_4)+ 2.d0/n*(C23_4+C24_4)*(C14_4+C24_4)  
    rho2_5=2.d0/(2.d0-n)*(C13_5+C14_5)*(C13_5+C23_5)+ 2.d0/n*(C23_5+C24_5)*(C14_5+C24_5)  

    rho3_1=4.d0/(n*(2.d0-n))*(C11_1+C12_1)*(C12_1+C22_1)       
    rho3_2=4.d0/(n*(2.d0-n))*(C11_2+C12_2)*(C12_2+C22_2)      
    rho3_3=4.d0/(n*(2.d0-n))*(C11_3+C12_3)*(C12_3+C22_3)
    rho3_4=4.d0/(n*(2.d0-n))*(C11_4+C12_4)*(C12_4+C22_4)
    rho3_5=4.d0/(n*(2.d0-n))*(C11_5+C12_5)*(C12_5+C22_5)



    !Parameter in the self-consistent equation for theta, superconductivity
    xi_1= 2.d0/(2.d0-n)*(C11_1+C12_1)*(C13_1+C14_1) + 2.d0/n*(C22_1+C12_1)*(C24_1+C14_1)    
    
    !Expression of p using Roth minimization. p1 = NN, p2 = NNN ...
    p1=n**2/4.d0-(rho1_1+phi*rho2_1)/(1.0d0-phi**2)-(rho1_1+rho2_1)/(1.0d0-phi)-(rho3_1)/(1.0d0+phi)!- n**2/4
    p2=n**2/4.d0-(rho1_2+phi*rho2_2)/(1.0d0-phi**2)-(rho1_2+rho2_2)/(1.0d0-phi)-(rho3_2)/(1.0d0+phi)!- n**2/4
    p3=n**2/4.d0-(rho1_3+phi*rho2_3)/(1.0d0-phi**2)-(rho1_3+rho2_3)/(1.0d0-phi)-(rho3_3)/(1.0d0+phi)!- n**2/4
    p4=n**2/4.d0-(rho1_4+phi*rho2_4)/(1.0d0-phi**2)-(rho1_4+rho2_4)/(1.0d0-phi)-(rho3_4)/(1.0d0+phi)!- n**2/4
    p5=n**2/4.d0-(rho1_5+phi*rho2_5)/(1.0d0-phi**2)-(rho1_5+rho2_5)/(1.0d0-phi)-(rho3_5)/(1.0d0+phi)!- n**2/4


    !Self consistent equations

    func_vec(1)=1.d0-C011-C022-n/2.d0 -2.d0*C012     !self consistent equation for n

    func_vec(2)=p1 - xx(2)              !self consistent equations for p
    func_vec(3)=p2 - xx(3)          
    func_vec(4)=p3 - xx(4)
    func_vec(5)=p4 - xx(5)
    func_vec(6)=p5 - xx(6)

    func_vec(7)=C11_1-C22_1- xx(7)    !self consistent equation for e
    func_vec(8)=C11_2-C22_2 - xx(8)
    func_vec(9)=C11_3-C22_3 - xx(9)
    func_vec(10)=C11_4 - C22_4 - xx(10)
    func_vec(11)=C11_5 - C22_5 - xx(11)

    func_vec(12)=(xx(12) - 2.0d0*(xi_1)/(phi+1.d0))  !Self equation for theta (gap equation)

    !Print out at each iteration the current parameter and minimization
    write(*,*)'param: mu=',xx(1),  ', p1=',xx(2), ', p2=',xx(3),', p3=',xx(4),', p4=', xx(5), 'p5=', xx(6)
    write(*,*) 'e1=', xx(7), 'e2=', xx(8), 'e3=', xx(9), 'e4=', xx(10), 'e5=', xx(11), 'theta=', xx(12)
    write(*,*)
    write(*,*) 'func :' , func_vec(1), func_vec(2), func_vec(3), func_vec(4), func_vec(5), func_vec(6)
    write(*,*) func_vec(7), func_vec(8), func_vec(9), func_vec(10), func_vec(11), func_vec(12)
    write(*,*)
    write(*,*)


       END
    
    


    !--------------------------------------------------------------------------------------------------------------------!
      !Broydn functions from llapack



          SUBROUTINE broydn(x,n,check)
          implicit none
          INTEGER n,nn,NP,MAXITS
          real(8) :: x(n),fvec,EPS,TOLF,TOLMIN,TOLX,STPMX 
          LOGICAL check
          PARAMETER (NP=40,MAXITS=100000,EPS=1.d-6,TOLF=1.d-6,TOLMIN=1.d-6,TOLX=EPS,STPMX=100000.)
          COMMON /newtv/ fvec(NP),nn
          SAVE /newtv/
         
    !    USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
          INTEGER i,its,j,k
          real(8) :: den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP),g(NP),p(NP),&
          &qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP)
          real(8) :: fmin
          LOGICAL restrt,sing,skip
          nn=n
          f=fmin(x)
          test=0.
          do 11 i=1,n
            if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
    11    continue
          if(test.lt..01*TOLF)then
            check=.false.
            return
          endif
          sum=0.
          do 12 i=1,n
            sum=sum+x(i)**2
    12    continue
          stpmax=STPMX*max(dsqrt(sum),float(n))
          restrt=.true.
          do 42 its=1,MAXITS
            if(restrt)then
              call fdjac(n,x,fvec,NP,r)
              call qrdcmp(r,n,NP,c,d,sing)
              if(sing) print*, 'singular Jacobian in broydn'
              do 14 i=1,n
                do 13 j=1,n
                  qt(i,j)=0.
    13          continue
                qt(i,i)=1.
    14        continue
              do 18 k=1,n-1
                if(c(k).ne.0.)then
                  do 17 j=1,n
                    sum=0.
                    do 15 i=k,n
                      sum=sum+r(i,k)*qt(i,j)
    15              continue
                    sum=sum/c(k)
                    do 16 i=k,n
                      qt(i,j)=qt(i,j)-sum*r(i,k)
    16              continue
    17            continue
                endif
    18        continue
              do 21 i=1,n
                r(i,i)=d(i)
                do 19 j=1,i-1
                  r(i,j)=0.
    19          continue
    21        continue
            else
              do 22 i=1,n
                s(i)=x(i)-xold(i)
    22        continue
              do 24 i=1,n
                sum=0.
                do 23 j=i,n
                  sum=sum+r(i,j)*s(j)
    23          continue
                t(i)=sum
    24        continue
              skip=.true.
              do 26 i=1,n
                sum=0.
                do 25 j=1,n
                  sum=sum+qt(j,i)*t(j)
    25          continue
                w(i)=fvec(i)-fvcold(i)-sum
                if(dabs(w(i)).ge.EPS*(dabs(fvec(i))+dabs(fvcold(i))))then
                  skip=.false.
                else
                  w(i)=0.
                endif
    26        continue
              if(.not.skip)then
                do 28 i=1,n
                  sum=0.
                  do 27 j=1,n
                    sum=sum+qt(i,j)*w(j)
    27            continue
                  t(i)=sum
    28          continue
                den=0.
                do 29 i=1,n
                  den=den+s(i)**2
    29          continue
                do 31 i=1,n
                  s(i)=s(i)/den
    31          continue
                call qrupdt(r,qt,n,NP,t,s)
                do 32 i=1,n
                  if(r(i,i).eq.0.) print*, 'r singular in broydn'
                  d(i)=r(i,i)
    32          continue
              endif
            endif
            do 34 i=1,n
              sum=0.
              do 33 j=1,n
                sum=sum+qt(i,j)*fvec(j)
    33        continue
              p(i)=-sum
    34      continue
            do 36 i=n,1,-1
              sum=0.
              do 35 j=1,i
                sum=sum-r(j,i)*p(j)
    35        continue
              g(i)=sum
    36      continue
            do 37 i=1,n
              xold(i)=x(i)
              fvcold(i)=fvec(i)
    37      continue
            fold=f
            call rsolv(r,n,NP,d,p)
            call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
            test=0.
            do 38 i=1,n
              if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
    38      continue
            if(test.lt.TOLF)then
              check=.false.
              return
            endif
            if(check)then
              if(restrt)then
                return
              else
                test=0.
                den=max(f,.5*n)
                do 39 i=1,n
                  temp=dabs(g(i))*max(dabs(x(i)),1.)/den
                  if(temp.gt.test)test=temp
    39          continue
                if(test.lt.TOLMIN)then
                  return
                else
                  restrt=.true.
                endif
              endif
            else
              restrt=.false.
              test=0.
              do 41 i=1,n
                temp=(dabs(x(i)-xold(i)))/max(dabs(x(i)),1.)
                if(temp.gt.test)test=temp
    41        continue
              if(test.lt.TOLX)return
            endif
    42    continue
          print*, 'MAXITS exceeded in broydn'
          END
    !--------------------------------------------------------------------------------------------------------------------!
          SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
          implicit none
          INTEGER n
          LOGICAL check
          real(8) :: f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
          PARAMETER (ALF=1.e-4,TOLX=1.e-7)
          EXTERNAL func
    !    USES func
          INTEGER i
          real(8) :: a,am,am2,amin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam
          check=.false.
          sum=0.
          do 11 i=1,n
            sum=sum+p(i)*p(i)
    11    continue
          sum=dsqrt(sum)
          if(sum.gt.stpmax)then
            do 12 i=1,n
              p(i)=p(i)*stpmax/sum
    12      continue
          endif
          slope=0.
          do 13 i=1,n
            slope=slope+g(i)*p(i)
    13    continue
          if(slope.ge.0.) print*, 'roundoff problem in lnsrch'
          test=0.
          do 14 i=1,n
            temp=dabs(p(i))/max(dabs(xold(i)),1.)
            if(temp.gt.test)test=temp
    14    continue
          amin=TOLX/test
          am=1.
    1     continue
            do 15 i=1,n
              x(i)=xold(i)+am*p(i)
    15      continue
            f=func(x)
            if(am.lt.amin)then
              do 16 i=1,n
                x(i)=xold(i)
    16        continue
              check=.true.
              return
            else if(f.le.fold+ALF*am*slope)then
              return
            else
              if(am.eq.1.)then
                tmplam=-slope/(2.*(f-fold-slope))
              else
                rhs1=f-fold-am*slope
                rhs2=f2-fold-am2*slope
                a=(rhs1/am**2-rhs2/am2**2)/(am-am2)
                b=(-am2*rhs1/am**2+am*rhs2/am2**2)/(am-am2)
                if(a.eq.0.)then
                  tmplam=-slope/(2.*b)
                else
                  disc=b*b-3.*a*slope
                  if(disc.lt.0.)then
                    tmplam=.5*am
                  else if(b.le.0.)then
                    tmplam=(-b+dsqrt(disc))/(3.*a)
                  else
                    tmplam=-slope/(b+dsqrt(disc))
                  endif
                endif
                if(tmplam.gt..5*am)tmplam=.5*am
              endif
            endif
            am2=am
            f2=f
            am=max(tmplam,.1*am)
          goto 1
          END
    !---------------------------------------------------------------------------------------------------------------!
          SUBROUTINE fdjac(n,x,fvec,np,df)
          implicit none
          INTEGER n,np,NMAX
          real(8) :: df(np,np),fvec(n),x(n),EPS
          PARAMETER (NMAX=40,EPS=1.e-4)
    !    USES funcv
          INTEGER i,j
          real(8) :: h,temp,f(NMAX)
          do 12 j=1,n
            temp=x(j)
            h=EPS*dabs(temp)
            if(h.eq.0.)h=EPS
            x(j)=temp+h
            h=x(j)-temp
            call funcv(n,f)
            x(j)=temp
            do 11 i=1,n
              df(i,j)=(f(i)-fvec(i))/h
    11      continue
    12    continue
          return
          END
    !---------------------------------------------------------------------------------------------------------------!
          FUNCTION fmin(x)
          implicit none
          INTEGER n,NP
          real(8) :: fmin,x(*),fvec
          PARAMETER (NP=40)
          COMMON /newtv/ fvec(NP),n
          SAVE /newtv/
    !    USES funcv
          INTEGER i
          real(8) :: sum
          call funcv(n,fvec)
          sum=0.0d0
          do 11 i=1,n
            sum=sum+fvec(i)**2
    11    continue
          fmin=0.5d0*sum
          return
          END
    !---------------------------------------------------------------------------------------------------------------!
          SUBROUTINE qrdcmp(a,n,np,c,d,sing)
          implicit none
          INTEGER n,np
          real(8) :: a(np,np),c(n),d(n)
          LOGICAL sing
          INTEGER i,j,k
          real(8) :: scale,sigma,sum,tau
          sing=.false.
          do 17 k=1,n-1
            scale=0.0d0
            do 11 i=k,n
              scale=max(scale,dabs(a(i,k)))
    11      continue
            if(scale.eq.0.)then
              sing=.true.
              c(k)=0.0d0
              d(k)=0.0d0
            else
              do 12 i=k,n
                a(i,k)=a(i,k)/scale
    12        continue
              sum=0.
              do 13 i=k,n
                sum=sum+a(i,k)**2
    13        continue
              sigma=dsign(dsqrt(sum),a(k,k))
              a(k,k)=a(k,k)+sigma
              c(k)=sigma*a(k,k)
              d(k)=-scale*sigma
              do 16 j=k+1,n
                sum=0.0d0
                do 14 i=k,n
                  sum=sum+a(i,k)*a(i,j)
    14          continue
                tau=sum/c(k)
                do 15 i=k,n
                  a(i,j)=a(i,j)-tau*a(i,k)
    15          continue
    16        continue
            endif
    17    continue
          d(n)=a(n,n)
          if(d(n).eq.0.0d0)sing=.true.
          return
          END
    !---------------------------------------------------------------------------------------------------------------!
          SUBROUTINE qrupdt(r,qt,n,np,u,v)
          implicit none
          INTEGER n,np
          real(8) :: r(np,np),qt(np,np),u(np),v(np)
    !    USES rotate
          INTEGER i,j,k
          do 11 k=n,1,-1
            if(u(k).ne.0.)goto 1
    11    continue
          k=1
    1     do 12 i=k-1,1,-1
            call rotate(r,qt,n,np,i,u(i),-u(i+1))
            if(u(i).eq.0.)then
              u(i)=dabs(u(i+1))
            else if(dabs(u(i)).gt.dabs(u(i+1)))then
              u(i)=dabs(u(i))*dsqrt(1.0d0+(u(i+1)/u(i))**2)
            else
              u(i)=dabs(u(i+1))*dsqrt(1.0d0+(u(i)/u(i+1))**2)
            endif
    12    continue
          do 13 j=1,n
            r(1,j)=r(1,j)+u(1)*v(j)
    13    continue
          do 14 i=1,k-1
            call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
    14    continue
          return
          END
    !---------------------------------------------------------------------------------------------------------------!
          SUBROUTINE rsolv(a,n,np,d,b)
          implicit none
          INTEGER n,np
          real(8) :: a(np,np),b(n),d(n)
          INTEGER i,j
          real(8) :: sum
          b(n)=b(n)/d(n)
          do 12 i=n-1,1,-1
            sum=0.0d0
            do 11 j=i+1,n
              sum=sum+a(i,j)*b(j)
    11      continue
            b(i)=(b(i)-sum)/d(i)
    12    continue
          return
          END
    !------------------------------------------------------------------------------------------------------------------!
          SUBROUTINE rotate(r,qt,n,np,i,a,b)
          implicit none
          INTEGER n,np,i
          real(8) :: a,b,r(np,np),qt(np,np)
          INTEGER j
          real(8) :: c,fact,s,w,y
          if(a.eq.0.)then
            c=0.
            s=dsign(1.0d0,b)
          else if(dabs(a).gt.dabs(b))then
            fact=b/a
            c=dsign(1.0d0/dsqrt(1.0d0+fact**2),a)
            s=fact*c
          else
            fact=a/b
            s=dsign(1.0d0/dsqrt(1.0d0+fact**2),b)
            c=fact*s
          endif
          do 11 j=i,n
            y=r(i,j)
            w=r(i+1,j)
            r(i,j)=c*y-s*w
            r(i+1,j)=s*y+c*w
    11    continue
          do 12 j=1,n
            y=qt(i,j)
            w=qt(i+1,j)
            qt(i,j)=c*y-s*w
            qt(i+1,j)=s*y+c*w
    12    continue
          return
          END
    !------------------------------------------------------------------------------------------------------------------------!
    !--------------------------------------------------------------------------------------------------------------------!
    subroutine jacobi(mat,v,d,siz,nrot)
    implicit none
    real(8), allocatable :: b(:) 
    real(8) :: mat(siz,siz),v(siz,siz),d(siz),z(siz),sm,er,h,t,s,c,g,theta,tau
    integer :: i,j,ip,iq,siz,nrot
    allocate(b(siz))
    do i=1,siz
        do j=1,siz
            if(i==j) then
            v(i,j)=1.0d0;
            else
            v(i,j)=0.0d0;
            endif
        enddo
    d(i)=mat(i,i);
    b(i)=d(i)
    z(i)=0.0d0;
    enddo
    nrot=0;
    do i=1,50
        !print*, i
        sm=0.0d0
        do ip=1,siz-1
            do iq=ip+1,siz
            sm=sm+dabs(mat(ip,iq))
            enddo
        enddo
        if(sm==0) RETURN
        if(i<4) then
        er=((0.2d0*sm)/(siz**2))
        else
        er=0.0d0
        endif
        do ip=1,siz-1
            do iq=ip+1,siz
                g=100.0*dabs(mat(ip,iq))
                if((i>4).and.(dabs(d(ip))+g==abs(d(ip))).and.(dabs(d(iq))+g==dabs(d(iq)))) then
                    mat(ip,iq)=0.0d0
                else if(dabs(mat(ip,iq))>er) then
                    h=d(iq)-d(ip)
                    if(dabs(h)+g==dabs(h))then
                        t=mat(ip,iq)/h
                    else
                        theta=0.5*h/mat(ip,iq)
                        t=1.0d0/(dabs(theta)+dsqrt(1.0d0+theta**2))
                        if(theta.lt.0.) t=-t
                    endif
                    c=1.0d0/dsqrt(1.0d0+t**2)
                    s=t*c
                    tau=s/(1.0d0+c)
                    h=t*mat(ip,iq)
                    z(ip)=z(ip)-h
                    z(iq)=z(iq)+h
                    d(ip)=d(ip)-h
                    d(iq)=d(iq)+h
                    mat(ip,iq)=0.0d0 
                    do j=1,ip-1
                        g=mat(j,ip)
                        h=mat(j,iq)
                        mat(j,ip)=g-s*(h+g*tau)
                        mat(j,iq)=h+s*(g-h*tau)
                    enddo
                    do j=ip+1,iq-1
                        g=mat(ip,j)
                        h=mat(j,iq)
                        mat(ip,j)=g-s*(h+g*tau)
                        mat(j,iq)=h+s*(g-h*tau)
                    enddo
                    do j=iq+1,siz
                        g=mat(ip,j)
                        h=mat(iq,j)
                        mat(ip,j)=g-s*(h+g*tau)
                        mat(iq,j)=h+s*(g-h*tau)
                    enddo
                    do j=1,siz
                        g=v(j,ip)
                        h=v(j,iq)
                        v(j,ip)=g-s*(h+g*tau)
                        v(j,iq)=h+s*(g-h*tau)
                    enddo
                    nrot=nrot+1
                endif
            enddo
        enddo
        do ip=1,siz
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.0d0
        enddo
    enddo
    deallocate(b)
    return
    end subroutine jacobi
    
    
