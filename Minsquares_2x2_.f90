Program MinSquares
    !*****************************************************************
    !***************************************************************
    !*************************************************************
    use parameters2x2
    use matrices2x2
    use Cij2x2
    implicit None

    integer :: nn2
    !For storing the datas in a nice way.
    character(len=40) ::tsvFormat 
    !The self-consistent equations
    real(8), allocatable :: func(:)
    ! For measuring execution time
    real :: start, finish  
    ! For checking if minsquares converged
    LOGICAL check          
    call cpu_time(start)
    !Number of unknown parameters
    nn2=3  

  
    allocate(func(nn2),xx(nn2))

!_____________These are the  initial guess for solving the equations-----------------!      

!x1=mu,x2=p,x3=e,x4=chi           
x01=1.05882445d0
x02=0.13946512d0
x03=-0.04403503d0

!_________________________________________________!
!Initialization of the unknown parameters.
xx(1)=x01
xx(2)=x02
xx(3)=x03

! minimization using the Broydn method.
call broydn(xx,nn2,check)   
call funcv(nn2,func)

!Convergence check
if (check) then
  write(*,*) '!!!!!!!!!!!!!!Convergence problems!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) 'Reinitilize the program with the outputs from this run as new input guess'
endif
call cpu_time(finish) 

!Display results
write(*,*) '----------------------------Set parameters-------------------------------------------'
write(*,*) 'n',n,'U',U,'T',temp
write(*,*) 't',t,'Nk', Nav
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) '                        mu          ','                p              ','             e               '
write(*,*) 'Initial guess : ', x01,x02,x03
write(*,*) 'Final output  : ',xx(1), ",",xx(2),",",xx(3)
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'functions :', func(1), ",", func(2), ",", func(3)
write(*,*) 'time to run the code  :', finish-start

!This part can be uncommented to create a text file storing datas, with a tab separator.

!open(unit=1,file="datas_2x2U8_funcn.txt",position="append", action="write") 
!tsvFormat = '(*(G0.17,:,"'//achar(9)//'"))'
!write(1,tsvFormat) n, U, temp,  xx(1),xx(2),xx(3), func(1), func(2), func(3), C011,C012,&
!&C022, C11_1, C12_1, C22_1
!close(1)
      
deallocate(func,xx)

end Program MinSquares

!--------------------------------------------------------------------------------------------------------------------!

SUBROUTINE funcv(nn3,func_vec)
!Here are defined the self-consistent equations that will give the new values of the parameters
!as a function of correlation functions.
  use matrices2x2
  use Cij2x2
  implicit none

  integer::nn3,i
  real(8) :: func_vec(nn3)
  real(8) :: Dd, rho1_1, rho3_1
  real(8) :: na_new,nb_new,e1_new,e2_new, rho_delta, rho_S, phi, p1

  !initializing the functions.
  do i=1,nn3
    func_vec(i)=0.0d0
  enddo

  !Return the momentum averaged correlation functions for the previous set of parameters.
  call averages_num
  !Initialize the variables involved in the self-consistent equation of the p parameter.
  Dd=n/2.d0-C022-C012
  phi=-2.d0/(2.d0-n)*(C011+C012)+2.d0/n*(C012+C022)
  rho1_1=2.d0/(2.d0-n)*(C11_1+C12_1)**2+2.d0/n*(C22_1+C12_1)**2
  rho3_1=4.d0/(n*(2.d0-n))*(C11_1+C12_1)*(C12_1+C22_1)
  p1=n**2/4.d0-rho1_1/(1-phi**2)-rho1_1/(1-phi)-rho3_1/(1+phi)


  func_vec(1)=1.d0-C011-C022-2*C012-n/2.d0   !self consistent equation for n parameter.
  func_vec(2)=p1-xx(2)  !p (Roth scheme).
  !func_vec(2)=C012     !Uncomment for using Pauli scheme instead.
  func_vec(3)=C11_1-C22_1-xx(3)    !self consistent equation for e parameter.
  
  !Value of the parameters for the current iteration
  write(*,*)'param', xx(1), xx(2), xx(3)
  !Value of the errors |xx_new - xx_old|
  write(*,*) 'func:', func_vec(1), func_vec(2), func_vec(3)
  write(*,*)

END


    !--------------------------------------------------------------------------------------------------------------------!
       !Broydn functions (from llapack)

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
    
    

