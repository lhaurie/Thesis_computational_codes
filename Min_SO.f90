Program MinSquares
    use parameters_so
    use matrices_so
    use Cij_so
      implicit None
      integer :: nn2=9     !Number of unknown parameters
      real :: start, finish  ! For measuring time
      character(len=40) ::tsvFormat  !For saving the datas in a text file properly if necessary
      double precision :: func(9)    !functions to minimize
      LOGICAL check          ! For checking if minsquares converged

      call cpu_time(start)
    

    !_____________These are the  initial guess for the solving the equations-----------------!
!x1=mu, x2=delta_n ,x3=p_x, x4=p_y, x5=e_x, x6=e_y, x7=exy, x8=eyx, x9=pxy

!Example of initial guess : OSMP at n=1.85, U=9.8, lambda=0.59
x01=0.87089480392482466      
x02=-0.15039502864855667       
x03=6.9469054762172314E-002  
x04=0.11711846590413325      
x05=-0.10162455748810664       
x06=-9.4801007426971885E-004  
x07=-3.1249767855757948E-002  
x08=-4.0304162437661224E-002  
x09=0.15147868252527452

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

!Minimization
call broydn(xx,nn2,check) 
call funcv(nn2,func)


if (check) then   !display an error message if we haven't converged properly.
    write(*,*) '!!!!!!!!!!!!!!Convergence problems!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) 'Reinitilize the program with the outputs from this run as new input guess'
endif
call cpu_time(finish)


!output when the code is done (converged or not)
write(*,*) 'Nav',Nav
write(*,*) '----------------------------Set parameters-------------------------------------------'
write(*,*) "tx",tx,'ty', ty, 'lambda', lambda, 'U',U,'T',1/beta
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Init guess         '
write(*,*)   x01,",",x02,",", x03,",", x04,",",x05,",",x06,",",x07,",",x08, ",", x09
write(*,*)
write(*,*)
write(*,*) 'Final output        '
write(*,*)   xx(1),",",xx(2),",", xx(3),",", xx(4),",",xx(5),",",xx(6),",",xx(7),",",xx(8),",",xx(9)
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Function value after the minimization'
write(*,*) func(1), func(2), func(3), func(4), func(5), func(6), func(7), func(8), func(9)
write(*,*)
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'C0 after the minimization'
write(*,*) C011, C012, C013, C014
write(*,*) C021, C022, C023, C024
write(*,*) C031, C032, C033, C034
write(*,*) C041, C042, C043, C044
write(*,*) "n",n ,'nx, ny', (n+xx(2))/2.d0, (n-xx(2))/2.d0, "Es",-8.d0*(C11+C12+C21+C22)-8.d0*(C33+C34+C43+C44)&
    &-4.d0*lambda*(C031+C032+C041+C042)+U*(C012+C022)+xx(1)*n
write(*,*) '-------------------------------------------------------------------------------------'
write(*,*) 'Time to run the program'
write(*,*) finish-start,'seconds',(finish-start)/60,'minutes'

!We can save the datas in a text file by uncommenting the line below.

!if ((abs(func(1))+abs(func(2))+abs(func(3))+abs(func(4))+abs(func(5))+abs(func(6))+abs(func(7))+abs(func(8)))/8.d0<=5.d-5) then
!open(unit=1,file="datas_2x2U0_funcn.txt",position="append", action="write") 
!tsvFormat = '(*(G0.17,:,"'//achar(9)//'"))'
!write(1,tsvFormat) n, lambda, tx, ty, temp,  xx(1),xx(2),xx(3),xx(4), xx(5), xx(6), xx(7), xx(8),& change it accordingly to Anurag datas
!&xx(9), func(1), func(2), func(3), func(4), func(5), func(6), func(7), func(8), func(9), C011,C012,C013, C014,&
!&C021,C022,C023,C024,C031, C032,C033, C034,C041, C042, C043, C044, C11, C12, C13, C14, C21, C22, C23,&
!& C24, C31,C32, C33,C34,C41,C42,C43, C44
!close(1)
!endif
end Program MinSquares
    
    !--------------------------------------------------------------------------------------------------------------------!
    
    
SUBROUTINE funcv(nn3,func_vec)
use matrices_so
use Cij_so
implicit none
integer:: nn3
real(8) :: px, py, pxy, DDxy, SSxy, nnxy
real(8) :: rhopx, rhopy, rhox, rhoy, phix, phiy, Omega0xx, Omega0yy, Omegax, Omegay, Omegapx, Omegapy
real(8) :: func_vec(nn3)
double precision :: alphax, betax, alphay, betay

do i=1,nn3
  func_vec(i)=0.0d0
enddo
!correlation functions averaged over momentum space
call averages(Nav)

!parameters nx and ny
nx=(n+xx(2))/2.d0
ny=(n-xx(2))/2.d0

!elements of px and py
alphax=2.d0/(2.d0-nx)           
betax=2.d0/nx
alphay=2.d0/(2.d0-ny)
betay=2.d0/ny

rhopx=alphax*(C11+C12)**2 + betax*(C12+C22)**2
rhopy=alphay*(C33+C34)**2 + betay*(C34+C44)**2

rhox=(alphax+betax)*(C11+C12)*(C22+C12)
rhoy=(alphay+betay)*(C33+C34)*(C44+C34)

phix=-alphax*(C011+C021)+betax*(C022+C012)
phiy=-alphay*(C033+C043)+betay*(C044+C034)

px=(nx)**2/4.d0 - rhopx/(1-phix**2) - rhopx/(1-phix) - rhox/(1+phix)
py=(ny)**2/4.d0 - rhopy/(1-phiy**2) - rhopy/(1-phiy) - rhoy/(1+phiy)


!Elements of pxy
Omega0xx=ny/2*(alphax*(C011+C021)-1)*(-1+phix)
Omega0yy=nx/2*(alphay*(C033+C043)-1)*(-1+phiy)
Omegax=  alphax*(C031+C041)**2 + betax*(C032+C042)**2 
Omegay=  alphay*(C013+C023)**2 + betay*(C014+C024)**2 
Omegapy= (alphay+betay)*(C013+C023)*(C024+C014)
Omegapx= (alphax+betax)*(C031+C041)*(C042+C032)

DDxy=  1/2.d0*(Omegapx/(1+phiy)+ Omegapy/(1+phix))
SSxy= - 1/2.d0*(Omegax/(1-phiy)+Omegay/(1-phix))
nnxy=1/2.d0*( (Omega0xx-Omegay)/(1-phix**2) + (Omega0yy-Omegax)/(1-phiy**2) )


pxy=nnxy + SSxy - DDxy


!x1=mu_x, x2=mu_y ,x3=p_x, x4=p_y, x5=e_x, x6=e_y, x7=u, x8=v
func_vec(1)=2.d0*(1.d0-C011-C022-C012-C021   +   1.d0-C033-C044-C043-C034) - n   !n=nx+ny = ntot, between 0 and 4 (x/y degen)
func_vec(2)=2.d0*(1.d0-C011-C022-C012-C021   -  (1.d0-C033-C044-C043-C034)) - xx(2)  !delta n
func_vec(3)=px - xx(3)  
func_vec(4)=py - xx(4)
func_vec(5)=C11-C22 - xx(5) !ex
func_vec(6)=C33-C44 - xx(6) !ey
func_vec(7)=C031 + C041 - C023 - C024 - xx(7)  !exy
func_vec(8)=C013 + C023 - C041 - C042 - xx(8)    !eyx
func_vec(9)=pxy - xx(9)    !pxy



!Print out at each iteration the current parameter and errors
write(*,*)'param:'
write(*,*)'mu=',xx(1),', delta_n=',xx(2), ', px=',xx(3),', py=',xx(4)
write(*,*)', ex=',xx(5), ', ey=',xx(6), ', exy=', xx(7), ', eyx=', xx(8)
write(*,*) ', pxy=',xx(9)
write(*,*) 'func :' 
write(*,*) func_vec(1), func_vec(2), func_vec(3), func_vec(4), func_vec(5)
write(*,*) func_vec(6), func_vec(7), func_vec(8), func_vec(9)
write(*,*)
write(*,*) '---------------------------------------------------------------'

END




!--------------------------------------------------------------------------------------------------------------------!
  !PAST THIS POINT, THERE ARE ONLY MINIMIZATION FUNCTION WE USE IN BROYDN THAT ANURAG COPIED, NOTHING UNDER THIS LINE IS INTERESTING
  !AND SHOULD BE MODIFIED !



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


