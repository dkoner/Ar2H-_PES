!#####################################################################
!#                Ar2H+ POTENTIAL ENERGY SURFACE                     # 
!#####################################################################
!=====================================================================
!                                                                    |
! To initialize, first call the subroutine:                          |
!-------------------------                                           |   
! call triatomic_coeff                                               |
!-------------------------                                           |
! Now, call the main surface routine:                                |
!----------------------------------------------------------          |
! call ar2hpes(rarh1,rarh2,rar2,ener,dvdr)                           |
!----------------------------------------------------------          |
! rarh1 = Ar1---H   distance in a.u.                                 |
! rarh2 = Ar2---H   distance in a.u.                                 |
! rar2  = Ar1---Ar2 distance in a.u.                                 |
! ener  = interaction energy in kcal/mol                             |
! (zero is set at ArH^+ + Ar assymtote)                              |
! dvdr  = an array with three derivatives of potential w.r.t. bond   |
! distances. The sequnce is dv/drarh1, dv/drarh2, dv/drar2           |
! 13940 CCSD(T)/aug-cc-pVQZ energies for the triatom                 |
! Overall RMSE 0.057 kcal/mol                                        | 
! By Dr. Debasish Koner                                              |
!=====================================================================

module triatomic
implicit none
integer, parameter :: nc = 255
real*8,dimension(:) :: a(nc)
integer,dimension(:,:) :: indx(3,nc)
end module

subroutine triatomic_coeff
use triatomic
integer::i, dummy

open(unit=1122,file="coeff.dat")

do i=1,nc
!========================TRIATOMIC COEFFICIENTS=======================
 read(1122,*)dummy, indx(:,i), a(i)
end do

end subroutine

subroutine ar2hpes(rarh1,rarh2,rar2,ener,dvdr)
implicit none
real*8, intent(in) :: rarh1,rarh2,rar2
real*8, intent(out) :: ener
real*8,dimension(:), intent(out) :: dvdr(3)
real*8 ::  earh1, earh2, ear2, tbie, der1, der2, der3
real*8,dimension(:) :: der(3)

earh1 = 0.0d0
earh2 = 0.0d0
ear2 = 0.0d0
tbie = 0.0d0
der=0.0d0
dvdr=0.0d0
ener=0.0d0

call ArHPLUS_PES(rarh1,earh1,der1)
call ArHPLUS_PES(rarh2,earh2,der2)
call Ar2_PES(rar2,ear2,der3)
call Ar2H_potential(rar2,rarh1,rarh2,tbie,der)
ener = earh1 + earh2 + ear2 + tbie 

if (rar2 < 4.0d0) ener=50.0d0 !to avoid unphysical large negetive numbers

dvdr(1) = der(2) + der1
dvdr(2) = der(3) + der2
dvdr(3) = der(1) + der3

end subroutine

subroutine Ar2H_potential(r1,r2,r3,ymod,der)
use triatomic
implicit none
real*8, intent(in) :: r1, r2, r3
real*8, intent(out) :: ymod
real*8,dimension(:), intent(out) :: der(3)
integer :: i, j, k, s
real*8 :: exarh1, exarh2, exar2, barh, bar2, exp1, exp2, exp3

exarh1=0.0d0
exarh2=0.0d0
exar2=0.0d0
ymod=0.0d0
der=0.0d0

bar2=0.74877703001999962d0
barh=0.65032919103537346d0

exp1=exp(-bar2*r1)
exp2=exp(-barh*r2)
exp3=exp(-barh*r3)

exar2=r1*exp1
exarh1=r2*exp2
exarh2=r3*exp3

do s=1,nc
  i=indx(1,s) 
  j=indx(2,s) 
  k=indx(3,s)
 
  ymod=ymod+a(s)*(exar2**i)*(exarh1**j)*(exarh2**k)
  der(1)=der(1)+a(s)*i*(exar2**(i-1))*(exp1-bar2*exar2)*(exarh1**j)*(exarh2**k)
  der(2)=der(2)+a(s)*j*(exarh1**(j-1))*(exp2-barh*exarh1)*(exar2**i)*(exarh2**k)
  der(3)=der(3)+a(s)*k*(exarh2**(k-1))*(exp3-barh*exarh2)*(exarh1**j)*(exar2**i)
end do

end subroutine

!#####################################################
!# 426 data, RMSE = 3.3042371021757581E-003 kcal/mol #
!#####################################################
subroutine ArHPLUS_PES(x,ener,der1)
implicit none
real*8, intent(in) :: x
real*8, intent(out) :: ener, der1
real*8::ylong,yshort, x1,sum,exp1, fst, dx1dx
integer::i,j
integer,parameter::nc=16
real*8,dimension(:)::a(nc)
real*8,parameter::req=2.424356d0,rlong=20.0d0

sum=0.0d0
ener=0.0d0
ylong=0.0d0
yshort=0.0d0
a=0.0d0

a(1)=  9710.9175238821917d0
a(2)=  1.9794384644666676d0
a(3)= -0.96352235587778734d0
a(4)= -84.343982258623029d0
a(5)= -4003.0865610344745d0
a(6)= -95307.921037497421d0
a(7)=  1255426.2537212772d0
a(8)= -6687558.3259664821d0
a(9)=  19648362.817514554d0
a(10)= -33491723.718043514d0
a(11)=  31186956.766249385d0
a(12)= -12331956.121231772d0
a(13)= 0.73387194777361686d0
a(14)=  6937.1932356924553d0
a(15)=  22550.101121644293d0
a(16)=  7287.0505608211088d0

x1=x+rlong*dexp(-(x-req))
exp1=exp(-a(13)*x)
sum=x*exp1
fst=a(1)*exp(-a(2)*x)/x
yshort=fst
do i=1,10
yshort=yshort+a(i+2)*(sum**i)
end do

ylong=-0.5d0*a(14)/x1**4-0.5d0*a(15)/x1**6-a(16)/x1**6
ener=ylong+yshort

dx1dx=1.0d0-rlong*dexp(-(x-req))
der1=-fst*(1.0/x+a(2))
do i=1,10
der1=der1+a(i+2)*i*(sum**(i-1))*(exp1-a(13)*sum)
end do

der1=der1+dx1dx*(2.0d0*a(14)/x1**5+3.0d0*a(15)/x1**7+6.0d0*a(16)/x1**7)

end subroutine

!#####################################################
!# 337 data, RMSE   4.9234655962841420E-004 kcal/mol #
!#####################################################
subroutine Ar2_PES(x,ener,der1)
implicit none
real*8, intent(in) :: x
real*8, intent(out) :: ener, der1
integer::i,j,n2
real*8::x1, f2n,yshort,ylong,sum, exp1, fst, dx1dx, dylong, df2n, con1
real*8,parameter::req=7.182001d0,rlong=15.0d0
real*8,dimension(:)::a(7),b(13)

ener=0.0d0
yshort=0.0d0
ylong=0.0d0
der1=0.0d0
sum=0.0d0
a=0.0d0
b=0.0d0

a(1)= 0.89846986149845076d0
a(2)=  63.917085642125301d0
a(3)=  1634.7713740771280d0
a(4)=  78382.360890962562d0
a(5)=  12776515.872891780d0
a(6)=  4205269768.5094700d0
a(7)=  1691855940029.9009d0
b(1)=  69519.441995375237d0
b(2)= 0.93469355317233127d0
b(3)= -2.4033288570738978d0
b(4)= -2128.1464814224214d0
b(5)=  3619.1359704237043d0
b(6)= -105865.55043249654d0
b(7)=  815622.07877514116d0
b(8)= -4062281.5581580843d0
b(9)=  12887602.124851679d0
b(10)= -25681977.394532215d0
b(11)=  29297639.883367710d0
b(12)= -14981327.190350780d0
b(13)= 0.63827892349879545d0

x1=x+rlong*dexp(-(x-req))
exp1=exp(-b(13)*x)
sum=x*exp1
fst=b(1)*exp(-b(2)*x)/x
yshort=fst
do i=1,10
yshort=yshort+b(i+2)*(sum**i)
end do

dx1dx=1.0d0-rlong*dexp(-(x-req))
der1=-fst*(1.0/x+b(2))
do i=1,10
der1=der1+b(i+2)*i*(sum**(i-1))*(exp1-b(13)*sum)
end do

ylong=0.0d0
dylong=0.0d0
do j=3,8
  n2=j*2
  con1=f2n(x1,n2,a(1))*a(j-1)/x1**(n2-1)
  ylong=ylong+con1/x1
  dylong=dylong+df2n(x1,n2,a(1))*a(j-1)/x1**n2-dfloat(n2)*con1/x1**2
end do
ener=yshort-ylong*627.509d0

der1=der1-dx1dx*dylong*627.509d0

end subroutine

function f2n(x,n2,smlb)
implicit none
real*8, intent(in) :: x, smlb
integer, intent(in) :: n2
real*8::f2n,sum1,fact
integer::i

sum1=0.0d0
do i=0,n2
sum1=sum1+(smlb*x)**i/fact(i)
end do
f2n=1.0d0-dexp(-smlb*x)*sum1
end function

function df2n(x,n2,smlb)
implicit none
real*8, intent(in) :: x, smlb
integer, intent(in) :: n2
real*8::df2n,sum1,fact,sum2
integer::i

sum1=0.0d0
do i=1,n2
sum1=sum1+dfloat(i)*smlb**i*x**(i-1)/fact(i)
end do

sum2=0.0d0
do i=0,n2
sum2=sum2+(smlb*x)**i/fact(i)
end do

df2n=-dexp(-smlb*x)*sum1+smlb*dexp(-smlb*x)*sum2
end function

function fact(x)
implicit none
integer, intent(in) :: x
real*8::fact
integer::i

fact=1.0d0
if (x>0)then
do i=1,x
fact=fact*dfloat(i)
end do
end if

end function
