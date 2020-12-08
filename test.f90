program ar2hpestest
implicit none
real*8:: rarh1,rarh2,rar2,ener,theta,abener
real*8,dimension(:)::dvdr(3)

call triatomic_coeff

rarh1=2.7d0
rarh2=2.7d0
rar2=5.2d0

call ar2hpes(rarh1,rarh2,rar2,ener,dvdr)

print*,ener,dvdr
!results
!  -101.93280575241496        3.8229350203716912        3.8229350203633921       -28.610974220097866     
end program
