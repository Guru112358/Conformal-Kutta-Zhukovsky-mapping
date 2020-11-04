program Kutta_Zhukovsky
implicit none
!variables involved
!z=x+i*y
!zeta=ksi+i*eta
!The K-Z transformation is given by:
!z=zeta+(1/zeta)=x+iy
complex::z,zeta
real,parameter::pi=4*ATAN(1.)
real,allocatable,dimension(:)::x,y,ksi,eta,theta    
real::theta0,thetamax,dtheta,r,ksi0,eta0
integer::niter,i
theta0=0.0
thetamax=2*pi
dtheta=0.01
niter=INT((thetamax-theta0)/dtheta)
write(*,*)
write(*,*)
write(*,*)"!------------------------------------------------------------------------------------------------------!"
write(*,*)"|  This program performs conformal mapping(Kutta_Zhukovsky) from  a circle to an airfoil contour...... |"
write(*,*)"|   or a  potato shape,depending on coordinates                                                        |"
write(*,*)"|   Try circle coordinates of (-0.2,0.2) with a radius of 1.2                                          |"
write(*,*)"!------------------------------------------------------------------------------------------------------!"
write(*,*)
write(*,*)
write(*,*)"Enter the x and y coordinates of the center of the circle to be transformed separated by comma:- "
write(*,*)
write(*,*)
read *,ksi0,eta0

write(*,*)"Enter the radius of the circle:-"
read *,r

open(1,file="circle.dat",status='replace')
open(2,file="kzfoil.dat",status='replace')
open(3,file="kzplot.plt",status='replace')
open(4,file="circle.plt",status='replace')

allocate(ksi(niter),eta(niter),theta(niter),x(niter),y(niter))

do i=1,niter
theta(i)=i*dtheta
ksi(i)=ksi0+r*cos(theta(i))
eta(i)=eta0+r*sin(theta(i))
zeta=complex(ksi(i),eta(i))
z=zeta+(1/zeta)
x(i)=real(z)
y(i)=imag(z)
!writing data to the files for plotting
write(2,*)x(i),y(i)
write(1,*)ksi(i),eta(i)
end do

!writing gnuplot script
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(4,*)"set xlabel 'ksi'"
write(4,*)"set ylabel 'eta'"
write(4,*)"set size ratio 1"
write(4,*)"set grid"
write(4,*)"set title 'Circle before transformation'"
write(4,*)"plot 'circle.dat' with line lt rgb 'red'"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(3,*)"set xlabel 'x'"
write(3,*)"set ylabel 'y'"
write(3,*)"set size ratio 0.1"
write(3,*)"set grid"
write(3,*)"set title 'contour after transformation'"
write(3,*)"plot 'kzfoil.dat' with line lt rgb 'blue'"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL SYSTEM('gnuplot -p circle.plt')
CALL SYSTEM('gnuplot -p kzplot.plt')

deallocate(ksi,eta,theta,x,y)


 close(1)
 close(2)
 close(3)
 close(4)

end program Kutta_Zhukovsky
