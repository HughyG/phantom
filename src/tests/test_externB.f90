!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testexternB
!
! Unit tests of the polynomial solver modules
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, quartic, testutils
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_externB
 
 private

contains
    !--------------------------------------------
    !+
    !  Unit tests of the Bexternal Module
    !+
    !--------------------------------------------
subroutine test_externB(ntests,npass)
 integer, intent(inout) :: ntests,npass
 
 if (id==master) write(*,"(/,a)") '--> TESTING external magnetic field'
 
 call test_Bexternal(ntests,npass)
 
 if (id==master) write(*,"(/,a)") '<-- magnetic field TESTS COMPLETE'

end subroutine test_externB
    
    !--------------------------------------------
    !+
    !  Unit tests of the Bexternal Module
    !+
    !--------------------------------------------
subroutine test_Bexternal(ntests,npass)
 use extern_Bfield, only:Bexternal,get_externalB_force,gradBexternal
 use units,         only:set_units
 use physcon,       only:au,solarm,pi
 use vectorutils,   only:cross_product3D
 use options,           only:divergence_advection
 integer, intent(inout) :: ntests,npass
 real :: rhoi,xi,yi,zi,eps,r0
 real :: fextx,fexty,fextz
 real :: B(3),curlB(3),JcrossB(3),Bchk(3)
 real :: dBxdx,dBxdy,dBxdz,sqrtscaleFactor
 real :: dBydx,dBydy,dBydz,rad
 real :: dBzdx,dBzdy,dBzdz
 real :: dBextxdx,dBextxdy,dBextxdz
 real :: dBextydx,dBextydy,dBextydz
 real :: dBextzdx,dBextzdy,dBextzdz
 integer :: nfail(24),i !ierr removed
 real, parameter :: tol = 1.e-6
 
 if (id==master) write(*,"(/,a)") '--> checking external B field'

 call set_units(dist=au,mass=solarm,G=1.d0)
 
 divergence_advection = .false.
 xi = 0.3
 yi = 0.5
 zi = 0.6
 eps = 1.e-10
 r0 = 0.25
 sqrtscaleFactor = 1.0/sqrt(4*pi)

 do i = 0,6
  xi = 1.12 - 0.2*i
  yi = -0.942 + 0.2*i
  zi = 0.5*((-1)**i)
  if (i==6) then
    xi = 0.1
    yi = -0.1
    zi = 0.05
  endif
  B(1) = Bexternal(xi,yi,zi,1)
  B(2) = Bexternal(xi,yi,zi,2)
  B(3) = Bexternal(xi,yi,zi,3)
 
  dBxdx = (Bexternal(xi+eps,yi,zi,1) - B(1))/eps
  dBxdy = (Bexternal(xi,yi+eps,zi,1) - B(1))/eps
  dBxdz = (Bexternal(xi,yi,zi+eps,1) - B(1))/eps
  dBydx = (Bexternal(xi+eps,yi,zi,2) - B(2))/eps
  dBydy = (Bexternal(xi,yi+eps,zi,2) - B(2))/eps
  dBydz = (Bexternal(xi,yi,zi+eps,2) - B(2))/eps
  dBzdx = (Bexternal(xi+eps,yi,zi,3) - B(3))/eps
  dBzdy = (Bexternal(xi,yi+eps,zi,3) - B(3))/eps
  dBzdz = (Bexternal(xi,yi,zi+eps,3) - B(3))/eps
  rhoi = 1.
 
  curlB(1) = dBzdy - dBydz
  curlB(2) = dBxdz - dBzdx
  curlB(3) = dBydx - dBxdy
 
  call cross_product3D(curlB,B,JcrossB)
  call get_externalB_force(xi,yi,zi,rhoi,fextx,fexty,fextz)
 
  print *, JcrossB(1), fextx
  print *, JcrossB(2), fexty
  print *, JcrossB(3), fextz
 
  call gradBexternal(xi,yi,zi,dBextxdx,dBextxdy,dBextxdz, &
         dBextydx,dBextydy,dBextydz,dBextzdx,dBextzdy,dBextzdz)

  call checkval(dBxdx,dBextxdx,tol,nfail(1),'dBextxdx numerical vs analytical')
  call checkval(dBxdy,dBextxdy,tol,nfail(2),'dBextxdy numerical vs analytical')
  call checkval(dBxdz,dBextxdz,tol,nfail(3),'dBextxdz numerical vs analytical')
  call checkval(dBydx,dBextydx,tol,nfail(4),'dBextydx numerical vs analytical')
  call checkval(dBydy,dBextydy,tol,nfail(5),'dBextydy numerical vs analytical')
  call checkval(dBydz,dBextydz,tol,nfail(6),'dBextydz numerical vs analytical')
  call checkval(dBzdx,dBextzdx,tol,nfail(7),'dBextzdx numerical vs analytical')
  call checkval(dBzdy,dBextzdy,tol,nfail(8),'dBextzdy numerical vs analytical')
  call checkval(dBzdz,dBextzdz,tol,nfail(9),'dBextzdz numerical vs analytical')
  call update_test_scores(ntests,nfail(1:9),npass)
 
  call checkval(fextx,JcrossB(1),tol,nfail(10),'fextx = (curl B x B)_x')
  call checkval(fexty,JcrossB(2),tol,nfail(11),'fexty = (curl B x B)_y')
  call checkval(fextz,JcrossB(3),tol,nfail(12),'fextz = (curl B x B)_z')
  call update_test_scores(ntests,nfail(10:12),npass)

  divergence_advection = .true.
 
 
  B(1) = Bexternal(xi,yi,zi,1)
  B(2) = Bexternal(xi,yi,zi,2)
  B(3) = Bexternal(xi,yi,zi,3)

  rad = sqrt(xi**2 + yi**2 + zi**2)
  
  if (rad < r0) then
    Bchk(1) = sqrtscaleFactor*((rad/r0)**8 - 2*(rad/r0)**4 + 1)
    Bchk(2) = 0.
    Bchk(3) = sqrtscaleFactor
  else 
    Bchk(1) = 0.
    Bchk(2) = 0.
    Bchk(3) = 0.
  endif

  call checkval(B(1),Bchk(1),tol,nfail(13),'Bx')
  call checkval(B(2),Bchk(2),tol,nfail(14),'By')
  call checkval(B(3),Bchk(3),tol,nfail(15),'Bz')
  call update_test_scores(ntests,nfail(13:15),npass)

  dBxdx = (Bexternal(xi+eps,yi,zi,1) - B(1))/eps
  dBxdy = (Bexternal(xi,yi+eps,zi,1) - B(1))/eps
  dBxdz = (Bexternal(xi,yi,zi+eps,1) - B(1))/eps
  dBydx = (Bexternal(xi+eps,yi,zi,2) - B(2))/eps
  dBydy = (Bexternal(xi,yi+eps,zi,2) - B(2))/eps
  dBydz = (Bexternal(xi,yi,zi+eps,2) - B(2))/eps
  dBzdx = (Bexternal(xi+eps,yi,zi,3) - B(3))/eps
  dBzdy = (Bexternal(xi,yi+eps,zi,3) - B(3))/eps
  dBzdz = (Bexternal(xi,yi,zi+eps,3) - B(3))/eps

  call gradBexternal(xi,yi,zi,dBextxdx,dBextxdy,dBextxdz, &
             dBextydx,dBextydy,dBextydz,dBextzdx,dBextzdy,dBextzdz)

  call checkval(dBxdx,dBextxdx,tol,nfail(16),'dBextxdx numerical vs analytical')
  call checkval(dBxdy,dBextxdy,tol,nfail(17),'dBextxdy numerical vs analytical')
  call checkval(dBxdz,dBextxdz,tol,nfail(18),'dBextxdz numerical vs analytical')
  call checkval(dBydx,dBextydx,tol,nfail(19),'dBextydx numerical vs analytical')
  call checkval(dBydy,dBextydy,tol,nfail(20),'dBextydy numerical vs analytical')
  call checkval(dBydz,dBextydz,tol,nfail(21),'dBextydz numerical vs analytical')
  call checkval(dBzdx,dBextzdx,tol,nfail(22),'dBextzdx numerical vs analytical')
  call checkval(dBzdy,dBextzdy,tol,nfail(23),'dBextzdy numerical vs analytical')
  call checkval(dBzdz,dBextzdz,tol,nfail(24),'dBextzdz numerical vs analytical')
  call update_test_scores(ntests,nfail(16:24),npass)  
enddo

end subroutine test_Bexternal
    
end module testexternB
    