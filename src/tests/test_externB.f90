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
    !  Unit tests of the polynomial solvers
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
    !  Unit tests of the quartic solver
    !+
    !--------------------------------------------
    subroutine test_Bexternal(ntests,npass)
     use extern_Bfield, only:Bexternal,get_externalB_force
     use units,         only:set_units
     use physcon,       only:au,solarm
     use vectorutils,   only:cross_product3D
     integer, intent(inout) :: ntests,npass
     real :: rhoi,xi,yi,zi,eps
     real :: fextx,fexty,fextz
     real :: B(3),curlB(3),JcrossB(3)
     real :: dBxdx,dBxdy,dBxdz
     real :: dBydx,dBydy,dBydz
     real :: dBzdx,dBzdy,dBzdz
     integer :: nfail(3) !ierr removed
     real, parameter :: tol = 1.e-6
    
     if (id==master) write(*,"(/,a)") '--> checking external B field'
    
     call set_units(dist=au,mass=solarm,G=1.d0)
     xi = 0.3
     yi = 0.5
     zi = 0.6
     eps = 1.e-10
     hi = 1.

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
     call checkval(fextx,JcrossB(1),tol,nfail(1),'fextx = (curl B x B)_x')
     call checkval(fexty,JcrossB(2),tol,nfail(2),'fexty = (curl B x B)_y')
     call checkval(fextz,JcrossB(3),tol,nfail(3),'fextz = (curl B x B)_z')
     call update_test_scores(ntests,nfail,npass)
    
    end subroutine test_Bexternal
    
end module testexternB
    