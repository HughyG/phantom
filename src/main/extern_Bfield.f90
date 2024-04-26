!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_Bfield
!
! This routine implements external forces related to various external
!   magnetic field configurations
!
! :References: None
!src/main/extern_Bfield.f90
! :Owner: Daniel Price
!
! :Runtime parameters:

   !   - Bstar_cgs      : *Magnetic field strength of the star in the disc at 1 AU in Gauss *

!
! :Dependencies: infile_utils, io, physcon, units
!
 implicit none
 private

 public :: get_externalB_force, externBfield, Bexternal
 public :: write_options_externB, read_options_externB, check_externB_settings

 ! default values for runtime options
 
 real,    private :: Bstar_cgs     = 0.0001

contains

!--------------------------------------------------------
!+
!    This subroutine and associated routines handles
!    everything to do with external B fields
!    (with non-zero spatial derivatives)
!
!    Originally written by D. Price & C. Toniolo, 2006
!    Rewritten for Phantom by D. Price 2015
!    Rewritten by H. Griffiths 2024
!+
!--------------------------------------------------------

subroutine externBfield(xi,yi,zi,vxi,vyi,vzi,rhoi, &
                        Bintx, Binty, Bintz, &
                        currJintx, currJinty, currJintz, &
                        Bextx, Bexty, Bextz, &
                        fextx, fexty, fextz, &
                        vdotgradBx, vdotgradBy, vdotgradBz, &
                        string)

 use io,  only:warning,error
 use physcon, only:pi
 use units, only:unit_Bfield
 real, intent(in) :: xi,yi,zi,vxi,vyi,vzi,rhoi
 real, intent(in) :: Bintx,Binty,Bintz
 real, intent(in) :: currJintx,currJinty,currJintz
 character(len=*), intent(in) :: string
 real, intent(out) :: Bextx,Bexty,Bextz,fextx,fexty,fextz
 real, intent(out) :: vdotgradBx,vdotgradBy,vdotgradBz

!   variables
! real :: currJextx,currJexty,currJextz
!  real :: currJintBextx,currJintBexty,currJintBextz
!  real :: currJextBintx,currJextBinty,currJextBintz
 real :: magMomx,magMomy,magMomz
 real :: rad,drdx,drdy,drdz,scaleFactor
 real :: dBextxdx,mdotr,drhoi
 real :: dBextxdy,dBextxdz,dBextydx,dBextydy
 real :: dBextydz,dBextzdx,dBextzdy,dBextzdz

!
!--get 1/rho
!
 if (rhoi > epsilon(rhoi)) then
    drhoi = 1./rhoi
 else
    drhoi = 0.
    call warning('externB',' rho <= tiny in externalBfield !!!')
    return
 endif
 !
 ! To stop uninitialized warnings
 !
 dBextxdx = 0.
 dBextxdy = 0.
 dBextxdz = 0.
 dBextydx = 0.
 dBextydy = 0.
 dBextydz = 0.
 dBextzdx = 0.
 dBextzdy = 0.
 dBextzdz = 0.
 drdx = 0.
 drdy = 0.
 drdz = 0.

!
!     find the radius of the current point and set the scale factor
!  
 rad = sqrt(xi**2+yi**2+zi**2)
 scaleFactor = 1./(4.*pi)
!
!     Set the magnetic moment of the dipole field
!  
 magMomx = 0.
 magMomy = 0.
 magMomz = 4.*pi*Bstar_cgs/(unit_Bfield)
 mdotr = xi*magMomx+yi*magMomy+zi*magMomz
 if (rad < tiny(0.)) then
   Bextx = 0.
   Bexty = 0.
   Bextz = 0.
 else
   Bextx = scaleFactor*(((3.*mdotr*xi)/(rad**5))-(magMomx/(rad**3)))
   Bexty = scaleFactor*(((3.*mdotr*yi)/(rad**5))-(magMomy/(rad**3)))
   Bextz = scaleFactor*(((3.*mdotr*zi)/(rad**5))-(magMomz/(rad**3)))
 endif
!============================================================!
 
   select case(trim(string))
   case('Bfield')
!
!--in this case return only the external B field
!  is returned in cartesian co-ordinates
!
!    
!     this is a dipole field centred about the origin in cartesian coordinates designed for testing
!

!
!
!
   
   case('fext')
!
!--get J_ext x B_ext as the external force
!
      call get_fext(magMomx,magMomy,magMomz,xi,yi,zi, &
            currJintx,currJinty,currJintz,Bextx,Bexty, &
            Bextx,Bintx,Binty,Bintz,scaleFactor,fextx,fexty,fextz,drhoi,rad,mdotr)

   case('all')
!
!  in this case return the external force
!  including mixed Jint x Bext terms and also
!  the term -v.grad Bext needed in the B evolution equation
!
      call get_fext(magMomx,magMomy,magMomz,xi,yi,zi, &
            currJintx,currJinty,currJintz,Bextx,Bexty, &
            Bextx,Bintx,Binty,Bintz,scaleFactor,fextx,fexty,fextz,drhoi,rad,mdotr)

      vdotgradBx = vxi*dBextxdx + vyi*dBextxdy + vzi*dBextxdz
      vdotgradBy = vxi*dBextydx + vyi*dBextydy + vzi*dBextydz
      vdotgradBz = vxi*dBextzdx + vyi*dBextzdy + vzi*dBextzdz

   case default
!
!        Default: no external B field; do nothing and return
!
          Bextx = 0.
          Bexty = 0.
          Bextz = 0.
          fextx = 0.
          fexty = 0.
          fextz = 0.
          vdotgradBx = 0.
          vdotgradBy = 0.
          vdotgradBz = 0.
      call error('externB','unknown string in call to externBfield')
   end select
 
end subroutine externBfield
!--------------------------------------------------------
!+
!  This subroutine computes the fext term
!  for the cartesian case
!+
!--------------------------------------------------------
subroutine get_fext(magMomx,magMomy,magMomz,xi,yi,zi,currJintx, &
            currJinty,currJintz,Bextx,Bexty,Bextz,Bintx,Binty, &
            Bintz,scaleFactor,fextx,fexty,fextz,drhoi,rad,mdotr)
 real, intent(in) :: magMomx,magMomy,magMomz
 real, intent(in) :: xi,yi,zi,scaleFactor,rad
 real, intent(in) :: currJintx,currJinty,currJintz
 real, intent(in) :: Bextx,Bexty,Bextz,drhoi
 real, intent(in) :: Bintx,Binty,Bintz,mdotr
 real, intent(out) :: fextx,fexty,fextz


 real :: dBextxdx,dBextxdy,dBextxdz
 real :: dBextydx,dBextydy,dBextydz
 real :: dBextzdx,dBextzdy,dBextzdz
 real :: CurrJextx,CurrJexty,CurrJextz
! real :: currJintBextx
! real :: currJintBexty,currJextBintx
! real :: currJintBextz,currJintBintx
 real :: currJextBextz,currJextBexty,currJextBextx
! real :: currJintBinty,currJextBinty
! real :: currJintBintz,currJextBintz

!
!  calculating partial derivatives of r
!
!  drdx = xi/rad
!  drdy = yi/rad
!  drdz = zi/rad
!
!  calculating quantities for the derivative
!
!  dmdotrdr = magMomx+magMomy+magMomz
!  drdotrdr = 2*(xi+yi+zi)
!
!  calculate partial derivatives of B with respect to r
! !
!  dBxdr = 3*scaleFactor*(((rad**2)*(dmdotrdr*xi+mdotr)-5*mdotr*0.5*drdotrdr*xi)/(rad**7) + (magMomx*0.5*drdotrdr)/(rad**5))
!  dBydr = 3*scaleFactor*(((rad**2)*(dmdotrdr*yi+mdotr)-5*mdotr*0.5*drdotrdr*yi)/(rad**7) + (magMomy*0.5*drdotrdr)/(rad**5))
!  dBydr = 3*scaleFactor*(((rad**2)*(dmdotrdr*zi+mdotr)-5*mdotr*0.5*drdotrdr*zi)/(rad**7) + (magMomz*0.5*drdotrdr)/(rad**5))
!  dBxdr = 3*scaleFactor*(magMomx/(rad**4)-5*magMomx*(xi**2)/(rad**6))
!  dBydr = 3*scaleFactor*(magMomy/(rad**4)-5*magMomy*(yi**2)/(rad**6))
!  dBzdr = 3*scaleFactor*(magMomz/(rad**4)-5*magMomz*(zi**2)/(rad**6))
!
!  calculate partial derivatives of B with respect to x,y and z
!
 dBextxdx = scaleFactor*(-6*magMomx*(xi**3)+9*((yi**2)+(zi**2))+3*(magMomy*yi+magMomz*zi)  &
   *(-4*(xi**2)+(yi**2)+(zi**2)))/(rad**7)
 dBextxdy = 3*scaleFactor*(-5*magMomz*(xi*yi*zi)+magMomy*xi*((xi**2)-4*(yi**2)+(zi**2))   &
   +magMomx*yi*(-4*(xi**2)+(yi**2)+(zi**2)))/(rad**7)
 dBextxdz = 3*scaleFactor*(-5*magMomy*(xi*yi*zi)+magMomz*xi*((xi**2)+(yi**2)-4*(zi**2))   &
   +magMomx*zi*(-4*(xi**2)+(yi**2)+(zi**2)))/(rad**7)
 dBextydx = 3*scaleFactor*(-5*magMomz*(xi*yi*zi)+magMomy*xi*((xi**2)-4*(yi**2)+(zi**2))   &
   +magMomx*yi*(-4*(xi**2)+(yi**2)+(zi**2)))/(rad**7)
 dBextydy = 3*scaleFactor*((magMomx*xi+magMomz*zi)*((xi**2)-4*(yi**2)+(zi**2))            &
   +magMomy*yi*(3*(xi**2)-2*(yi**2)+3*(zi**2)))/(rad**7)
 dBextydz = 3*scaleFactor*(-5*magMomx*(xi*yi*zi)+magMomz*yi*((xi**2)+(yi**2)-4*(zi**2))   &
   +magMomy*zi*((xi**2)-4*(yi**2)+(zi**2)))/(rad**7)
 dBextzdx = 3*scaleFactor*(-5*magMomy*(xi*yi*zi)+magMomz*xi*((xi**2)+(yi**2)-4*(zi**2))   &
   +magMomx*zi*(-4*(xi**2)+(yi**2)+(zi**2)))/(rad**7)
 dBextzdy = 3*scaleFactor*(-5*magMomx*(xi*yi*zi)+magMomz*yi*((xi**2)+(yi**2)-4*(zi**2))   &
   +magMomy*zi*((xi**2)-4*(yi**2)+(zi**2)))/(rad**7)
 dBextzdz = 3*scaleFactor*((magMomx*xi+magMomy*yi)*((xi**2)+(yi**2)-4*(zi**2))              &
   +magMomz*zi*(3*(xi**2)+3*(yi**2)-2*(zi**2)))/(rad**7)
!
!  calculate J_ext
!
 CurrJextx = dBextzdy-dBextydz
 CurrJexty = dBextxdz-dBextzdx
 CurrJextz = dBextydx-dBextxdy
! print *, dBextydx, dBextxdy,CurrJextz
!
!  calculate J_ext x B_ext
!
 currJextBextx = CurrJexty*Bextz - CurrJextz*Bexty
 currJextBexty = CurrJextz*Bextx - CurrJextx*Bextz
 currJextBextz = CurrJextx*Bexty - CurrJexty*Bextx 
!
!  calculate J_int x B_ext
!
!  currJintBextx = currJinty*Bextz - CurrJintz*Bextyn
!  currJintBexty = CurrJintz*Bextx - CurrJintx*Bextz
!  currJintBextz = CurrJintx*Bexty - CurrJinty*Bextx 
!
!  calculate J_ext x B_int
!
!  currJextBintx = currJexty*Bintz - CurrJextz*Binty
!  currJextBinty = CurrJextz*Bintx - CurrJextx*Bintz
!  currJextBintz = CurrJextx*Binty - CurrJexty*Bintx 
!
!  calculate J_int x B_int
!
!  currJintBintx = currJinty*Bintz - CurrJintz*Binty
!  currJintBinty = CurrJintz*Bintx - CurrJintx*Bintz
!  currJintBintz = CurrJintx*Binty - CurrJinty*Bintx 
!
!  calculate F = (Curl of B x B)/rho
! !
 fextx = dBextxdx
 fexty = dBextydy
 fextz = dBextzdz
!  fextx = 0*(currJextBextx)*(drhoi)
!  fexty = 0*(currJextBexty)*(drhoi)
!  fextz = 0*(currJextBextz)*(drhoi)
 return
end subroutine get_fext

!------------------------------------------------------------
!+
!  This function acts as an interface to externBfield
!  cutting out the dummy arguments needed on the first call
!+
!------------------------------------------------------------
real function Bexternal(xcoord,ycoord,zcoord,icomponent)
 use io, only:fatal
 integer, intent(in) :: icomponent
 real,    intent(in) :: xcoord,ycoord,zcoord
 real :: dumx,dumy,dumz,Bextx,Bexty,Bextz
 real :: dumrhoi,dumgx,dumgy,dumgz

 dumrhoi = 1.
 call externBfield(xcoord,ycoord,zcoord,0.,0.,0.,dumrhoi, &
                   0.,0.,0.,0.,0.,0., &
                   Bextx,Bexty,Bextz, &
                   dumx,dumy,dumz,dumgx,dumgy,dumgz,'Bfield')

 select case(icomponent)
 case(1)
    Bexternal= Bextx
 case(2)
    Bexternal= Bexty
 case(3)
    Bexternal= Bextz
 case default
    Bexternal = 0.
    call fatal('Bexternal','error in Bexternal call')
 end select

end function Bexternal

!------------------------------------------------------------
!+
!  This subroutine acts as an interface to externBfield
!  which returns only the external force component
!  (i.e. can be called from externf for doing
!   hydro relaxation runs in the external tokamak potential)
!+
!------------------------------------------------------------
subroutine get_externalB_force(xi,yi,zi,rhoi,fextx,fexty,fextz)
 real, intent(in)  :: xi,yi,zi,rhoi
 real, intent(out) :: fextx,fexty,fextz
 real :: dumx,dumy,dumz,dumgx,dumgy,dumgz

 call externBfield(xi,yi,zi,0.,0.,0.,rhoi, &
                   0.,0.,0.,0.,0.,0.,dumx,dumy,dumz, &
                   fextx,fexty,fextz,dumgx,dumgy,dumgz,'fext')

end subroutine get_externalB_force


!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_externB(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to force from external magnetic field'
 call write_inopt(Bstar_cgs,'Bstar','magnetic field strength in Gauss at stellar surface',iunit)

end subroutine write_options_externB

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externB(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use physcon, only:pi
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_externB'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('Bstar')
    ngot = ngot + 1
    read(valstring,*,iostat=ierr) Bstar_cgs
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_externB

!-----------------------------------------------------------------------
!+
!  checks that input file options are reasonable
!+
!-----------------------------------------------------------------------
subroutine check_externB_settings(ierr)
 use io, only:error
 integer, intent(out) :: ierr

 ierr = 0

end subroutine check_externB_settings

end module extern_Bfield
