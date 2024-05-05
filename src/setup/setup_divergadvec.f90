!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for divergence advection tests in 3D
!
! :References:
!    Tricco, T. S. and D. J. Price: 2012, J. Comp. Phys. 231, 7214–7236
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - gamma   : *adiabatic index*
!   - iselect : * which wave test to run*
!   - nx      : *resolution (number of particles in x) for -xleft < x < xshock*
!   - r0      : *radius of the perturbation*
!
! :Dependencies:
!
 implicit none
 public :: setpart
!
! runtime options and default settings
!
 integer :: iselect
 integer :: nx

 integer, parameter :: maxwaves = 3
 character(len=*), parameter :: wavetype(0:maxwaves-1) = &
      (/'no cleaning                                ', &
        'hyperbolic cleaning                        ', &
        'hyperbolic/parabolic cleaning              '/)

 public :: set_perturbation,cons_to_prim! to avoid compiler warnings

 private

contains

!----------------------------------------------------------------
!+
!  setup for MHD wave tests
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use dim,            only:maxvxyzu
 use setup_params,   only:rhozero,ihavesetupB
 use unifdis,        only:set_unifdis,rho_func
 use boundary,       only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,&
                          dxbound,dybound,dzbound
 use options,        only:iexternalforce
 use externalforces, only:iext_externB
 use part,           only:Bxyz,mhd,periodic,igas
 use io,             only:master
 use slab,         only:set_slab
 use prompting,      only:prompt
 use mpiutils,       only:bcast_mpi
 use physcon,        only:pi
 use geometry,       only:igeom_rotated,igeom_cartesian,&
                          set_rotation_angles,coord_transform
 use timestep,       only:tmax,dtmax
 use mpidomain,      only:i_belong
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 real :: deltax,totmass
 integer :: i,ierr,igeom
 real :: przero,uuzero,Bvec(3),vvec(3),Bzero(3),vzero(3)
 real :: uui,Bxi,r0
 real :: gam1,scaleFactor
 real :: drho,dv(3),dB(3),du
 character(len=len(fileprefix)+6) :: setupfile
!
!--general parameters
!
 time = 0.
 gamma = 5./3.
!
!--default settings
!
 iselect = 0
 r0 = 0.25
 !
 ! read setup parameters from the .setup file.
 ! if file does not exist, then ask for user input
 !
 setupfile = trim(fileprefix)//'.setup'
 call read_setupfile(setupfile,gamma,ierr)
 if (ierr /= 0) then
    if (id==master) then
       call interactive_setup()
       call write_setupfile(setupfile,gamma)
       print*,' Edit '//trim(setupfile)//' and rerun phantomsetup'
    endif
    stop
 endif
 gamma = 5./3.
!
!--setup parameters
!
 igeom = igeom_cartesian
 r0 = 0.25
 rhozero  = 1.
 du = 0.
 drho = 0.
 dv = 0.
 dB = 0.
 gam1 = gamma - 1.
 scaleFactor = 1./sqrt(4*pi)
 przero = 6
 vzero  = (/1,1,0/)
 uuzero = przero/(gam1*rhozero)

 call print_amplitudes(rhozero,drho,vzero,dv,Bzero,dB,uuzero,du)

 if (maxvxyzu < 4) then
    polyk = przero/rhozero**gamma
 else
    polyk = 0.
 endif

 call bcast_mpi(nx)
!
!--boundaries
!
!  call set_slab(id,master,nx,-0.75,1.5,-0.75,1.5,deltax,hfact,npart,xyzh)

 call set_boundary(-0.75,1.5,-0.75,1.5,-0.75,0.75)
 deltax = dxbound/nx
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,.false.,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(igas)


 do i=1,npart
    call set_perturbation(xyzh(1,i),xyzh(2,i),r0,scaleFactor,Bxi)
    
    vxyzu(1,i) = 1
    vxyzu(2,i) = 1
    vxyzu(3,i) = 0
    bvec = Bzero + (/Bxi,0.,scaleFactor/)
    uui  = uuzero
  

    if (maxvxyzu >= 4) vxyzu(4,i) = uui
    if (mhd) Bxyz(1:3,i) = Bvec
 enddo

 if (mhd) ihavesetupB = .true.

 tmax = 1.
 dtmax = 0.1*tmax

contains

end subroutine setpart

!-------------------------------------------------
!+
!  simple stretchmapping routine to implement
!  sinusoidal density perturbation
!+
!-------------------------------------------------
subroutine set_perturbation(xi,yi,r0,scaleFactor,Bx)
 real, intent(in)  :: xi,yi,scaleFactor,r0
 real, intent(out) :: Bx
 real    :: rad

 rad = sqrt(xi**2+yi**2)

 if (rad<r0) then
   Bx = scaleFactor*((rad/r0)**8 - 2*(rad/r0)**4 + 1)
 else
   Bx = 0.
 endif

end subroutine set_perturbation


!-----------------------------------------------------
!+
!  nice printout of amplitudes for various quantities
!+
!-----------------------------------------------------
subroutine print_amplitudes(rho,drho,v,dv,B,dB,u,du)
 real, intent(in) :: rho,drho,v(3),dv(3),B(3),dB(3),u,du

 write(*,"('    |',8(a10,1x,'|'))") 'rho','v1','v2','v3','B1','B2','B3','u'
 write(*,"(' q0 |',8(1pg10.2,1x,'|'))") rho,v,B,u
 write(*,"(' dq |',8(1pg10.2,1x,'|'))") drho,dv,dB,du

end subroutine print_amplitudes

!------------------------------------------------
!+
!  conservative to primitive variable transform
!+
!------------------------------------------------
pure subroutine cons_to_prim(q,rho,v,B,u)
 real, intent(in)  :: q(8)
 real, intent(out) :: rho,v(3),B(3),u

 rho = q(1)
 v   = q(2:4)/rho
 B   = q(5:7)
 u   = q(8)/rho - 0.5*dot_product(v,v)

end subroutine cons_to_prim

!------------------------------------------
!+
!  Prompt user for setup options
!+
!------------------------------------------
subroutine interactive_setup()
 use prompting, only:prompt
 integer :: i

 print "(5(/,i2,' : ',a))",(i,trim(wavetype(i)),i=0,maxwaves-1)
 call prompt('Select which problem to run ',iselect,0,maxwaves-1)

 nx = 50
 call prompt('Enter resolution (number of particles in x)',nx,8)

end subroutine interactive_setup

!------------------------------------------
!+
!  Write setup parameters to input file
!+
!------------------------------------------
subroutine write_setupfile(filename,gamma)
 use infile_utils, only:write_inopt
 use dim,          only:tagline,maxvxyzu
 character(len=*), intent(in) :: filename
 real,             intent(in) :: gamma
 integer,          parameter  :: lu = 20
 integer                      :: ierr1

 write(*,"(a)") ' Writing '//trim(filename)//' with setup info'
 open(unit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom MHD linear wave test setup'

 write(lu,"(/,a)") '# MHD wave tests'
 call write_inopt(iselect,'iselect',' which wave test to run',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing iselect'

 write(lu,"(/,a)") '# resolution'
 call write_inopt(nx,'nx','resolution (number of particles in x) for -xleft < x < xshock',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing nx'

 write(lu,"(/,a)") '# Equation-of-state properties'
 call write_inopt(gamma,'gamma','adiabatic index',lu,ierr1)

 close(unit=lu)

end subroutine write_setupfile

!------------------------------------------
!+
!  Read setup parameters from input file
!+
!------------------------------------------
subroutine read_setupfile(filename,gamma,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'setup_wave: reading setup options from ',trim(filename)

 nerr = 0
 call read_inopt(nx,'nx',db,min=8,errcount=nerr)
 call read_inopt(iselect,'iselect',db,min=0,errcount=nerr)
 call read_inopt(gamma,'gamma',db,min=1.,errcount=nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_wave: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup
