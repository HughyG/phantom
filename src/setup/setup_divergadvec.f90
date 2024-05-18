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
!   - iBfield : * which external field to use*
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
 integer :: iBfield

 integer, parameter :: maxwaves = 3
 character(len=*), parameter :: wavetype(0:maxwaves-1) = &
      (/'no cleaning                                ', &
        'hyperbolic cleaning                        ', &
        'hyperbolic/parabolic cleaning              '/)
 
 integer, parameter :: maxBtype = 2
 character(len=*), parameter :: Btype(0:maxBtype-1) = &
      (/'Internal B fiel                           ', &
         'External B field                          '/)

 private

contains

!----------------------------------------------------------------
!+
!  setup for MHD wave tests
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use dim,            only:maxvxyzu,periodic
 use setup_params,   only:rhozero,ihavesetupB
 use boundary,       only:dxbound,dybound,dzbound
 use options,        only:iexternalforce
 use cons2prim,      only:cons2prim_everything
 use externalforces, only:iext_externB
 use part,           only:Bxyz,mhd,periodic,igas
 use io,             only:master,fatal
 use slab,           only:set_slab
 use prompting,      only:prompt
 use mpiutils,       only:bcast_mpi
 use physcon,        only:pi
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
 integer :: i,ierr
 real :: przero,uuzero,Bvec(3),Bzero(3),vzero(3)
 real :: uui,Bxi,r0,rad
 real :: gam1,scaleFactor
 real(kind=4) :: dumdvdx(4,1),dumalphaind(4,1)
 real :: dumrad(4,1),dumeos(4,1),dumrprop(4,1),dumbevol(4,1)
 real :: dumdustevol(4,1),dumdustfrac(4,1)
 character(len=len(fileprefix)+6) :: setupfile

 if (.not.periodic) call fatal('setup_divBadvect','need to compile with PERIODIC=yes')
!
!--general parameters
!
 time = 0.
 gamma = 5./3.
!
!--default settings
!
 iselect = 0
 iBfield = 0
 r0 = 0.25
 !
 ! read setup parameters from the .setup file.
 ! if file does not exist, then ask for user input
 !
 setupfile = trim(fileprefix)//'.setup'
 call read_setupfile(setupfile,ierr)
 if (ierr /= 0) then
    if (id==master) then
       call interactive_setup()
       call write_setupfile(setupfile)
       print*,' Edit '//trim(setupfile)//' and rerun phantomsetup'
    endif
    stop
 endif
!
!--setup parameters
!
 r0 = 0.25
 rhozero  = 1.
 Bzero = 0.
 gam1 = gamma - 1.
 scaleFactor = 1./sqrt(4*pi)
 przero = 6
 vzero  = (/1,1,0/)
 uuzero = przero/(gam1*rhozero)

 if (maxvxyzu < 4) then
    polyk = przero/rhozero**gamma
 else
    polyk = 0.
 endif

 call bcast_mpi(nx)
!
!--boundaries
!
 call set_slab(id,master,nx,-2.,3.,-2.,3.,deltax,hfact,npart,xyzh)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(igas)

 if (iBfield == 1) then
   iexternalforce = iext_externB
 call cons2prim_everything(npart,xyzh,vxyzu,dumdvdx,dumrad,dumeos, &
         dumrprop,Bxyz(1:3,:),dumbevol,dumdustevol,dumdustfrac,dumalphaind)
   iexternalforce = 0
 endif

 do i=1,npart
    rad = sqrt(xyzh(1,i)**2+xyzh(2,i)**2)
    if (rad<r0) then
       Bxi = scaleFactor*((rad/r0)**8 - 2*(rad/r0)**4 + 1)
    else
       Bxi = 0.
    endif
    vxyzu(1,i) = 1
    vxyzu(2,i) = 1
    vxyzu(3,i) = 0
    if (iBfield == 0) then
       Bvec = Bzero + (/Bxi,0.,scaleFactor/)
       Bxyz(1:3,i) = Bvec
    endif
    uui  = uuzero

    if (maxvxyzu >= 4) vxyzu(4,i) = uui
 enddo

 if (mhd) ihavesetupB = .true.

 tmax = 1.
 dtmax = 0.1*tmax

end subroutine setpart

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

 print "(5(/,i2,' : ',a))",(i,trim(Btype(i)),i=0,maxBtype-1)
 call prompt('Select which problem to run ',iBfield,0,maxBtype-1)

 nx = 50
 call prompt('Enter resolution (number of particles in x)',nx,8)

end subroutine interactive_setup

!------------------------------------------
!+
!  Write setup parameters to input file
!+
!------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use dim,          only:tagline,maxvxyzu
 character(len=*), intent(in) :: filename
 integer,          parameter  :: lu = 20
 integer                      :: ierr1

 write(*,"(a)") ' Writing '//trim(filename)//' with setup info'
 open(unit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom div B advection test'

 write(lu,"(/,a)") '# div B advection tests'
 call write_inopt(iselect,'iselect',' which test to run',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing iselect'

 write(lu,"(/,a)") '# B field type'
 call write_inopt(iBfield,'iBfield',' which external field to use',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing iBfield'

 write(lu,"(/,a)") '# resolution'
 call write_inopt(nx,'nx','resolution (number of particles in x) for -xleft < x < xshock',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing nx'

 close(unit=lu)

end subroutine write_setupfile

!------------------------------------------
!+
!  Read setup parameters from input file
!+
!------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'setup_wave: reading setup options from ',trim(filename)

 nerr = 0
 call read_inopt(nx,'nx',db,min=8,errcount=nerr)
 call read_inopt(iselect,'iselect',db,min=0,errcount=nerr)
 call read_inopt(iBfield,'iBfield',db,min=0,errcount=nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_divBadvect: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup
