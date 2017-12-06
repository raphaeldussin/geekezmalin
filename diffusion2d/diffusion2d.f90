PROGRAM diffusion2d

  !! Author : Raf Dussin

  USE NETCDF ! use netcdf I/O library

  ! old fortran used to assume variables starting with i = integers,...
  ! override this behavior with :
  IMPLICIT NONE

  ! we then need to declare the type of all variables

  !! Time integration
  REAL(KIND=8)                           :: dt                ! time step
  INTEGER                                :: nt                ! number of time steps
  INTEGER                                :: nwrite            ! write every nwrite steps
  INTEGER                                :: kt                ! time loop index

  !! Anomaly and physical parameters
  REAL(KIND=8)                           :: zx0, zy0, r0      ! location/radius of anomaly
  REAL(KIND=8)                           :: kdiff             ! diffusivity

  !! Domain
  REAL(KIND=8)                           :: Lx, Ly            ! horizontal dimensions
  REAL(KIND=8)                           :: dxy               ! horizontal resolution
  INTEGER                                :: nx, ny, nxd, nyd  ! size of domain
  INTEGER                                :: isc, iec,jsc,jec  ! start/end of compute domain
  INTEGER                                :: nhalo=1           ! number of halo points

  !! Namelist reading
  CHARACTER(256)                         :: namlist           ! namelist name
  INTEGER                                :: numnam=20, iost=0 ! namelist I/O
  INTEGER                                :: narg, iargc       ! command line I/O
  LOGICAL                                :: isnam

  !! NETCDF store file, dimensions and variables IDs as integers
  !! and needs a filename and variable name
  CHARACTER(256)                         :: cfileout          ! netcdf file name
  CHARACTER(256)                         :: cvarout='tracer'  ! netcdf variable name
  INTEGER                                :: ncid1             ! netcdf file index
  INTEGER                                :: dimid1, dimid2    ! netcdf dimension index
  INTEGER                                :: dimid3            ! netcdf dimension index
  INTEGER                                :: varid1, varid2    ! netcdf variable index
  INTEGER                                :: varid3, varid4    ! netcdf variable index
  INTEGER                                :: start(3),count(3) ! netcdf size of variable

  ! Dynamically allocated array for 2d variables
  REAL(KIND=8), DIMENSION(:,:), &
              & ALLOCATABLE              :: znew, zold, zrhs
  ! Dynamically allocated array for 1d variables
  REAL(KIND=8), DIMENSION(:), &
              & ALLOCATABLE              :: zx, zy

  INTEGER                                :: ji, jj            ! loop indices

  ! declare what parameters to read from the namelist ( = parameter file )
  NAMELIST/naminput/ Lx, Ly, dxy, dt, nt, nwrite, zx0, zy0, r0, kdiff, cfileout

  !------------------------------------------------------------
  !---- Namelist reading --------------------------------------

  narg = iargc()

  ! If no namelist is passed in argument we stop
  IF (narg /= 1) THEN
      PRINT *, 'this program takes a namelist as single argument' ; STOP
  ELSE
      CALL getarg(1,namlist)
  ENDIF

  PRINT *, 'The namelist used is : ' , namlist

  ! Check existence of the namelist
  INQUIRE(FILE=namlist , EXIST=isnam )

  IF (isnam .eqv. .FALSE.) THEN
     PRINT *, 'Namelist not found' ; STOP
  ENDIF

  ! then Read it
  OPEN( UNIT=numnam, FILE=namlist, FORM='FORMATTED', &
      & ACCESS='SEQUENTIAL', ACTION= 'READ' , STATUS='OLD' , IOSTAT=iost )

  REWIND( numnam )
  READ( numnam, naminput )

  CLOSE(UNIT=numnam)

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Define the grid for the experiment : here we are going to 
  ! work with a bi-periodic grid and put one halo point on each
  ! side so we can compute the 5pt laplacian everywhere

  ! number of points in the computational domain
  nxd = INT(Lx) / INT(dxy)
  nyd = INT(Ly) / INT(dxy)
  ! number of points with halo
  nx = nxd + 2*nhalo
  ny = nyd + 2*nhalo
  ! start/end indices of the computational domain
  isc = 1 + nhalo
  iec = nx - nhalo
  jsc = 1 + nhalo
  jec = ny - nhalo

  PRINT *, 'size of computational domain (nxd,nyd) = ', nxd, nyd
  PRINT *, 'size of domain with halo (nx,ny) = ', nx, ny
  PRINT *, 'computational domain starts/end at i = ', isc,iec
  PRINT *, 'computational domain starts/end at j = ', jsc,jec

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Now that we have the dimensions, we can allocate the arrays for
  ! our computation, and populate them
  ALLOCATE( zx(nx), zy(ny) )
  ALLOCATE( zold(nx,ny) , znew(nx,ny), zrhs(nx,ny) )

  ! create x-coordinate from 0 to Lx
  DO ji=1,nx
     zx(ji) = (ji-1-nhalo) * dxy
  ENDDO

  ! create y-coordinate from 0 to Ly
  DO jj=1,ny
     zy(jj) = (jj-1-nhalo) * dxy
  ENDDO

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Prepare output file

  ! Create file
  CALL CHECK( NF90_CREATE(trim(cfileout), nf90_clobber, ncid1) )

  ! Define dimensions
  CALL CHECK( NF90_DEF_DIM(ncid1, 'x' , nxd,  dimid1) )
  CALL CHECK( NF90_DEF_DIM(ncid1, 'y' , nyd,  dimid2) )
  CALL CHECK( NF90_DEF_DIM(ncid1, 'time' , NF90_UNLIMITED, dimid3) )

  ! Define variables
  CALL CHECK( NF90_DEF_VAR(ncid1, 'x',  NF90_DOUBLE, (/ dimid1 /), varid1) )
  CALL CHECK( NF90_DEF_VAR(ncid1, 'y',  NF90_DOUBLE, (/ dimid2 /), varid2) )
  CALL CHECK( NF90_DEF_VAR(ncid1, 'time',  NF90_DOUBLE, (/ dimid3 /), varid3) )
  CALL CHECK( NF90_DEF_VAR(ncid1, trim(cvarout), NF90_DOUBLE, (/ dimid1, dimid2, dimid3 /), varid4) )

  ! Assign units attributes to coordinate variables.
  CALL CHECK( NF90_PUT_ATT(ncid1, varid1, "units", "km") )
  CALL CHECK( NF90_PUT_ATT(ncid1, varid2, "units", "km") )

  ! End define mode. Now we can write data in the file.
  CALL CHECK( NF90_ENDDEF(ncid1) )

  ! Write coordinates to output file
  CALL CHECK( NF90_PUT_VAR(ncid1, varid1, zx(isc:iec)) )
  CALL CHECK( NF90_PUT_VAR(ncid1, varid2, zy(jsc:jec)) )

  ! Create array for start and count, required to write data
  ! for tracer in the right record with the right size
  start = (/ 1, 1, 1 /)
  count = (/ nxd, nyd, 1 /)

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Initial condition = gaussian centered in (zx0,zy0), radius r0

  DO jj=jsc,jec
     DO ji=isc,iec
        znew(ji,jj) = EXP(- 0.5 * (zx(ji) - zx0)**2 / (r0*r0) ) * &
 &                    EXP(- 0.5 * (zy(jj) - zy0)**2 / (r0*r0) )
     ENDDO
  ENDDO

  ! write initial condition to file
  CALL CHECK( NF90_PUT_VAR(ncid1, varid4, znew(isc:iec,jsc:jec), &
 &                         start = start, count = count) )

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Stability test for Euler Scheme

  IF (kdiff*dt/(dxy*dxy) >= 0.5) THEN
     PRINT *, 'stability condition requires kdiff * dt / dxy*dxy < 0.5 , here : ', kdiff*dt/(dxy*dxy)
     STOP
  ENDIF

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Time Loop

  DO kt=1,nt

     print *, ' time step = ',kt

     ! swap arrays
     zold = znew

     ! update halo values using the bi-periodic properties of the domain
     CALL halo_update(zold,nx,ny,isc,iec,jsc,jec)

     ! compute the diffusive term
     CALL laplacian_5pt(zold,dxy,zrhs,nx,ny,isc,iec,jsc,jec)

     ! integrate forward in time using the Euler time-stepping scheme
     znew = zold + dt * (kdiff * zrhs )

     ! Every nwrite steps, write the state of the model into the output file
     IF (MOD(kt,nwrite) == 0) THEN
        start(3) = start(3) + 1
        CALL CHECK( NF90_PUT_VAR(ncid1, varid4, znew(isc:iec,jsc:jec), &
 &                               start = start, count = count) )
     ENDIF

  ENDDO ! end time loop

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Run is finished, closing netcdf file and cleaning up

  CALL CHECK( NF90_CLOSE(ncid1) )
  DEALLOCATE( zx, zy, zold, znew, zrhs )

CONTAINS

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Laplacian for diffusive term using a 5 points stencil
  SUBROUTINE laplacian_5pt(fieldin,hres,fieldout,nx,ny,isc,iec,jsc,jec)
    IMPLICIT NONE
    REAL(8), DIMENSION(nx,ny), INTENT(in) :: fieldin
    REAL(8), DIMENSION(nx,ny), INTENT(out) :: fieldout
    REAL(8), INTENT(in) :: hres
    INTEGER, INTENT(in) :: nx,ny,isc,iec,jsc,jec
    INTEGER :: ji, jj
    REAL(8) :: scfac

    scfac = (1 / (hres*hres))
    DO jj=jsc,jec
       DO ji=isc,iec
          fieldout(ji,jj) = scfac * (fieldin(ji-1,jj) + fieldin(ji+1,jj) + fieldin(ji,jj-1) + &
  &                                  fieldin(ji,jj+1) - 4*fieldin(ji,jj) )
       ENDDO
    ENDDO

  END SUBROUTINE

  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Update for halos, using bi-periodicity
  SUBROUTINE halo_update(field,nx,ny,isc,iec,jsc,jec)
    IMPLICIT NONE
    REAL(8), DIMENSION(nx,ny), INTENT(inout) :: field
    INTEGER, INTENT(in) :: nx,ny,isc,iec,jsc,jec
    INTEGER :: ji, jj

    DO jj=jsc,jec
       field(1,jj)  = field(iec,jj)
       field(nx,jj) = field(isc,jj)
    ENDDO

    DO ji=isc,iec
       field(ji,1) = field(ji,jec)
       field(ji,ny) = field(ji,jsc)
    ENDDO

    ! corners are not needed by 5pt stencil so we can safely put them to zero
    ! otherwise it would need a more clever way to set the value
    field(1,1) = 0.
    field(1,ny) = 0.
    field(nx,1) = 0.
    field(nx,ny) = 0.

  END SUBROUTINE
  
  !------------------------------------------------------------
  !------------------------------------------------------------
  ! This subroutine returns potential errors in a call to netcdf
  ! library as a human readable error (in english) rather than an
  ! error code (e.g. 34)
  SUBROUTINE check(status)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: status

    IF(status /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(status))
      STOP 2
    ENDIF
  END SUBROUTINE check

END PROGRAM diffusion2d
