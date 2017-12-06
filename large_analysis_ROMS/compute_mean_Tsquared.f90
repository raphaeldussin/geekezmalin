PROGRAM compute_mean_Tsquared
!------------------------------------------------------
!------------------------------------------------------
!
! Compute the mean of T*T in ROMS output
! Raf Dussin, 2015

USE NETCDF

   IMPLICIT NONE

   INTEGER                                       :: narg, iargc               ! command line arguments
   INTEGER                                       :: ijarg, nfil
   INTEGER                                       :: start(4), count(4)        ! netcdf diemnsions
   ! netcdf file id
   INTEGER                                       :: ncid_in1, ncid_in, ncid_out
   ! netcdf dimensions id for input
   INTEGER                                       :: id_xi_rho, id_xi_u, id_xi_v, id_xi_psi
   INTEGER                                       :: id_eta_rho, id_eta_u, id_eta_v, id_eta_psi
   INTEGER                                       :: id_s_rho, id_s_w
   INTEGER                                       :: id_ocean_time
   ! netcdf dimensions id for output
   INTEGER                                       :: outid_xi_rho, outid_xi_u, outid_xi_v, outid_xi_psi
   INTEGER                                       :: outid_eta_rho, outid_eta_u, outid_eta_v, outid_eta_psi
   INTEGER                                       :: outid_s_rho, outid_s_w
   INTEGER                                       :: outid_ocean_time
   ! netcdf variables ids for input
   INTEGER                                       :: id_var1
   ! netcdf variables ids for output
   INTEGER                                       :: outid_var1
   ! grid dimension
   INTEGER                                       :: Lm, Mm, Lm1, Mm1, Lm2, Mm2, Sr, Sw
   ! loop indexes
   INTEGER                                       :: jf, jk
   REAL(4)                                       :: spval=1.E+37
   ! arrays
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE   :: zT2, zsum
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE   :: zwrk1, zwrk2
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE   :: zmask

   CHARACTER(LEN=1024), DIMENSION(:),ALLOCATABLE :: cf_list                   ! list of input files
   CHARACTER(LEN=1024)                           :: cldum, cf1                ! temp string, name of first file
   CHARACTER(LEN=1024)                           :: cfilein                   ! name of current input file
   CHARACTER(LEN=1024), PARAMETER                :: cfileout="mean_T2.nc"     ! name of output file
   CHARACTER(LEN=64)                             :: cvar_out                  ! name of output variable

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read command line
   narg= iargc()

   IF ( narg == 0 ) THEN
      PRINT *,' Usage : ./compute_mean_Tsquared list_of_model_files       '
      PRINT *,' output in mean_T2.nc                                      '
      PRINT *,' --------------------------------------------------------- ' 
      STOP
   ENDIF

   ALLOCATE ( cf_list(narg) )
   ijarg = 1
   nfil = 0
   DO WHILE ( ijarg <= narg )
      CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
      SELECT CASE ( cldum )
      CASE DEFAULT         ! then the argument is a file
         nfil          = nfil + 1
         cf_list(nfil) = TRIM(cldum)
      END SELECT
   END DO

   PRINT *, 'Number of files to process = ', nfil

   cvar_out = 'T2'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read dimensions and lon/lat/depth from first file

   cf1 = cf_list(1)

   CALL check( NF90_OPEN(trim(cf1), nf90_nowrite, ncid_in1) )

   CALL check( NF90_INQ_DIMID(ncid_in1, 'xi_rho', id_xi_rho) )
   CALL check( NF90_INQ_DIMID(ncid_in1, 'xi_u',   id_xi_u)   )
   CALL check( NF90_INQ_DIMID(ncid_in1, 'xi_v',   id_xi_v)   )
   !CALL check( NF90_INQ_DIMID(ncid_in1, 'xi_psi', id_xi_psi) )

   CALL check( NF90_INQ_DIMID(ncid_in1, 'eta_rho', id_eta_rho) )
   CALL check( NF90_INQ_DIMID(ncid_in1, 'eta_u',   id_eta_u)   )
   CALL check( NF90_INQ_DIMID(ncid_in1, 'eta_v',   id_eta_v)   )
   !CALL check( NF90_INQ_DIMID(ncid_in1, 'eta_psi', id_eta_psi) )

   CALL check( NF90_INQ_DIMID(ncid_in1, 's_rho', id_s_rho) )
   CALL check( NF90_INQ_DIMID(ncid_in1, 's_w',   id_s_w)   )

   CALL check( NF90_INQ_DIMID(ncid_in1, 'ocean_time', id_ocean_time) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_in1, id_xi_rho, len=Lm2) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_in1, id_eta_rho,len=Mm2) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_in1, id_s_rho,  len=Sr)  )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_in1, id_s_w,    len=Sw)  )

   CALL check( NF90_CLOSE(ncid_in1) )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Compute array sizes and allocate them

   Lm = Lm2-2 ; Lm1 = Lm+1
   Mm = Mm2-2 ; Mm1 = Mm+1

   PRINT *, 'Dimensions are Lm2 x Mm2 x Sr = ', Lm2, Mm2, Sr
   ALLOCATE( zT2(Lm2,Mm2,Sr), zsum(Lm2,Mm2,Sr) )

   ALLOCATE( zwrk1(Lm2,Mm2) , zwrk2(Lm2,Mm2) )
   ALLOCATE( zmask(Lm2,Mm2) )

   count = (/ Lm2, Mm2, 1, 1 /)
   start = (/ 1, 1, 1, 1 /)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Reopen first file to get the mask

   CALL check( NF90_OPEN(trim(cf1), nf90_nowrite, ncid_in1) )
   ! define land sea mask
   zmask(:,:) = 1.0
   CALL check( NF90_INQ_VARID(ncid_in1, 'temp', id_var1) )
   CALL check( NF90_GET_VAR(ncid_in1, id_var1, zwrk1, start = start, &
        &                 count = count) )

   WHERE ( zwrk1(:,:) == spval ) zmask(:,:) = 0.0
   CALL check( NF90_CLOSE(ncid_in1) )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  create the netcdf output file

   CALL check( nf90_create(trim(cfileout), nf90_clobber, ncid_out) )

   CALL check( nf90_def_dim(ncid_out, 'xi_rho',  Lm2, outid_xi_rho)  )
   CALL check( nf90_def_dim(ncid_out, 'eta_rho', Mm2, outid_eta_rho) )
   CALL check( nf90_def_dim(ncid_out, 'xi_u',  Lm1, outid_xi_u)  )
   CALL check( nf90_def_dim(ncid_out, 'eta_u', Mm2, outid_eta_u) )
   CALL check( nf90_def_dim(ncid_out, 'xi_v',  Lm2, outid_xi_v)  )
   CALL check( nf90_def_dim(ncid_out, 'eta_v', Mm1, outid_eta_v) )
   !CALL check( nf90_def_dim(ncid_out, 'xi_psi',  Lm1, outid_xi_psi)  )
   !CALL check( nf90_def_dim(ncid_out, 'eta_psi', Mm1, outid_eta_psi) )
   CALL check( nf90_def_dim(ncid_out, 's_rho',   Sr,  outid_s_rho) )
   CALL check( nf90_def_dim(ncid_out, 's_w',     Sw,  outid_s_w)   )
   CALL check( nf90_def_dim(ncid_out, 'ocean_time', NF90_UNLIMITED, outid_ocean_time) )

   CALL check( nf90_def_var(ncid_out, trim(cvar_out), NF90_FLOAT, &
   & (/ outid_xi_rho, outid_eta_rho, outid_s_rho, outid_ocean_time /), outid_var1) )

   ! Assign units attributes to coordinate variables.
   CALL check( nf90_put_att(ncid_out, outid_var1, "units", "degC.degC") )
   CALL check( nf90_put_att(ncid_out, outid_var1, "_FillValue", spval) )

   ! End define mode.
   CALL check( nf90_enddef(ncid_out) )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Initialize

   zT2(:,:,:)  = 0.0d0
   zsum(:,:,:) = 0.0d0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Loop on files

   DO jf=1,nfil

      cfilein = cf_list(jf)

      PRINT *, 'Processing file ', trim(cfilein)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Check existence of variable in the current input netcdf file

      CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_in) )
      CALL check( NF90_INQ_VARID(ncid_in, 'temp', id_var1) )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Loop on vertical coordinate (fortran is faster with 2d arrays than 3d)
      DO jk=1,Sr

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! read arrays from nc file

         ! read temperature at level jk
         start(3) = jk
         CALL check( NF90_GET_VAR(ncid_in, id_var1, zwrk1, start = start, &
                  &                 count = count) )

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! compute some derived quantity (here T*T)

         zwrk2(:,:) = zwrk1(:,:) * zwrk1(:,:) * zmask(:,:)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! add to a sum array
         zsum(:,:,jk) = zsum(:,:,jk) + zwrk2(:,:)

      ENDDO ! z loop

   ! close current input file
   CALL check( NF90_CLOSE(ncid_in) )

   ENDDO ! loop on files

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! loop on files done now compute mean
   !! divide by number of files -> mean

   zT2(:,:,:) = zsum(:,:,:) / nfil

   !! mask final array
   WHERE (zsum(:,:,:) == 0.0) zT2(:,:,:) = spval

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! write 3d final array to file, loop on levels because ifort makes a seg
   !! fault when dumping all 3d array

   DO jk=1,Sr

      count = (/ Lm2, Mm2, 1, 1 /)
      start = (/ 1, 1, jk, 1 /)

   CALL check( NF90_PUT_VAR(ncid_out, outid_var1, REAL(zT2(:,:,jk)), start = start, &
               &            count = count) )

   ENDDO
   ! close output file
   CALL check( nf90_close(ncid_out) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS

      SUBROUTINE check(status)
        INTEGER, INTENT(in) :: status
    
        IF(status /= nf90_noerr) THEN
          PRINT *, trim(nf90_strerror(status))
          STOP 2
        ENDIF
      END SUBROUTINE check

END PROGRAM
