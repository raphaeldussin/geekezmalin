!! Water masses computations with f2py
!! (Raf Dussin)

  SUBROUTINE volume_watermass_from_ts(dx,dy,thickness,temp,salt,tmin,tmax,smin,smax,volume)

  IMPLICIT NONE

  REAL(4), DIMENSION(:,:)      :: dx, dy
  REAL(4), DIMENSION(:,:,:)    :: thickness
  REAL(4), DIMENSION(:,:,:)    :: temp,salt
  REAL(4)                      :: tmin,tmax,smin,smax
  INTEGER                      :: nx,ny,nz
  REAL(8)                      :: volume
!f2py INTENT(in)               :: dx, dy, thickness, temp,salt
!f2py INTENT(in)               :: tmin,tmax,smin,smax
!f2py INTENT(out)              :: volume
  INTEGER                      :: ji,jj,jk

  nx = SIZE(temp,1)
  ny = SIZE(temp,2)
  nz = SIZE(temp,3)

  volume=0.0

  DO jk=1,nz
    DO jj=1,ny
       DO ji=1,nx

         IF ( temp(ji,jj,jk) .ge. tmin ) THEN
           IF ( temp(ji,jj,jk) .le. tmax ) THEN
             IF (salt(ji,jj,jk) .ge. smin ) THEN
               IF ( salt(ji,jj,jk) .le. smax ) THEN

                  volume = volume + thickness(ji,jj,jk) * dy(ji,jj) * dx(ji,jj)

               ENDIF
             ENDIF
           ENDIF
         ENDIF

      ENDDO
    ENDDO
  ENDDO

  END SUBROUTINE

  SUBROUTINE volume_watermass_from_ts_v2(dx,dy,thickness,temp,salt,tmin,tmax,smin,smax,nx,ny,nz,volume)

  IMPLICIT NONE

  REAL(4), DIMENSION(nx,ny)    :: dx, dy
  REAL(4), DIMENSION(nx,ny,nz) :: thickness
  REAL(4), DIMENSION(nx,ny,nz) :: temp,salt
  REAL(4)                      :: tmin,tmax,smin,smax
  INTEGER                      :: nx,ny,nz
  REAL(8)                      :: volume
!f2py INTENT(in)               :: dx, dy, thickness, temp,salt
!f2py INTENT(in)               :: tmin,tmax,smin,smax
!f2py INTENT(in)               :: nx,ny,nz
!f2py INTENT(out)              :: volume
  INTEGER                      :: ji,jj,jk

  volume=0.0

  DO jk=1,nz
    DO jj=1,ny
       DO ji=1,nx

         IF ( temp(ji,jj,jk) .ge. tmin ) THEN
           IF ( temp(ji,jj,jk) .le. tmax ) THEN
             IF (salt(ji,jj,jk) .ge. smin ) THEN
               IF ( salt(ji,jj,jk) .le. smax ) THEN

                  volume = volume + thickness(ji,jj,jk) * dy(ji,jj) * dx(ji,jj)

               ENDIF
             ENDIF
           ENDIF
         ENDIF

      ENDDO
    ENDDO
  ENDDO

  END SUBROUTINE

  SUBROUTINE volume_watermass_from_ts_v3(dx,dy,thickness,temp,salt,tmin,tmax,smin,smax,nx,ny,nz,volume)

  IMPLICIT NONE

  REAL(4), DIMENSION(ny,nx)    :: dx, dy
  REAL(4), DIMENSION(nz,ny,nx) :: thickness
  REAL(4), DIMENSION(nz,ny,nx) :: temp,salt
  REAL(4)                      :: tmin,tmax,smin,smax
  INTEGER                      :: nx,ny,nz
  REAL(8)                      :: volume
!f2py INTENT(in)               :: dx, dy, thickness, temp,salt
!f2py INTENT(in)               :: tmin,tmax,smin,smax
!f2py INTENT(in)               :: nx,ny,nz
!f2py INTENT(out)              :: volume
  INTEGER                      :: ji,jj,jk

  volume=0.0

  DO ji=1,nx
    DO jj=1,ny
       DO jk=1,nz

         IF ( temp(jk,jj,ji) .ge. tmin ) THEN
           IF ( temp(jk,jj,ji) .le. tmax ) THEN
             IF (salt(jk,jj,ji) .ge. smin ) THEN
               IF ( salt(jk,jj,ji) .le. smax ) THEN

                  volume = volume + thickness(jk,jj,ji) * dy(jj,ji) * dx(jj,ji)

               ENDIF
             ENDIF
           ENDIF
         ENDIF

      ENDDO
    ENDDO
  ENDDO

  END SUBROUTINE
