!! Water masses computations with f2py
!! (Raf Dussin)

  SUBROUTINE volume_watermass_from_ts(temp,salt,thickness,dx,dy,tmin,tmax,smin,smax,nx,ny,nz,volume)

  IMPLICIT NONE

  REAL(4), DIMENSION(nz,ny,nx), INTENT(in) :: temp,salt
  REAL(4), DIMENSION(nz,ny,nx), INTENT(in) :: thickness
  REAL(4), DIMENSION(ny,nx), INTENT(in)    :: dx, dy
  REAL(4), INTENT(in)                      :: tmin,tmax,smin,smax
  INTEGER, INTENT(in)                      :: nx,ny,nz
  REAL(8), INTENT(out)                     :: volume

  INTEGER                                  :: ji,jj,jk

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

  SUBROUTINE volume_watermass_from_rho(rho,thickness,dx,dy,rhomin,rhomax,nx,ny,nz,volume)

  IMPLICIT NONE

  REAL(4), DIMENSION(nz,ny,nx), INTENT(in) :: rho
  REAL(4), DIMENSION(nz,ny,nx), INTENT(in) :: thickness
  REAL(4), DIMENSION(ny,nx), INTENT(in)    :: dx, dy
  REAL(4), INTENT(in)                      :: rhomin,rhomax
  INTEGER, INTENT(in)                      :: nx,ny,nz
  REAL(8), INTENT(out)                     :: volume

  INTEGER                                  :: ji,jj,jk

  volume=0.0

  DO ji=1,nx
    DO jj=1,ny
      DO jk=1,nz
 
         IF ( rho(jk,jj,ji) .ge. rhomin ) THEN
           IF ( rho(jk,jj,ji) .le. rhomax ) THEN

               volume = volume + thickness(jk,jj,ji) * dy(jj,ji) * dx(jj,ji)

           ENDIF
         ENDIF

      ENDDO
    ENDDO
  ENDDO

  END SUBROUTINE
