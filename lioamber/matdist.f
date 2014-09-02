!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE matdist(NN,MA,MB,dist)
!------------------------------------------------------------------------------!
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: NN
       REAL*8,INTENT(IN)  :: MA(NN,NN),MB(NN,NN)
       REAL*8,INTENT(OUT) :: dist
       INTEGER            :: ii,jj
!------------------------------------------------------------------------------!
       dist=0.0d0
       DO ii=1,NN
       DO jj=1,NN
         dist=dist+(MA(ii,jj)-MB(ii,jj))**2
       ENDDO
       ENDDO
       dist=sqrt(dist)
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
