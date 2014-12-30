!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine rmmputP(M,Pmat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in) :: M
       real*8,intent(in)  :: Pmat(M,M)
       integer            :: rmmp0,idx0,ii,jj
!
       rmmp0=0
       do jj=1,M
         idx0=(2*M-jj)*(jj-1)/2
         RMM(rmmp0+idx0+jj)=Pmat(jj,jj)
         if (jj.lt.M) then
         do ii=jj+1,M
           RMM(rmmp0+idx0+ii)=Pmat(ii,jj)*(2.0d0)
         enddo
         endif
       enddo
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine rmmputSF(M,SFmat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in) :: M
       real*8,intent(in)  :: SFmat(M,M)
       integer            :: rmmp0,idx0,ii,jj
!
       rmmp0=M*(M+1)
       do jj=1,M
         idx0=(2*M-jj)*(jj-1)/2
         do ii=jj,M
           RMM(rmmp0+idx0+ii)=SFmat(ii,jj)
         enddo
       enddo
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
