!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine rmmgetP(M,Pmat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in) :: M
       real*8,intent(out) :: Pmat(M,M)
       integer            :: rmmp0,idx0,ii,jj
!
       rmmp0=0
       do jj=1,M
         idx0=(2*M-jj)*(jj-1)/2
         Pmat(jj,jj)=RMM(rmmp0+idx0+jj)
         if (jj.lt.M) then
         do ii=jj+1,M
           Pmat(ii,jj)=RMM(rmmp0+idx0+ii)*(0.5d0)
           Pmat(jj,ii)=RMM(rmmp0+idx0+ii)*(0.5d0)
         enddo
         endif
       enddo

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine rmmgetSF(M,SFmat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in) :: M
       real*8,intent(out) :: SFmat(M,M)
       integer            :: rmmp0,idx0,ii,jj
!
       rmmp0=M*(M+1)
       do jj=1,M
         idx0=(2*M-jj)*(jj-1)/2
         do ii=jj,M
           SFmat(ii,jj)=RMM(rmmp0+idx0+ii)
           SFmat(jj,ii)=RMM(rmmp0+idx0+ii)
         enddo
       enddo
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
