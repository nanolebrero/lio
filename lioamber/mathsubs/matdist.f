!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function matdist_dd(MA,MB)
     > result(dist)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       real*8                 :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj

       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!--------------------------------------------------------------------!
       function matdist_zd(MA,MB)
     > result(dist)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*16             :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj

       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!--------------------------------------------------------------------!
       function matdist_dz(MA,MB)
     > result(dist)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       complex*16,intent(in)  :: MB(:,:)
       complex*16             :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj

       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!--------------------------------------------------------------------!
       function matdist_zz(MA,MB)
     > result(dist)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       complex*16,intent(in)  :: MB(:,:)
       complex*16             :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj
       
       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!--------------------------------------------------------------------!
       function matdist_cd(MA,MB)
     > result(dist)
       implicit none
       complex*8,intent(in)   :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*8              :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj

       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!--------------------------------------------------------------------!
       function matdist_dc(MA,MB)
     > result(dist)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       complex*8,intent(in)   :: MB(:,:)
       complex*8              :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj

       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!--------------------------------------------------------------------!
       function matdist_cc(MA,MB)
     > result(dist)
       implicit none
       complex*8,intent(in)   :: MA(:,:)
       complex*8,intent(in)   :: MB(:,:)
       complex*8              :: diff
       real*8                 :: dist
       integer                :: ni,nj,ii,jj

       ni=min(size(MA,1),size(MB,1))
       nj=min(size(MA,2),size(MB,2))
       dist=0.0d0
       do jj=1,nj
       do ii=1,ni
         diff=MA(ii,jj)-MB(ii,jj)
         dist=dist+(abs(diff))**2
       enddo
       enddo
       dist=sqrt(dist)
       return;end function
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
