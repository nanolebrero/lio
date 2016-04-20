!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine fqbias_calc(M,current_time,sqsmat,fockio)
  implicit none
  integer,intent(in)   :: M
  real*8,intent(in)    :: current_time
  real*8,intent(in)    :: sqsmat(M,M)
  real*8,intent(inout) :: fockio(M,M)

  real*8,allocatable  :: forb_qpot(:)
  real*8              :: Qpot_t
  real*8              :: newterm
  integer             :: ii,jj,kk

  call gaussian_shape(current_time,0.0d0,0.7648d0,Qpot_0,Qpot_t)

  allocate(forb_qpot(M))
  do kk=1,M
    forb_qpot(kk)=forb_factor(kk)*Qpot_t
  enddo

  do ii=1,M
  do jj=1,M
     newterm=0.0d0
     do kk=1,M
        newterm=newterm+sqsmat(ii,kk)*forb_qpot(kk)*sqsmat(kk,jj)
     enddo
     fockio(ii,jj)=fockio(ii,jj)+newterm
  enddo
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
