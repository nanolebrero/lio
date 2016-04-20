!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine fqbias_lowdinpop(M,N,DensMtrx,sqsmat,Qatom,file_name)
!--------------------------------------------------------------------!
  implicit none
  integer,intent(in)          :: M,N
!  complex*16,intent(in)       :: DensMtrx(M,M)
  real*8,intent(in)           :: DensMtrx(M,M)
  real*8,intent(in)           :: sqsmat(M,M)
  integer,intent(in)          :: Qatom(N)
  character(len=*),intent(in) :: file_name

  integer            :: ii,jj,kk,nn
  real*8             :: newterm
  real*8,allocatable :: group_pop(:)
  real*8,allocatable :: DensReal(:,:)


  allocate(DensReal(M,M))
  allocate(group_pop(3))

!  do ii=1,M
!  do jj=1,M
!    DensReal(ii,jj)=REAL(DensMtrx(ii,jj))
!  enddo
!  enddo
  DensReal=DensMtrx

  group_pop=0.0d0

  do kk=1,N
     nn=group_of_atom(kk)
     group_pop(nn)=group_pop(nn)+real(Qatom(kk))
  enddo

  do ii=1,M
  do jj=1,M
  do kk=1,M
     nn=group_of_base(kk)
     newterm=DensReal(ii,jj)*sqsmat(jj,kk)*sqsmat(kk,ii)
     group_pop(nn)=group_pop(nn)-newterm
  enddo
  enddo
  enddo
  

! Output
!--------------------------------------------------------------------!
  open(unit=1001,file=file_name,position='APPEND')
  write(unit=1001,fmt=100) ' Group      Population'
  do kk=1,3
    write(unit=1001,fmt=101) kk,group_pop(kk)
  enddo
  write(unit=1001,fmt=100) ''
  close(unit=1001)
100 format(A)
101 format(I5,4x,f16.8)

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
