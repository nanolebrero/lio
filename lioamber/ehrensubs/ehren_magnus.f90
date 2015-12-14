!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  function ehren_magnus(Nsize,Norder,Fmat,Rold,dt) &
  result(Rnew)
!
! Fmat,Rold,Rnew => All in the ON basis
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)    :: Nsize,Norder
  real*8,intent(in)     :: Fmat(Nsize,Nsize)
  complex*16,intent(in) :: Rold(Nsize,Nsize)
  real*8,intent(in)     :: dt
  real*8                :: Rnew(Nsize,Nsize)

  integer :: kk
  real*8  :: factinv
  complex*16,allocatable :: Omega1(:,:)
  complex*16,allocatable :: ConmutAcum(:,:)
  complex*16,allocatable :: TermPos(:,:),TermNeg(:,:)


! Initializations
!------------------------------------------------------------------------------!
  allocate(Omega1(Nsize,Nsize))
  allocate(ConmutAcum(Nsize,Nsize))
  allocate(TermPos(Nsize,Nsize),TermNeg(Nsize,Nsize))
  Rnew=Rold
  ConmutAcum=Rold
  Omega1=DCMPLX(0.0d0,-1.0d0)*(Fmat)*(dt)
  factinv=1.0d0


! Calculations
!------------------------------------------------------------------------------!
  do kk=1,Norder
    TermPos=matmul(Omega1,ConmutAcum)
    TermNeg=matmul(ConmutAcum,Omega1)
    ConmutAcum=TermPos-TermNeg
    factinv=factinv/dble(kk)
    Rnew=Rnew+factinv*ConmutAcum
  enddo


  deallocate(Omega1,ConmutAcum,TermPos,TermNeg)
  return;end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
