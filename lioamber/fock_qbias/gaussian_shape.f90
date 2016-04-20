!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine gaussian_shape(x,x0,sigma,y0,y)
!--------------------------------------------------------------------!
  implicit none
  real*8,intent(in)  :: x
  real*8,intent(in)  :: x0
  real*8,intent(in)  :: sigma
  real*8,intent(in)  :: y0
  real*8,intent(out) :: y

  real*8 :: expo

  expo=(x-x0)/sigma
  expo=(-1.0)*expo*expo
  y=y0*exp(expo)

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
