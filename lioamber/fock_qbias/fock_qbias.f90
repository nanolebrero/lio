!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  module fock_qbias
!--------------------------------------------------------------------!
!
!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
  implicit none
  integer             :: Nforb
  real*8              :: Qpot_0
  integer,allocatable :: group_of_base(:)
  integer,allocatable :: group_of_atom(:)
  real*8,allocatable  :: forb_factor(:)
  contains
!
!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
  include 'fqbias_setup.f90'
  include 'read_atom_prop.f90'
  include 'align_ab_prop.f90'

  include 'fqbias_calc.f90'
  include 'gaussian_shape.f90'

  include 'fqbias_lowdinpop.f90'
  end module
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
