!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine fqbias_setup(Natom,Nbase,atom_of_base,file_name,Qpot)
!--------------------------------------------------------------------!
  implicit none
  integer,intent(in)          :: Natom,Nbase
  integer,intent(in)          :: atom_of_base(Nbase)
  character(len=*),intent(in) :: file_name
  real*8,intent(in)           :: Qpot
  integer                     :: kk

  if (allocated(forb_factor)) deallocate(forb_factor)
  allocate(forb_factor(Nbase))

  if (allocated(group_of_base)) deallocate(group_of_base)
  allocate(group_of_base(Nbase))

  if (allocated(group_of_atom)) deallocate(group_of_atom)
  allocate(group_of_atom(Natom))



  call read_atom_prop &
  (Natom,group_of_atom,file_name)

  call align_ab_prop &
  (Natom,Nbase,group_of_atom,atom_of_base,group_of_base)



  Nforb=Nbase
  Qpot_0=Qpot
  do kk=1,Nbase
    forb_factor(kk)= 0.0d0
    if (group_of_base(kk).eq.1) forb_factor(kk)=-1.0d0
    if (group_of_base(kk).eq.2) forb_factor(kk)= 1.0d0
  enddo



  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
