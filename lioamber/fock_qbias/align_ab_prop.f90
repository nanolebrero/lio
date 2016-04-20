!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine &
  align_ab_prop(Natom,Nbase,prop_of_atom,atom_of_base,prop_of_base)
!
! Align atom - base property
!--------------------------------------------------------------------!
  implicit none
  integer,intent(in)  :: Natom,Nbase
  integer,intent(in)  :: prop_of_atom(Natom)
  integer,intent(in)  :: atom_of_base(Nbase)
  integer,intent(out) :: prop_of_base(Nbase)
  integer :: kk

  do kk=1,Nbase
    prop_of_base(kk)=prop_of_atom(atom_of_base(kk))
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
