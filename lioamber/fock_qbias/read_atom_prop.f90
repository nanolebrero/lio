!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine read_atom_prop(Natom,prop_of_atom,file_name)
!--------------------------------------------------------------------!
  implicit none
  integer,intent(in)          :: Natom
  integer,intent(out)         :: prop_of_atom(Natom)
  character(len=*),intent(in) :: file_name
  integer                     :: kk

  open(unit=1001,file=file_name)
  do kk=1,Natom
    read(unit=1001,fmt=*) prop_of_atom(kk)
  enddo
  close(unit=1001)

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
