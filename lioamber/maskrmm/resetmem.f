!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE resetmem()
       use garcha_mod, only:kkind,kkinds,cool,cools
       implicit none
!------------------------------------------------------------------------------!
       if (allocated(kkind))  deallocate(kkind)
       if (allocated(kkinds)) deallocate(kkinds)
       if (allocated(cool))   deallocate(cool)
       if (allocated(cools))  deallocate(cools)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
