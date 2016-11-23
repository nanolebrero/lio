!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! TIME
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE calc_current_time(propagator,lpfrg_steps, initial_step
     > , istep, tdstep,t)
       IMPLICIT NONE
       integer, intent(in) :: propagator,lpfrg_steps, initial_step 
       real*8, intent(inout) :: t
       integer, intent(in) :: istep
       real*8, intent(in) :: tdstep 
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
              if (istep.gt.1) then         
              if ((propagator.eq.2) .and. 
     >           (initial_step+istep.le.lpfrg_steps)) then
                 t=t+tdstep*0.002419
              else
                 t=t+tdstep*0.02419
              endif
              write(*,*) 'evolution time (fs)  =', t
              endif
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

