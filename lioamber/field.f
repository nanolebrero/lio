            subroutine EFIELD(istep,fxx,fyy,fzz,fieldtype,E1,M11,MM)            
!!!!!!!!  FIELD TYPES: GAUSSIAN,INITSTAT,INITDUMP
            use garcha_mod
            REAL*8,intent(out) :: fxx,fyy,fzz
            REAL*8 :: fac 
            CHARACTER(LEN=8)   :: FIELDTYPE
            REAL*8, intent(out) :: E1
            REAL*8 :: ux,uy,uz
            call dip(ux,uy,uz)
            write(*,*) 'ux,uy,uz =', ux,uy,uz
            if (exter) then
                g=1.0D0
                fac=2.54D0
                E1=0
                IF (FIELDTYPE.EQ.'GAUSSIAN') THEN
                   write(*,*) 'FIELD TYPE: GAUSSIAN'
                   fxx=fx*exp(-0.2*(real(istep-50))**2)
                   fyy=fy*exp(-0.2*(real(istep-50))**2)
                   fzz=fz*exp(-0.2*(real(istep-50))**2)
                ELSE IF (FIELDTYPE.EQ.'INITDUMP') THEN
                   write(*,*) 'FIELD TYPE: INITIALLY POLARIZED - DUMPED
     > DEPOLARIZATION'
                   fxx=fx*exp(-0.00000005*((istep)**2))
                   fyy=fy*exp(-0.00000005*((istep)**2))
                   fzz=fz*exp(-0.00000005*((istep)**2))
                   write(*,*) Fxx,Fyy,Fzz
                ELSE IF (FIELDTYPE.EQ.'INITSTAT') THEN
                   write(*,*) 'FIELD TYPE: INITIALLY POLARIZED'
                   fxx=0
                   fyy=0
                   fzz=0
                ELSE IF (FIELDTYPE.EQ.'SELFCONS') THEN
                   write(*,*) 'FIELD TYPE: SELF CONSISTENT FIELD'      
                   fxx=fx
                   fyy=fy
                   fzz=fz
                ENDIF
            else
                g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
                Fx=ux/2.54D0
                Fy=uy/2.54D0
                Fz=uz/2.54D0
                fac=(2.54D0*2.00D0)
            endif
            call dip2(g,Fxx,Fyy,Fzz)
            write(*,*) 'epsilon =', epsilon
            write(*,*) 'a0 =',a0
            write(*,*) 'fac =', fac
            write(*,*) 'Qc2 =', Qc2
            E1=-1.00D0*g*(Fxx*ux+Fyy*uy+Fzz*uz)/fac -
     >      0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0 
            do k=1,MM
               E1=E1+RMM(k)*RMM(M11+k-1)
            enddo
            write(*,*) 'E1=', E1
            return
            end












