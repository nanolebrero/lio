!ppl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE TDop()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c REAL TIME-TDDFT
c 
c Dario Estrin, 1992
c Nano, Dario, Uriel, Damian 2012
c
c  This subrutine takes the converged density matrix from an SCF calculation
c  and evolves it in time. In the input file the total number of propagation
c  steps is specified (nstep) as well as the time of each evolution step 
c  (tdstep). 
c  This implementation has two alternatives to evolve the density in time. The 
c  first one (propagator=1) is the Verlet algorithm that uses a convination of 
c  Liouville von Newmann expresion for the time derivative of the density matrix 
c  and a first order Taylor expansion of the density matrix. The second one 
c  (propagator=2) is the Magnus propagation scheme that uses Backer Campbell
c  Hausdorff (BCH) formula. For this reason when Magnus is used the number of 
c  total conmutators in the BCH espansion has to be specified (NBCH, default=10). 
c  A narrow gaussian type electric field can be introduced during the time 
c  evolution in order to excite all electronic frequencies with the same intensity.
c  Once this perturbation is turned on (Field=t, exter=t) each component of the
c  external electric field has to be specified in the input file (Fx,Fy,Fz).
c  In each step of the propagation the cartesian components of the sistem's dipole
c  are stored in files x.dip, y.dip, z.dip.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c       USE latom
       USE garcha_mod
       USE mathsubs
       use general_module ! FFR-fock_qbias
       use fock_qbias ! FFR-fock_qbias
#ifdef CUBLAS
        use cublasmath
#endif
       IMPLICIT REAL*8 (a-h,o-z)

       INTEGER :: istep
       REAL*8 :: t,E2
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: 
     >   xnano2,Y,fock_a,fock_b,
     >   F1a_a,F1a_b,F1b_a,F1b_b,overlap,rhoscratch
       real*8, dimension (:,:), ALLOCATABLE :: elmu
       DIMENSION q(natom)
       REAL*8,dimension(:),ALLOCATABLE :: factorial
#ifdef TD_SIMPLE
       COMPLEX*8 :: Im,Ix
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) ::
     >   rho_a,rho_b,rhonew_a,rhonew_b,rhold_a,rhold_b,xnano,
     >   rho1
#else
       COMPLEX*16 :: Im,Ix
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) ::
     >   rho_a,rho_b,rhonew_a,rhonew_b,rhold_a,rhold_b,xnano,
     >   rho1
#endif
!!------------------------------------!!
!! FFR ADD
       INTEGER ::
     >   pert_steps,lpfrg_steps,chkpntF1a,chkpntF1b
       REAL*8 ::
     >   dt_magnus,dt_lpfrg
!! CUBLAS
#ifdef CUBLAS
      integer sizeof_real
      parameter(sizeof_real=8)
      integer sizeof_complex
#ifdef TD_SIMPLE
      parameter(sizeof_complex=8)
#else
      parameter(sizeof_complex=16)
#endif
      integer stat
      integer*8 devPtrX, devPtrY,devPtrXc
      external CUBLAS_INIT, CUBLAS_SET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX
#endif
!!
!!   GROUP OF CHARGES
       LOGICAL             :: groupcharge
       INTEGER             :: ngroup
       INTEGER,ALLOCATABLE :: group(:)
       REAL*8,ALLOCATABLE  :: qgr(:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! FFR-fock_qbias : Variable declaration
       logical :: apply_qbias
       real*8  :: weight
       real*8  :: current_time

       real*8,dimension(:),allocatable   :: Dvec
       real*8,dimension(:,:),allocatable :: Vmat
       real*8,dimension(:,:),allocatable :: sqsm
       real*8,dimension(:,:),allocatable :: Xmat,Xtrp,Ymat,Ytrp
       real*8,dimension(:,:),allocatable :: mirho_a,mirho_b,mirho
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       call g2g_timer_start('TD Open Shell')
       call g2g_timer_start('inicio')
       just_int3n = false
       ALLOCATE(factorial(NBCH))
!!------------------------------------!!
! Mulliken
       ipop=1
! Group of charges
       groupcharge=.false.
!!------------------------------------!!
#ifdef CUBLAS
       write(*,*) 'USING CUBLAS'
#endif
#ifdef TD_SIMPLE
        write(*,*) 'simple presition complex'
#else
        write(*,*) 'double presition complex'
#endif
       if(propagator.eq.2) then
          dt_magnus=tdstep
          dt_lpfrg=tdstep*0.10D0
          factorial(1)=1.0D0
#ifdef CUBLAS
       DO ii=1,NBCH
          factorial(ii)=1.0D0/ii
       ENDDO
#else     
       DO ii=2,NBCH
         factorial(ii)=factorial(ii-1)/ii
       ENDDO
#endif
       endif
       if(propagator.eq.1) then
          dt_lpfrg=tdstep
       endif
!!------------------------------------!!
!! FFR ADD:
       pert_steps=195
       lpfrg_steps=200
       chkpntF1a=185
       chkpntF1b=195
!--------------------------------------------------------------------!
! Pointers -
       Ndens=1
       E=0.0D0
       E1=0.0D0
       En=0.0D0
       E2=0.0D0
       idip=1
       ngeo=ngeo+1
       Im=(0.0D0,2.0D0)
       sq2=sqrt(2.D0)
       MM=M*(M+1)/2 
       MM2=M**2
       MMd=Md*(Md+1)/2
       Md2=2*Md
       M2=2*M
!
       ALLOCATE(xnano(M,M),xnano2(M,M),
     >   rhold_a(M,M),Y(M,M),
     >   rho1(M,M),rho_a(M,M),rho_b(M,M),rhonew_a(M,M),rhonew_b(M,M),
     >   rhold_b(M,M),fock_a(M,M),fock_b(M,M))
!
      if(propagator.eq.2) allocate (F1a_a(M,M),F1a_b(M,M),F1b_a(M,M),
     > F1b_b(M,M))   
!--------------------------------------------------------------------!
           call spunpack_rtc('L',M,rhoalpha,rho_a)
           call spunpack_rtc('L',M,rhobeta,rho_b)
!------------------------------------------------------------------------------!
c first i
            M1=1
c now Fold
            M3=M1+MM
c now S, F also uses the same position after S was used
            M5=M3+MM
c now G
            M7=M5+MM
c now Gm
            M9=M7+MMd
c now H
            M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
            M13=M11+MM
c aux ( vector for ESSl)
            M15=M13+M
c Least squares
            M17=M15+MM
c vectors of MO
            M18=M17+MMd
c weights (in case of using option )
            M19=M18+M*NCO
c RAM storage of two-electron integrals (if MEMO=T)
            M20 = M19 + natom*50*Nang   
c
            Nel=2*NCO+Nunp
c Initializations/Defaults
       write(*,*) ' TD CALCULATION  '
!--------------------------------------!
           niter=0
           D1=1.D0
           D2=1.D0
!--------------------------------------!
           Qc=0.0D0
           do i=1,natom
             Qc=Qc+Iz(i)
           enddo
           Qc=Qc-Nel
           Qc2=Qc**2
!------------------------------------------------------------------------------!
! Two electron integral with neighbor list.
!
            do i=1,natom
              natomc(i)=0
!
              do j=1,natom
                d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >            (r(i,3)-r(j,3))**2
                zij=atmin(i)+atmin(j)
                ti=atmin(i)/zij
                tj=atmin(j)/zij
                alf=atmin(i)*tj
                rexp=alf*d(i,j)
                if (rexp.lt.rmax) then
                  natomc(i)=natomc(i)+1
                  jatc(natomc(i),i)=j
                endif 
              enddo
            enddo
            do ii=nshell(0),1,-1
              nnps(nuc(ii))=ii
            enddo
            do ii=nshell(0)+nshell(1),nshell(0)+1,-1
              nnpp(nuc(ii))=ii
            enddo
            do ii=M,nshell(0)+nshell(1)+1,-1
              nnpd(nuc(ii))=ii
            enddo
!------------------------------------------------------------------------------!
! H H core, 1 electron matrix elements
            call int1(En)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! FFR-fock_qbias : Setup
       allocate(sqsm(M,M))
       allocate(Vmat(M,M),Dvec(M))
       allocate(Xmat(M,M),Xtrp(M,M),Ymat(M,M),Ytrp(M,M))
       allocate(mirho_a(M,M),mirho_b(M,M),mirho(M,M))

       call spunpack('L',M,RMM(M5),Smat)
       call sdiag_canonical(Smat,Dvec,Vmat,Xmat,Xtrp,Ymat,Ytrp)
       sqsm=matmul(Vmat,Ytrp)
       apply_qbias=.true.
       if (apply_qbias) then
         weight=0.1d0
         call fqbias_setup(natom,M,nuc,'atomgroup',weight)
       endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!--------------------------------------!
! SOLVENT CASE                        !OJO QUE ACA ANTES NO HABIA UN IF
            if(nsol.gt.0) then          
              call g2g_timer_start('intsol')
              call intsol(E1s,Ens,.true.)
              call g2g_timer_stop('intsol')
            endif
            E1=0.D0
            do k=1,MM
              E1=E1+RMM(k)*RMM(M11+k-1)
            enddo
            if(ipop.eq.1)then
              allocate(overlap(M,M),rhoscratch(M,M))
              call spunpack('L',M,RMM(M5),overlap)
            endif
!--------------------------------------!
c Diagonalization of S matrix, after this is not needed anymore
c s is in RMM(M13,M13+1,M13+2,...,M13+MM)
!--------------------------------------!
! ESSL OPTION
#ifdef essl
       call DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2)              
#endif
!--------------------------------------!
! LAPACK OPTION
#ifdef pack
       call dspev('V','L',M,RMM(M5),RMM(M13),X,M,RMM(M15),info)      
#endif
!--------------------------------------!
! Here, we obtain the transformation matrices X and Y for converting 
! from the atomic orbital to a molecular orbital basis (truncated
! during linear dependency elimination). 
! S is the overlap matrix
! s is the diagonal eigenvalue matrix of S
! U is the eigenvector matrix of S
! X=U s^(-1/2)
! matrix X's dimension is M*3M. In the first M*M terms it contains
! the transformation matrices and in the other M*2M terms it contains
! auxiliar matrices.
            call g2g_timer_start('inicio1')
            do j=1,M
              if (RMM(M13+j-1).lt.1.0D-06) then
                write(*,*) 'LINEAR DEPENDENCY DETECTED'
                do i=1,M
                  X(i,j)=0.0D0
                  Y(i,j)=0.0D0
                enddo
              else
                do i=1,M
                  Y(i,j)=X(i,j)*sqrt(RMM(M13+j-1))
                  X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))       ! hay que poner a las matrices de transformacion en el garcha_mod !!!!!!
                enddo
              endif
            enddo
!! CUBLAS ---------------------------------------------------------------------!
#ifdef CUBLAS
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=cmplx(X(i,j),0.0D0)
               ENDDO
            ENDDO
            stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
            stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
            stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
            if (stat.NE.0) then
            write(*,*) "X and/or Y memory allocation failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_real,X,M,devPtrX,M)
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrXc,M)
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=cmplx(Y(i,j),0.0D0)
               ENDDO
            ENDDO
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrY,M)
            if (stat.NE.0) then
            write(*,*) "X and/or Y setting failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
            rho1=0
#endif
!------------------------------------------------------------------------------!
! External Electric Field components
!
!       write(*,*) 'fx =', fx
!       write(*,*) 'fy =', fy
!       write(*,*) 'fz =', fz
!------------------------------------------------------------------------------!
! Rho is transformed to the orthonormal basis
#ifdef CUBLAS
       call g2g_timer_start('complex_rho_ao_to_on-cu')
       rho_a=basechange_cublas(M,rho_a,devPtrY,'dir')
       rho_b=basechange_cublas(M,rho_b,devPtrY,'dir')
       call g2g_timer_stop('complex_rho_ao_to_on-cu')
#else
       rho_a=matmul(ytrans,rho_a)
       rho_a=matmul(rho_a,y)
       rho_b=matmul(ytrans,rho_b)
       rho_b=matmul(rho_b,y)
#endif
            call g2g_timer_start('int2')
            call int2()
            call g2g_timer_stop('int2')
            call g2g_timer_start('int3mmem')
            if (memo) then
              deallocate (kkind,kkinds)
              deallocate(cool,cools)
            endif
            call int3mem()
            call int3mems()
            call g2g_timer_stop('int3mmem')
!------------------------------------------------------------------------------!
#ifdef CUBLAS
            call CUBLAS_FREE(devPtrY)
#endif
!            deallocate(y)
!------------------------------------------------------------------------------!
            call g2g_timer_stop('inicio')
!##############################################################################!
! HERE STARTS THE TIME EVOLUTION
!##############################################################################!
            write(*,*) 'PROPAGATION'
            do 999 istep=1, ntdstep
!--------------------------------------!
              call g2g_timer_start('iteration')
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >                             .and.(.not.tdrestart)) then
                 t=(istep-1)*tdstep*0.1
              else
                 t=19.9*tdstep
                 t=t+(istep-200)*tdstep
              endif

              if (propagator.eq.1) then
                 t=(istep-1)*tdstep
              endif

              t=t*0.02419
              write(*,*) 'evolution time (fs)  =', t
!--------------------------------------!
              call int3lu(E2)
              call g2g_solve_groups(0,Ex,0)
              write(*,*) '! step & energy', istep,E
              E1=0.0D0


c ELECTRIC FIELD CASE - Type=gaussian (ON)
            if(istep.lt.pert_steps) then
               if (field) then
                 call dip(ux,uy,uz)
                 if (exter) then
                   g=1.0D0
                   fac=2.54D0
                   fxx=fx*exp(-0.2*(real(istep-50))**2)
                   fyy=fy*exp(-0.2*(real(istep-50))**2)
                   fzz=fz*exp(-0.2*(real(istep-50))**2)
                   write(*,*) fxx,fyy,fzz
!
                 else
                   g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
                   Fx=ux/2.54D0
                   Fy=uy/2.54D0
                   Fz=uz/2.54D0
                   fac=(2.54D0*2.00D0)
!
                 endif
                 call intfld(g,Fxx,Fyy,Fzz)
                 E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/fac -
     >        0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
              endif
            endif
!------------------------------------------------------------------------------!
! E1 includes solvent 1 electron contributions
            do k=1,MM
              E1=E1+RMM(k)*RMM(M11+k-1)
            enddo



!------------------------------------------------------------------------------!
! Here we obtain the fock matrix in the molecular orbital (MO) basis.
! where U matrix with eigenvectors of S , and s is vector with
! eigenvalues
           call g2g_timer_start('fock')
            call spunpack('L',M,RMM(M5),fock_a)
            call spunpack('L',M,RMM(M3),fock_b)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! FFR-fock_qbias : 
         if (apply_qbias) then
           current_time=t
           call fqbias_calc(M,current_time,sqsm,fock_a)
           call fqbias_calc(M,current_time,sqsm,fock_b)
         endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

#ifdef CUBLAS
            fock_a = basechange_cublas(M,fock_a,devPtrX,'dir')
            fock_b = basechange_cublas(M,fock_b,devPtrX,'dir')
#else
            fock_a=matmul(xtrans,fock_a)
            fock_a=matmul(fock_a,xmm)
            fock_b=matmul(xtrans,fock_b)
            fock_b=matmul(fock_b,xmm)
#endif
            call g2g_timer_stop('fock')
c Fock triangular matrix contained in RMM(M5,M5+1,M5+2,...,M5+MM) is copied to square matrix fock.
             call sprepack('L',M,RMM(M5),fock_a)
             call sprepack('L',M,RMM(M3),fock_b)
c Now fock is stored in molecular orbital basis.
c
!  stores F1a and F1b for magnus propagation
            if((propagator.eq.2) .and. (.not.tdrestart)) then
               if(istep.eq.chkpntF1a) then
                  F1a_a=fock_a
                  F1a_b=fock_b         
               endif
               if(istep.eq.chkpntF1b) then
                  F1b_a=fock_a
                  F1b_b=fock_b         
               endif         
            endif
!  stores F1a and F1b checkpoints to restart the dynamics      !ESTO HAY QUE ADAPTARLO A OPEN-SHELL!!!!!!!!!!
!            if(writedens .and. propagator.eq.2) then
!               kk=istep+5
!               ii=istep+15
!            if(mod (kk,500) == 0) then
!               open(unit=7624,file='F1b.restart')
!               rewind 7624
!               do i=1,M
!                  do j=1,M
!                     write(7624,*) fock(i,j)
!                  enddo
!               enddo
!               endif 
!               if(mod (ii,500) == 0) then
!                 open(unit=7625,file='F1a.restart')
!                 rewind 7625
!                 do i=1,M
!                    do j=1,M
!                       write(7625,*) fock(i,j)
!                    enddo
!                 enddo
!               endif
!            endif
            E=E1+E2+En
            if (sol) then
                E=E+Es
            endif
!--------------------------------------------------------------------!
            if ((propagator.eq.1).or.
     >      (((propagator.eq.2).and.(istep.lt.lpfrg_steps))
     >      .and. (.not.tdrestart))) then
           write(*,*) 'Verlet'
c In the first step of the propagation we extrapolate rho back in time
c using Verlet algorithm to calculate rhold.
c using matmul 
c           if(istep.eq.1) then
c             rhold=rho+(tdstep*Im*(matmul(fock,rho)))
c             rhold=rhold-(tdstep*Im*(matmul(rho,fock)))
c           endif
c using conmutc
              if(istep.eq.1) then

#ifdef CUBLAS
               call g2g_timer_start('cuconmut')
               rhold_a=commutator_cublas(fock_a,rho_a)
               rhold_a=rho_a+dt_lpfrg*(Im*rhold_a)
               rhold_b=commutator_cublas(fock_b,rho_b)
               rhold_b=rho_b+dt_lpfrg*(Im*rhold_b)
               call g2g_timer_stop('cuconmut')
#else
               rhold_a=commutator(fock_a,rho_a)
               rhold_a=rho_a+dt_lpfrg*(Im*rhold_a)
               rhold_b=commutator(fock_b,rho_b)
               rhold_b=rho_b+dt_lpfrg*(Im*rhold_b)
#endif

              endif
!####################################################################!
! DENSITY MATRIX PROPAGATION USING VERLET ALGORITHM
! using matmul:
c           rhonew=rhold-(tdstep*Im*(matmul(fock,rho)))
c           rhonew=rhonew+(tdstep*Im*(matmul(rho,fock)))
c--------------------------------------c
! using conmutc:
                call g2g_timer_start('Verlet')
#ifdef CUBLAS
               call g2g_timer_start('cuconmut')
               rhonew_a=commutator_cublas(fock_a,rho_a)
               rhonew_a=rhold_a-dt_lpfrg*(Im*rhonew_a)
               rhonew_b=commutator_cublas(fock_b,rho_b)
               rhonew_b=rhold_b-dt_lpfrg*(Im*rhonew_b)
               call g2g_timer_stop('cuconmut')
#else
               rhonew_a=commutator(fock_a,rho_a)
               rhonew_a=rhold_a-dt_lpfrg*(Im*rhonew_a)
               rhonew_b=commutator(fock_b,rho_b)
               rhonew_b=rhold_b-dt_lpfrg*(Im*rhonew_b)
#endif
               call g2g_timer_stop('Verlet')
c Density update (rhold-->rho, rho-->rhonew)
               rhold_a=rho_a
               rhold_b=rho_b
               rho_a=rhonew_a
               rho_b=rhonew_b
! END OF VERLET PROPAGATOR
!####################################################################!
              else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! DENSITY MATRIX PROPAGATION USING MAGNUS ALGORITHM
                 write(*,*) 'Magnus'
#ifdef CUBLAS
                call g2g_timer_start('cupredictor')
                call cupredictor_op(F1a_a,F1b_a,F1a_b,F1b_b,fock_a,
     >               fock_b,rho_a,rho_b,factorial,devPtrX,Fxx,
     >               Fyy,Fzz,g,devPtrXc)
                call g2g_timer_stop('cupredictor')
                call g2g_timer_start('cumagnus')
                call cumagnusfac(fock_a,rho_a,rhonew_a,M,NBCH,dt_magnus,
     >          factorial)
                call cumagnusfac(fock_b,rho_b,rhonew_b,M,NBCH,dt_magnus,
     >          factorial)
                call g2g_timer_stop('cumagnus')
#else
                call g2g_timer_start('predictor')
                call predictor_op(F1a_a,F1b_a,F1a_b,F1b_b,fock_a,fock_b,
     >          rho_a,rho_b,factorial,Fxx,Fyy,Fzz,g)
                call g2g_timer_stop('predictor')
                call g2g_timer_start('magnus')
                call magnus(fock_a,rho_a,rhonew_a,M,NBCH,dt_magnus,
     >          factorial)
                call magnus(fock_b,rho_b,rhonew_b,M,NBCH,dt_magnus,
     >          factorial)
                call g2g_timer_stop('magnus')
#endif
                 F1a_a=F1b_a
                 F1a_b=F1b_b
                 F1b_a=fock_a
                 F1b_b=fock_b
                 rho_a=rhonew_a
                 rho_b=rhonew_b
! END OF MAGNUS PROPAGATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
              endif
!####################################################################!
c Here we transform the density to the atomic orbital basis and take the real part of it. The imaginary part of the density 
c can be descarted since for a basis set of purely real functions the fock matrix is real and symetric and depends only on 
c the real part of the complex density matrix. (This won't be true in the case of hybrid functionals)
c with matmul:
#ifdef CUBLAS

             call g2g_timer_start('complex_rho_on_to_ao-cu')   !TODO: VER COMO HACER ESTO CON UNA SOLA MATRIZ AUXILIAR COPIANDO A LA POSICION CORRESPONDIENTE DEL RMM
             rho1=basechange_cublas(M,rho_a,devPtrXc,'inv')
             mirho_a=rho1 ! FFR-fock_qbias
             call sprepack_ctr('L',M,RMM,rho1)
             call sprepack_ctr('L',M,rhoalpha,rho1)
             rho1=basechange_cublas(M,rho_b,devPtrXc,'inv')
             mirho_b=rho1 ! FFR-fock_qbias
             call sprepack_ctr('L',M,rhobeta,rho1)
             DO i=1,MM
                RMM(i)=RMM(i)+rhobeta(i)
             ENDDO
             call g2g_timer_stop('complex_rho_on_to_ao-cu')

#else
             call g2g_timer_start('complex_rho_on_to_ao')
             rho1=matmul(x,rho_a)
             rho1=matmul(rho1,xtrans)
             mirho_a=rho1 ! FFR-fock_qbias
             call sprepack_ctr('L',M,RMM,rho1)
             call sprepack_ctr('L',M,rhoalpha,rho1)
             rho1=matmul(x,rho_b)
             rho1=matmul(rho1,xtrans)
             mirho_b=rho1 ! FFR-fock_qbias
             call sprepack_ctr('L',M,rhobeta,rho1)
             DO i=1,MM
                RMM(i)=RMM(i)+rhobeta(i)
             ENDDO
             call g2g_timer_stop('complex_rho_on_to_ao')
#endif
!       rho1=REAL(rho1)
c with matmulnanoc:
c          call matmulnanoc(rho,xtrans,rho1,M)
c          rho1 = REAL(rho1)
! Stores the density matrix each 500 steps to be able to restart the dynamics
!              if(writedens) then                                        !TODO: VAMOS A TENER QUE RESCRIBIR TODO ESTO
!                 if(mod (istep,500) == 0) then
!                     open(unit=5374,file='rho.restart')
!                     rewind 5374
!                     do j=1,M
!                        do k=1,M
!                           write(5374,*) rho1(j,k)
!                        enddo
!                     enddo
!                  endif
! In the last step density matrix is stored
!                  if (istep.eq.ntdstep) then
!                    open(unit=44,file='rholast')
!                    do j=1,M
!                       do k=1,M
!                          write(44,*) rho1(j,k)
!                       enddo
!                    enddo
!                  endif
!              endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! FFR-fock_qbias : Write Lowdin Group Population
       call fqbias_lowdinpop(M,natom,mirho_a,sqsm,Iz,'lowdina_gr.txt')
       call fqbias_lowdinpop(M,natom,mirho_b,sqsm,Iz,'lowdinb_gr.txt')
       mirho=mirho_a+mirho_b
       call fqbias_lowdinpop(M,natom,mirho,sqsm,Iz,'lowdin_gr.txt')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!###################################################################!
!# DIPOLE MOMENT CALCULATION
              if(istep.eq.1) then
                open(unit=134,file='x.dip')
                open(unit=135,file='y.dip')
                open(unit=136,file='z.dip')
        write(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
        write(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
        write(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
              endif
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >      .and. (.not.tdrestart)) then
                  if(mod ((istep-1),10) == 0) then
                     call g2g_timer_start('DIPOLE') 
                     call dip(ux,uy,uz)
                     call g2g_timer_stop('DIPOLE')
                     write(134,901) t,ux
                     write(135,901) t,uy
                     write(136,901) t,uz
                  endif
              else
                  call g2g_timer_start('DIPOLE')
                  call dip(ux,uy,uz)
                  call g2g_timer_stop('DIPOLE')
                  write(134,901) t,ux
                  write(135,901) t,uy
                  write(136,901) t,uz
              endif
c u in Debyes
!# END OF DIPOLE MOMENT CALCULATION
c------------------------------------------------------------------------------------
!       if (nopt.ne.3) then
!       write(*,300) niter,DAMP,E
!       endif
c-------------------------MULLIKEN CHARGES-----------------------------------------------!
                if(istep.eq.1) then
                  open(unit=1111111,file='Mulliken')
                  if (groupcharge) then
                      open(unit=678,file='Mullikin')
                      allocate(group(natom))
                      ngroup=0
                      do n=1,natom
                        read(678,*) kk
                        group(n)=kk
                        if (kk.gt.ngroup) ngroup=kk
                      enddo
                      allocate(qgr(ngroup))
                      close(unit=678)
                      open(unit=678,file='MullikenGroup')
                  endif
                endif
!
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >      .and. (.not.tdrestart)) then
                if(mod ((istep-1),10) == 0) then
                   call spunpack_rho('L',M,RMM,rhoscratch)
                   rhoscratch=matmul(overlap,rhoscratch)
                   do n=1,natom
                      q(n)=Iz(n)
                   enddo
                   do i=1,M
                      q(Nuc(i))=q(Nuc(i))-rhoscratch(i,i)
                   enddo
                   if(groupcharge) qgr=0.0d0
                   do n=1,natom
                      write(1111111,760) n,Iz(n),q(n)
                      if(groupcharge) then
                        qgr(group(n))=qgr(group(n))+q(n)
                      endif
                   enddo
                   if(groupcharge) then
                    do n=1,ngroup
                       write(678,761) n,n,qgr(n)
                    enddo
                    write(678,*) '------------------------------------'
                   endif
                 endif
              else
                 call spunpack_rho('L',M,RMM,rhoscratch)
                 rhoscratch=matmul(overlap,rhoscratch)
                 do n=1,natom
                    q(n)=Iz(n)
                 enddo
                 do i=1,M
                    q(Nuc(i))=q(Nuc(i))-rhoscratch(i,i)
                 enddo
                 if(groupcharge) qgr=0.0d0
                 do n=1,natom
                    write(1111111,760) n,Iz(n),q(n)
                    if(groupcharge) then
                       qgr(group(n))=qgr(group(n))+q(n)
                    endif
                 enddo
                 if(groupcharge) then
                   do n=1,ngroup
                       write(678,761) n,n,qgr(n)
                   enddo
                 endif
                 write(678,*) '------------------------------------'
              endif
!-----------------------------------------------------------------------------------!

c      write(*,*) 'Coulomb E',E2-Ex,Ex
               call g2g_timer_stop('iteration')
               write(*,*)
 999           continue
!
!##############################################################################!
! HERE FINISHES THE PROPAGATION
!##############################################################################!

 995   continue
c
c
#ifdef CUBLAS
         call CUBLAS_FREE ( devPtrX )
#endif
         if (memo) then
            deallocate (kkind,kkinds)
            deallocate(cool,cools)
         endif
         if(propagator.eq.2) then
           deallocate (F1a_a,F1b_a,F1a_b,F1b_b)
         endif
         if (GRAD) then
            if(nopt.eq.0) then
              write(*,*)
              write(*,600)
              write(*,610)
              write(*,620) E1,E2-Ex,En
              if (sol) then
                 write(*,615)
                 write(*,625) Es
              endif
            endif
            write(*,*)
            write(*,450) E
         else
            E=E-Ex
         endif
c calculation of energy weighted density matrix
c
          kk=0
          do 307 j=1,M
             do 307 i=j,M
                kk=kk+1
                RMM(M15+kk-1)=0.D0
                if(i.eq.j) then
                    ff=2.D0
                else
                    ff=4.D0
                endif
                do 309 k=1,NCO
                   RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1)
     >  *ff*X(i,M2+k)*X(j,M2+k)
 309  continue
 307   continue

c
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
c        if (icharge.eq.1) then
c          Q1=-(2*NCO+Nunp)
c         do n=1,natom
c          Q1=Q1+Iz(n)
c         enddo
c         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,RMM,map,Q1)
c        endif
c
c outputs final  MO ---------------------
!      do l=1,M
c      do n=1,NCO+3
!      do n=1,M
!        X(indexii(l),M+n)=X(l,M2+n)
!      enddo
!      enddo
!
!      do l=1,M
!        write(2,400) (X(l,M+n),n=1,NCO)
!      enddo
!--------------------------------------!
! Writes down MO coefficients and orbital energies
!       write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
!       do n=1,NCO
!         write(29,850) n,RMM(M13+n-1)
!         write(29,400) (X(l,M+n),l=1,M)
!       enddo
!       do n=NCO+1,M
!         write(29,851) n,RMM(M13+n-1)
!         write(29,400) (X(l,M+n),l=1,M)
!       enddo
!       close(29)
!--------------------------------------!
c
c
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
       call g2g_timer_stop('TD Open Shell')
#ifdef CUBLAS
       deallocate(xnano,rho_a,fock_a,rho_b,fock_b,rho1)
#else
       deallocate(xnano,rho_a,fock_a,rho_b,fock_b,rho1)
#endif
       DEALLOCATE(factorial)
       if(ipop.eq.1) deallocate(rhoscratch,overlap)
!------------------------------------------------------------------------------!
 500  format('SCF TIME ',I6,' sec')
 450  format ('FINAL ENERGY = ',F19.12)
 400  format(4(E14.7E2,2x))
 300  format(I3,E14.6,2x,F14.7)
 600  format('  ENERGY CONTRIBUTIONS IN A.U.')
 610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  format(2x,'SOLVENT')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 625  format(F14.7)
 760  format(I3,9x,I3,6x,F10.4)
 761  format(I3,9x,I3,6x,F14.8)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
 850  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
 851  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7,
     >    '(NON OCC.)')
 900  format(F15.9,2x,3(F15.9,2x),2x,F15.9)
 901  format(F15.9,2x,F15.9)
 777  format(4(F8.4,2x))
 776  format (3(F8.4,2x))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
