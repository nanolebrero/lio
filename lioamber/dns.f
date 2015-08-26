c function calculating density functionals
c standard, for local density functionals, recalculates
c everything on every iteration
c
      SUBROUTINE DNS(DENSX,Xi)    
      use garcha_mod
c
      implicit real*8 (a-h,o-z)
      dimension Xi(3)
      real*8, dimension(:), ALLOCATABLE :: ds,F,W
c
c now we should evaluate all same loops as the ones used for
c 1 electron matrix elements, but doing only products
c then, the particular density functional wanted is calculated
c
       allocate(ds(ntq),F(M),W(ng))

        do 221 k=1,natom
         ds(k)=(xi(1)-r(k,1))**2+(xi(2)-r(k,2))**2+(xi(3)-r(k,3))**2
 221     continue

      fc=1.D0
      if (NORM) then
       fc=1.D0/sqrt(3.D0)
      endif
c
      DENSX=0.D0
      DENSX1=0.0D0
c
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
c basis functions evaluated at r are calculated
c
      do 1 i=1,M
        W(i)=0.D0
 1      F(i)=0.D0
c
c --- s  case -------
      do 10 i=1,ns
c
      di=ds(Nuc(i))
c
      do 15 ni=1,ncont(i)
c
      rexp=a(i,ni)*di
c
      if (rexp.gt.30.D0) go to 16
      t=exp(-rexp)
      F(i)=F(i)+t*c(i,ni)
  16  continue
  15  continue
c
  10  continue
c
c--- p  case -------------
      do 20 i=ns+1,ns+np,3
c
      di=ds(Nuc(i))
c
      do 20 ni=1,ncont(i)

      rexp=a(i,ni)*di
      if (rexp.gt.30.D0) goto 21
      t=exp(-rexp)*c(i,ni)
      do 25 l1=1,3
c
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i+l1-1
c
      F(ii)=F(ii)+term
  25  continue
c
  21  continue
  20  continue
c
c-- d case  ------------
      do 40 i=ns+np+1,M,6
c
      di=ds(Nuc(i))
c
      do 40 ni=1,ncont(i)
c
      rexp=a(i,ni)*di
c
      if (rexp.gt.30.) goto 41
      t=exp(-rexp)*c(i,ni)
      do 45 l1=1,3
c
      t1=xi(l1)-r(Nuc(i),l1)
      do 45 l2=1,l1
      t2=xi(l2)-r(Nuc(i),l2)
      if (l1.eq.l2) then
       t2=t2*fc
      endif
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
  45  continue
c
  41  continue
  40  continue
c
c now calculation of vector W : density matrix scalar F
c
c 
      k=0
      do 50 j=1,M
c
       if (F(j).eq.0.0D0) then
        k=k+M-j+1
        goto 50
       endif
c
      do 51 i=j,M
       k=k+1
 51   W(j)=W(j)+RMM(k)*F(i)
 50   continue
c
      do 60 i=1,M
       DENSX=DENSX+F(i)*W(i)
  60  continue
c
       deallocate(ds,F,W)
      return

      end
