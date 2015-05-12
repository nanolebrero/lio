	subroutine intECP
	use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP,nECP,bECP, aECP,Zcore, Lmax, expnumbersECP
	implicit none

	integer z,l,t

!Escribe coeficientes como testeo
        do z=1,118
                do l=0, Lmax(Z)
                        do t=1, expnumbersECP(Z,l)
                                write(*,9018) Z,l,t,nECP(Z,l,t), bECP(Z,l,t), aECP(Z,l,t)
                        end do
                end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!terminos <Ai|V|Aj>
!termino local























        9018 format(/1x,'Z =',i4,2x, 'L =',i4,2x,'coefnumber =',i3,2x,  'n =',i2,2x, 'b =', f15.5,2x,'c =',f15.5)



	contains

	DOUBLE PRECISION FUNCTION VijlocalAAA
	implicit none
	double precision :: acumk, acuml, acumm
	integer :: k,l,m
	VijlocalAAA=0.d0
!	do 
	end function VijlocalAAA

!subrutinas mixtas
	DOUBLE PRECISION FUNCTION SAAA (k,alpha,Ccoef,na,ma,la,nb,mb,lb,necp)
	use ECP_mod, only :  Qnl
	implicit none
	integer, intent(in) :: na,ma,la,nb,mb,lb,necp
	double precision, intent (in) :: k,alpha,Ccoef
	integer :: lambda, lmax
	Qnl=0.d0
	SAAA=0.d0
	lmax=na+ma+la+nb+mb+lb
!lmax  = 0 para <s||s>, 1 para <s||p>, 2 para <s||d>, ... , 6 para <d||d>
!Ccoef = exponente basei + exponente base j + exponente ecp
	call Qtype1(K,Ccoef,lmax,necp)
		do lambda=0, lmax
!			SAAA=SAAA + Qnl(na+ma+la+nb+mb+lb+necp,lambda)*OMEGA1(K,lambda,na+nb,ma+mb,la+lb)
		end do
	return
	end function SAAA

!subrutinas angulares

	DOUBLE PRECISION FUNCTION OMEGA2(K,lambda,l,m,a,b,c)
!Esta funcion devuelve el valor de omega2 evaluado en el vector K
!OMEGA2(K,lambda,l,m,a,b,c) = Suma(o=-lambda a o=lambda) ( Y lambda o(K)*  int (x/r)^a * (y/r)^b * (z/r)^c Ylambda o() Ylm()  d angulo solido )

!no esta testeada a full ya q es muy engorroso
	use ECP_mod, only : intangular,angularint
        implicit none
        DOUBLE PRECISION, intent(in), dimension(3) :: K
        integer, intent(in) :: lambda,l,m,a,b,c
        DOUBLE PRECISION,dimension(3) :: Kun
        integer :: o,r,s,t,u,v,w
        Double precision :: SUM1, SUM2
	Kun=K/sqrt(K(1)**2 + K(2)**2 + K(3)**2)
	SUM1=0.d0
        SUM2=0.d0
	OMEGA2=0.d0
!	write(15,*) "lkkk"
!	write(15,*) "lambda,l,m,a,b,c",lambda,l,m,a,b,c
	do o=-lambda,lambda
		do r=0,lambda
			do s=0,lambda-r
				t=lambda-r-s
                                SUM1=SUM1+Ucoef(lambda,o,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
!				write(15,*) "sum1",sum1
				do u=0,l
					do v = 0,l-u
						w=l-u-v
			                        SUM2=SUM2+Ucoef(lambda,o,r,s,t)*Ucoef(l,m,u,v,w)*angularint(a+r+u,b+s+v,c+t+w)
!						write(15,*) "a+r+u,b+s+v,c+t+w",a+r+u,b+s+v,c+t+w
!						write(15,*) "Ucoef(lambda,o,r,s,t)*Ucoef(l,m,u,v,w)*angularint(a+r+u,b+s+v,c+t+w)",Ucoef(lambda,o,r,s,t),Ucoef(l,m,u,v,w),angularint(a+r+u,b+s+v,c+t+w)
!						write(15,*) "sum2",sum2,"u,v,w",u,v,w
					end do
				end do
			end do
		end do
!		write(15,*) "sum1",sum1,"sum2",sum2
		OMEGA2=OMEGA2+SUM1*SUM2
!		write(15,*) o,omega2
                SUM1=0.d0
                SUM2=0.d0
	end do
	end function OMEGA2

	DOUBLE PRECISION FUNCTION OMEGA1(K,l,a,b,c)
!Esta funcion devuelve el valor de omega evaluado en el vector K
!OMEGA1(K,l,a,b,c) = Suma(u=-l a u=l) ( Ylu(K)* int(x/r)^a * (y/r)^b * (z/r)^c Ylu() d angulo solido )
	use ECP_mod, only : intangular,angularint
	implicit none
	DOUBLE PRECISION, intent(in), dimension(3) :: K
	integer, intent(in) :: l,a,b,c 
	DOUBLE PRECISION,dimension(3) :: Kun
	integer :: r,s,t,u
	Double precision :: SUM1, SUM2
	SUM1=0.d0
	SUM2=0.d0
	OMEGA1=0.d0
	Kun=K/sqrt(K(1)**2 + K(2)**2 + K(3)**2)

!	write(*,*) "l",l
	do u=-l,l
!		write(*,*) "u",u
		do r=0,l
!			write(*,*) "r",r
			do s=0,l-r
!				write(*,*) "s",s
				t=l-r-s
!				write(*,*) "t",t
				SUM1=SUM1+Ucoef(l,u,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
				SUM2=SUM2+Ucoef(l,u,r,s,t)*angularint(a+r,b+s,c+t)
!Usa la funcion angular int en lugar del array int angular. Hay que agrandar el array
!				write(*,*) SUM1, SUM2, angularint(a+r,b+s,c+t),a+r,b+s,c+t
			end do
		end do
!		write(*,*) SUM1, SUM2, intangular(a+r,b+s,c+t),a+r,b+s,c+t
		OMEGA1=OMEGA1+SUM1*SUM2
		SUM1=0.d0
		SUM2=0.d0
	end do
	return
	end function OMEGA1


	DOUBLE PRECISION FUNCTION Ucoef(l,m,lx,ly,lz)
	use ECP_mod, only :  l0,l1,l2,l3
	implicit none
	integer, intent(in) :: l,m,lx,ly,lz
	if (l .eq. 0) then
	Ucoef=l0(1)
	elseif (l .eq. 1) then
	Ucoef=l1(2*lx+ly+1,m)
	elseif (l .eq. 2) then
	Ucoef=l2(0.5*lx*(7-lx)+ly+1,m)
        elseif (l .eq. 3) then
	Ucoef=l3(0.5*lx*(9-lx)+ly+1,m)
	else
	write(*,*) "ECP error l is grater than 3"
	end if
	return
	end function Ucoef



        end subroutine intECP


!subrutinas Radiales

        subroutine Qtype1(K,Ccoef,lmax,necp)
!agrega a la matriz Qnl los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl ya que hay q llamar a esta rutina por cada termino del pseudopotencial
!luego hay que optimizar la rutina para que calcule Qnl a partir de otros Qnl
	use ECP_mod, only :  alpha, betha, Bn, Cn, Qnl
	implicit none
	double precision, intent(in) :: K,Ccoef
!lmax  = 0 para <s||s>, 1 para <s||p>, 2 para <s||d>, ... , 6 para <d||d>	
!necp corresponde al exponente del paseudopotencial en r^n
	integer, intent(in) :: necp,lmax
!nmin y nmax dan el rango en que tiene q calcular Q
	integer :: nmin,nmax
!variables auxiliares
	integer :: n,l,i
        double precision :: acoef, gam, acum
	
!Caso 1-calcula todas las integrales
	nmin=necp
	nmax=necp+lmax

	acoef=K/(2.d0*Ccoef)
	gam=0.5d0*exp(K**2/(4*Ccoef))

	call ByC(acoef,ccoef,nmin,nmax,Bn,Cn)


	do n=nmin,nmax
	do l=0,lmax

	acum=0.d0
	do i=l,1,-2
		acum=acum+alpha(l,i)*Bn(n-i)/k**i
	end do
	do i=l+1,1,-2
		acum=acum+betha(l+1,i)*Cn(n-i)/k**i
	end do
		Qnl(n,l)=Qnl(n,l)+acum*gam
		acum=0.d0
	end do
	end do

!escribo todo como test
!	write(*,*) "entro en QNL"
!	do n=nmin,nmax
!        do l=0,lmax
!	write(*,*) n,l,Qnl(n,l)
!	end do
!	end do
        End subroutine Qtype1


        subroutine ByC(acoef,ccoef,nmin,nmax,Barray,Carray)
!calcula los coeficientes B y C 
	use ECP_mod, only : DAW,DAWERF,NEXTCOEF,pi,pi12
	IMPLICIT NONE
!acoef,ccoef son los coeficientes para el calculo de B y C
! A(ccoef,acoef)= int exp(-ccoef*x^2)(x+acoef)^n dx from -acoef to inf
!nmin,nmax delimitan el rango de n
!Bn=An(ccoef,acoef) + An(ccoef,-acoef)
!Cn=An(ccoef,acoef) - An(ccoef,-acoef)
	DOUBLE PRECISION, intent(in) :: acoef,ccoef
	integer, intent(in) :: nmin,nmax
	DOUBLE PRECISION, intent(inout), dimension (-12:14) :: Barray,Carray
!variables auxiliares
	DOUBLE PRECISION :: c0sq,ncos,ca
	integer :: i

	c0sq=sqrt(ccoef)
	Barray(0)=pi12/c0sq
	Carray(0)=Barray(0)*erf(acoef*c0sq)

	if (nmax>0) then
	        Barray(1)= exp(-ccoef*acoef**2)/ccoef + acoef*Carray(0)
	        Carray(1)= Barray(0)*acoef
		do i=2,nmax
			ncos=(i-1)/(2*ccoef)
			Barray(i)=ncos*Barray(i-2)+acoef*Carray(i-1)
			Carray(i)=ncos*Carray(i-2)+acoef*Barray(i-1)
		end do
	end if

        if (nmin<0) then
		Barray(-1)=2*pi12*DAWERF(acoef*c0sq)
		Carray(-1)=2*pi12*DAW(acoef*c0sq)
		if (nmin<-1) then
			ca=2*ccoef*acoef
			Barray(-2)=ca*Carray(-1)-2*ccoef*Barray(0)
			Carray(-2)=2*ca*exp(-ccoef*acoef**2)+ca*Barray(-1)-2*ccoef*Carray(0)
			do i=-3,nmin,-1
				Barray(i)=NEXTCOEF(1,i,ca,exp(-ccoef*acoef**2),ccoef,Carray(i+1), Barray(i+2))
				Carray(i)=NEXTCOEF(-1,i,ca,exp(-ccoef*acoef**2),ccoef,Barray(i+1), Carray(i+2))
			end do
		end if
	end if
        end subroutine ByC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!this routine obtain Q(l1,l2,n,k1,k2,a) where
!Q(l1,l2,n,k1,k2,a)= int r^n exp(-ar^2) Ml1(k1r) Ml2(k2r) dr form r=0 to inf
!where Ml are the modified spherical Bessel function of the first kind.

        subroutine Qtype2(Ka,Kb,Ccoef,l1max,l2max,necp)
!agrega a la matriz Qnl1l2 los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl1l2 ya que hay q llamar a esta rutina por cada termino del pseudopotencial
!luego hay que optimizar la rutina para que calcule Qnl1l2 a partir de otros Qnl1l2
	use ECP_mod, only :  alpha, betha, rho, tau, sigma, sigmaR, Qnl1l2
	implicit none
	double precision, intent(in) :: Ka,Kb,Ccoef
!l1max y l2max = 0 para s, 1 para p, 2 para d, etc	
!n corresponde al exponente del paseudopotencial en r^n
	integer, intent(in) :: necp,l1max, l2max
	integer :: nmin,nmax
!variables auxiliares
	integer :: i,j,n,l1,l2
	double precision :: alfok, betok, acum1
	
!Caso 1-calcula todas las integrales
	nmin=necp
	nmax=necp+l1max+l2max
!	write(*,*) nmin,nmax
	call integrals(Ka,Kb,Ccoef,necp-2-l1max-l2max,necp+l1max+l2max-2)
	acum1=0.d0

	do n=nmin,nmax
	do l1=0,l1max
	do l2=0,l2max
	do i=l1,1,-2
		alfok=alpha(l1,i)/Ka**i
		do j=l2,1,-2
			if (tau(n-i-j) == 0.d0) write(*,*) "Error, no calculo tau",n-i-j
			acum1=acum1+alpha(l2,j)*tau(n-i-j)/Kb**j
		end do
		do j=l2+1,1,-2
			acum1=acum1+betha(l2+1,j)*sigmaR(n-i-j)/Kb**j
			if (sigmaR(n-i-j) == 0.d0) write(*,*) "Error, no calculo sigmaR",n-i-j
		end do
		Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*alfok
		acum1=0.d0
	end do
	do i=l1+1,1,-2
		betok=betha(l1+1,i)/Ka**i
		do j=l2,1,-2
			acum1=acum1+alpha(l2,j)*sigma(n-i-j)/Kb**j
			if (sigma(n-i-j) == 0.d0) write(*,*) "Error, no calculo sigma",n-i-j
		end do
		do j=l2+1,1,-2
			acum1=acum1+betha(l2+1,j)*rho(n-i-j)/Kb**j
			if (rho(n-i-j) == 0.d0) write(*,*) "Error, no calculo rho",n-i-j
		end do
		Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*betok
                acum1=0.d0
	end do

	end do
	end do
	end do

!escribo todo como chekeo
	if (1) then
	write(*,*) "   n    l1    l2    Qnl1l2"
        do n=nmin,nmax
        do l1=0,l1max
        do l2=0,l2max
	write(*,*) n,l1,l2,Qnl1l2(n,l1,l2)
	end do
        end do
        end do
	end if

        End subroutine Qtype2

	subroutine integrals(Ka,Kb,Ccoef,nmin,nmax)
!obtain integrals rho, tau, sigma and sigmaR from n between nmin and nmax
! rho(n) = int exp(-cr^2) * sinh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! sigma(n) = int exp(-cr^2) * sinh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf
! sigmaR(n) = int exp(-cr^2) * cosh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! tau(n) = int exp(-cr^2) * cosh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf

	use ECP_mod, only : Bn1,Bn2,Cn1,Cn2,rho, tau, sigma, sigmaR
	implicit none
	integer, intent(in) :: nmin,nmax
	double precision, intent(in) :: Ka,Kb,Ccoef
	double precision, dimension(2) :: acoef,gammacoef
	double precision :: signo
	acoef(1)=0.5d0*(Ka+Kb)/Ccoef
	acoef(2)=0.5d0*abs(Ka-Kb)/Ccoef
	gammacoef(1)=0.25d0*exp(Ccoef*acoef(1)**2)
	gammacoef(2)=0.25d0*exp(Ccoef*acoef(2)**2)
!	call ByCdoble(acoef(1),acoef(2),Ccoef,nmin,nmax)
	call ByC(acoef(1),Ccoef,nmin,nmax,Bn1,Cn1)
	call ByC(acoef(2),Ccoef,nmin,nmax,Bn2,Cn2)
	rho=gammacoef(1)*Bn1-gammacoef(2)*Bn2
	tau=gammacoef(1)*Bn1+gammacoef(2)*Bn2
	sigma=gammacoef(1)*Cn1+sign(1d0,Ka-Kb)*gammacoef(2)*Cn2
	sigmaR=gammacoef(1)*Cn1-sign(1d0,Ka-Kb)*gammacoef(2)*Cn2

	end subroutine integrals

        subroutine ByCdoble(a1coef,a2coef,c0coef,nmin,nmax)
!esta subrutina ya no la usa!!!!!
!calcula los coeficientes B y C para las integrales rho tau sigma y sigma raya
	use ECP_mod, only : Bn1,Bn2,Cn1,Cn2,DAW,DAWERF,NEXTCOEF,pi,pi12
	IMPLICIT NONE
!a1coef,a2coef,c0coef son los coeficientes para el calculo de B y C
! Ai(c0coef,aicoef)= int exp(-c0coef*x^2)(x+aicoef)^n dx from -aicoef to inf
!nmin,nmax delimitan el rango de n
!Bni=An(c0coef,aicoef) + An(c0coef,-aicoef)
!Cni=An(c0coef,aicoef) - An(c0coef,-aicoef)
	DOUBLE PRECISION, intent(in) :: a1coef,a2coef,c0coef
	integer, intent(in) :: nmin,nmax
!variables auxiliares
	DOUBLE PRECISION :: c0sq,ncos,ca1,ca2
	integer :: i
	write(*,*) "byc",nmin,nmax
	c0sq=sqrt(c0coef)
	Bn1(0)=pi12/c0sq
	Bn2(0)=Bn1(0)
	Cn1(0)=Bn1(0)*erf(a1coef*c0sq)
	Cn2(0)=Bn1(0)*erf(a2coef*c0sq)

	if (nmax>0) then
	        Bn1(1)= exp(-c0coef*a1coef**2)/c0coef + a1coef*Cn1(0)
	        Bn2(1)= exp(-c0coef*a2coef**2)/c0coef + a2coef*Cn2(0)
	        Cn1(1)= Bn1(0)*a1coef
	        Cn2(1)= Bn1(0)*a2coef
		do i=2,nmax
			ncos=(i-1)/(2*c0coef)
			Bn1(i)=ncos*Bn1(i-2)+a1coef*Cn1(i-1)
			Bn2(i)=ncos*Bn2(i-2)+a2coef*Cn2(i-1)
			Cn1(i)=ncos*Cn1(i-2)+a1coef*Bn1(i-1)
			Cn2(i)=ncos*Cn2(i-2)+a2coef*Bn2(i-1)
		end do
	end if

        if (nmin<0) then
		Bn1(-1)=2*pi12*DAWERF(a1coef*c0sq)
		Bn2(-1)=2*pi12*DAWERF(a2coef*c0sq)
		write(*,*) "DAWERF(a2coef*c0sq)",DAWERF(a2coef*c0sq)
		Cn1(-1)=2*pi12*DAW(a1coef*c0sq)
		Cn2(-1)=2*pi12*DAW(a2coef*c0sq)
		if (nmin<-1) then
			ca1=2*c0coef*a1coef
			ca2=2*c0coef*a2coef
			Bn1(-2)=ca1*Cn1(-1)-2*c0coef*Bn1(0)
			Bn2(-2)=ca2*Cn2(-1)-2*c0coef*Bn2(0)
			Cn1(-2)=2*ca1*exp(-c0coef*a1coef**2)+ca1*Bn1(-1)-2*c0coef*Cn1(0)
			Cn2(-2)=2*ca2*exp(-c0coef*a2coef**2)+ca2*Bn2(-1)-2*c0coef*Cn2(0)
			do i=-3,nmin,-1
				Bn1(i)=NEXTCOEF(1,i,ca1,exp(-c0coef*a1coef**2),c0coef,Cn1(i+1), Bn1(i+2))
				Bn2(i)=NEXTCOEF(1,i,ca2,exp(-c0coef*a2coef**2),c0coef,Cn2(i+1), Bn2(i+2))
				Cn1(i)=NEXTCOEF(-1,i,ca1,exp(-c0coef*a1coef**2),c0coef,Bn1(i+1), Cn1(i+2))
				Cn2(i)=NEXTCOEF(-1,i,ca2,exp(-c0coef*a2coef**2),c0coef,Bn2(i+1), Cn2(i+2))
			end do
		end if
	end if

        end subroutine ByCdoble

	subroutine ReasignZ
!cambia la carga del nucleo sacandole la carga del core y guarda las cargas originales en IzECP
	use garcha_mod, only :Iz,nuc,nshell, natom
	use ECP_mod, only : ZlistECP,IzECP,Zcore,ecptypes
	implicit none
	integer :: i,j
	allocate (IzECP(natom))
	do i=1,nshell(0)+nshell(1)+nshell(2) !sumas mas si se agregan funciones f,g,etc
		IzECP(nuc(i))=Iz(nuc(i))
		write(*,*) "carga antes", Iz(i)
		do j=1, ecptypes
			if (Iz(i) .eq. ZlistECP(j)) then
				Iz(i)=Iz(i)-Zcore(j)
				write(*,*) "carga nueva" , Iz(i)
				write(*,*) "caraga guardada", IzECP(nuc(i))
			end if
		end do
	end do
	end subroutine ReasignZ

!	9018 format(/1x,'Z =',i4,2x, 'L =',i4,2x,'coefnumber =',i3,2x, 'n =',i2,2x, 'b =', f15.5,2x,'c =',f15.5)
!        end subroutine intECP
