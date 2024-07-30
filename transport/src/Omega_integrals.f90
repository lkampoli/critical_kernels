!
! In this module, Omega-integrals and their ratios AA, BB, CC
! are calculated using the Lennard-Jones potential for
! moderate temperatures and the repulsive (Born-Meyer) potential
! for high temperatures.
! It uses the module constant.f90 containing main constants and
! variables definition
!
! Input variable: T

MODULE OMEGA_INTEGRALS

USE CONSTANT

IMPLICIT NONE

CONTAINS

  SUBROUTINE OMEGA

  ! Calculation of OMEGA-integrals for given T,
  ! repulsive potential

  INTEGER i, j
  REAL(dp) :: x11, eij, mij, sig_ij, sij, tx, ax, a10, a10r, a2, a4, a6, deriv1

  REAL(dp) :: RX(5,5), VX(5,5)

 !Definition of the repulsive potential parmeters; data given by
 ![1] J.Bzowski, J.Kestin, E.A.Mason, F.J. Uribe, J. Phys. Chem. Ref. Data,
 !19 (5), 1179-1232, 1990
 ![2] V.V. Riabov, Journ. Thermophys. Heat Transf, 10(2), 209-216, 1996

 !Matrix RX

 RX(1,1)=7.2489E-0002
 RX(2,1)=0.072477318
 RX(3,1)=7.3836E-0002
 RX(4,1)=7.9164E-0002
 RX(5,1)=0.087723146
 RX(2,2)=0.072952699
 RX(3,2)=0.073889995
 RX(4,2)=0.079763543
 RX(5,2)=0.088047267
 RX(3,3)=7.5143E-0002
 RX(4,3)=8.0500E-0002
 RX(5,3)=0.089387872
 RX(4,4)=8.7835E-0002
 RX(5,4)=0.096840576
 RX(5,5)=0.107903965

   DO I=1,5
	DO J=2,5
	  if(i<j) then
	  RX(i,j)=RX(j,i)
	  endif
	END DO
  END DO

 !Matrix VX

 VX(1,1)=4.9896E+0006
 VX(2,1)=3.0763E+06
 VX(3,1)=1.8722E+0006
 VX(4,1)=1.1712E+0006
 VX(5,1)=4.3296E+05
 VX(2,2)=1.6047E+06
 VX(3,2)=1.0972E+06
 VX(4,2)=5.8993E+05
 VX(5,2)=2.6006E+05
 VX(3,3)=7.4296E+0005
 VX(4,3)=4.3051E+0005
 VX(5,3)=1.9028E+05
 VX(4,4)=2.0455E+0005
 VX(5,4)=1.0162E+05
 VX(5,5)=5.0486E+04

   DO I=1,5
	DO J=2,5
	  if(i<j) then
	  VX(i,j)=VX(j,i)
	  endif
	END DO
  END DO

 !Parameters sigma(i) and eps(i) of the Lenard-Jones potential are
 !defined in module constant.f90. Data given by:
 !R.J.Kee, J.A.Miller, T.N. Jefferson, CHEMKIN: A General-Purpose,
 !Problem-Independent, Transportable, Fortran Chemical Kinetics Code Package,
 !Sandia National Laboratories, SAND80-8003, 1980

  DO i=1,5
	DO j=1,5
  		sig_ij=(sigma(i)+sigma(j))/2
		mij=mass(j)*mass(i)/(mass(i)+mass(j))
		eij=sqrt(eps(i)*eps(j)*(sigma(i)*1e10)**6 &
			*(sigma(j)*1e10)**6)/(sig_ij*1e10)**6
		sij=pi*(sig_ij)**2*sqrt(kb*t/2./pi/mij)

		tx=t/eij

		if(tx>10) then
			! Omega11, Omega22, B, C for high temperature
			! conditions; repulsive potential
			ax=log(vx(i,j))-log(tx)
			a10=log(vx(i,j)/10)
			a10r=a10*rx(i,j)
			a2=-267.00+1/(a10r)**2*(201.570+174.672/a10+ &
			(7.36916/a10)**2);
			a4=26700-1/(a10r)**2*(19.2265+27.6938/a10+ &
			(3.29559/a10)**2)*1e3;
			a6=-8.90e5+1/(a10r)**2*(6.31013+10.2266/a10+ &
			(2.33033/a10)**2)*1e5;
			omega11(i,j)=(ax*rx(i,j))**2*(0.89+a2/(tx)**2+ &
			a4/(tx)**4+a6/(tx)**6)*sij

			a2=-33.0838+1/(a10r)**2*(20.0862+72.1059/a10+ &
			(8.27648/a10)**2)
			a4=101.571-1/(a10r)**2*(56.4472+286.393/a10+ &
			(17.7610/a10)**2);
			a6=-87.7036+1/(a10r)**2*(46.3130+277.146/a10+ &
			(19.0573/a10)**2)
			omega22(i,j)=(ax*rx(i,j))**2*(1.04+a2/(log(tx))**2+ &
			a4/(log(tx))**3+a6/(log(tx))**4)*sij*2

			deriv1=1/omega11(i,j)*sij*((ax*rx(i,j))**2* &
			(-2.*a2/tx**3-4.*a4/tx**5-6.*a6/tx**7)-2.*rx(i,j)**2/tx* &
			ax*(0.89+a2/tx**2+a4/tx**4+a6/tx**6))

			cc(i,j)=1+tx/3.*deriv1

			bb(i,j)=1+3.*cc(i,j)-3.*cc(i,j)**2-tx**2/3.* &
			(1/omega11(i,j)*sij*rx(i,j)**2*(-4./tx*ax* &
			(-2.*a2/tx**3-4.*a4/tx**5-6.*a6/tx**7)+ax**2* &
			(6.*a2/tx**4+20.*a4/tx**6-42.*a6/tx**8)+ &
			2./tx**2*(ax+1)*(0.89+a2/tx**2+a4/tx**4+a6/tx**6))- &
			deriv1**2)

		else
			! Omega11, Omega12, Omega12, Omega13, B, C for
			! moderate temperatures; Lennard-Jones potential
			x11=log((tx))+1.4
			omega11(i,j)=1/(-0.16845-2.25768e-2/x11/x11+ &
			0.19779/x11+0.64373*x11-9.26718e-2*x11*x11+ &
			7.1131e-3*x11**3)*sij

			x11=log((tx))+1.5
			omega22(i,j)=1/(-0.40811-5.08552e-2/x11/x11+ &
			0.3401/x11+0.70375*x11-0.10699*x11*x11+ &
			7.62686e-3*x11**3)*sij*2.

			x11=log((tx))+1.1
			omega12(i,j)=1/(0.40785+9.25303e-4/x11/x11+ &
			2.79680e-4/x11+0.44739*x11-6.27242e-2*x11*x11+ &
			5.98567e-3*x11**3)*sij*3.

			x11=log((tx))+4.0
			omega13(i,j)=1/(25.04929+63.12444/x11/x11- &
			65.87398/x11-4.13758*x11+0.34999*x11*x11- &
			1.0096e-2*x11**3)*sij*12.

			bb(i,j)=(5./3.*omega12(i,j)-4./12.*omega13(i,j)) &
			/omega11(i,j)

			cc(i,j)=1/omega11(i,j)*omega12(i,j)/3.0

		endif

		aa(i,j)=1/omega11(i,j)*omega22(i,j)/2.0

	END DO
  END DO

  END SUBROUTINE OMEGA

END MODULE OMEGA_INTEGRALS
