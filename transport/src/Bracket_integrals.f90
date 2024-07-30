!
! In this module, bracket integrals in the transport linear systems are calculated.
!
! It uses modules:
! - constant.f90
! - Specific_heat.f90
! - Omega_integrals.f90
!
! Input variables: X(i), T

MODULE BRACKET_INTEGRALS

USE CONSTANT
USE OMEGA_INTEGRALS
USE SPECIFIC_HEAT

IMPLICIT NONE

CONTAINS

  SUBROUTINE BRACKET

	INTEGER i, j, k
   	REAL(dp) :: mij
	REAL(dp), DIMENSION(5) :: PHI=(/0., 0., 0., 0., 0./)

	DO i=1,5
		DO j=1,5
		mij=mass(j)*mass(i)/(mass(i)+mass(j))
		lambda(i,j)=75./64.*t*kb/mij/OMEGA22(i,j)*kb
		eta(i,j)=5./8.*kb*t/OMEGA22(i,j)
		END DO
	END DO

	PHI(1)=KB/23.73*(1.+PI**(3/2)/2.*SQRT(EPS(1)/T)+(PI*PI/4.+2.)*EPS(1)/T+ &
			(PI*EPS(1)/T)**(3/2))
	PHI(2)=KB/20.72*(1.+PI**(3/2)/2.*SQRT(EPS(2)/T)+(PI*PI/4.+2.)*EPS(2)/T+ &
			(PI*EPS(2)/T)**(3/2))
			!+KB*C_V_O2*PI*ETA(2,2)/4./(EXP(0.00116* &
			!SQRT(MASS(2)/AMU/2)*(en_O2(1)/KB)**1.333*(T**(-1./3.)- &
!			!0.015*(MASS(2)/AMU/2)**(0.25))-18.42)*1.013E5)
	PHI(3)=KB/9.16*(1.+PI**(3/2)/2.*SQRT(EPS(3)/T)+(PI*PI/4.+2.)*EPS(3)/T+&
			(PI*EPS(3)/T)**(3/2))
			!+KB*C_V_CO*PI*ETA(3,3)/4./(EXP(0.00116* &
			!SQRT(MASS(3)/AMU/2)*(en_CO(1)/KB)**1.333*(T**(-1./3.)- &
			!0.015*(MASS(3)/AMU/2)**(0.25))-18.42)*1.013E5)

	DO i=1,5
	lambda_int(i)=0
	END DO

	DO i=1,3
	beta0011(i)=0
	END DO

	DO i=1,5
		DO j=1,5
		mij=mass(j)*mass(i)/(mass(i)+mass(j))
		if(i==j) then
			lambda00(i,j)=0
			lambda01(i,j)=0
			lambda11(i,j)=(x(i))**2/lambda(i,i)
			h00(i,j)=(x(i))**2/eta(i,i)
			beta11(i,j)=4.*t/pi*x(i)**2/eta(i,i)*phi(i)
			DO k=1,5
				if(k.NE.i) then
				lambda00(i,j)=lambda00(i,j)+x(i)*x(k)/ &
				lambda(i,k)/2./aa(i,k)

				lambda01(i,j)=lambda01(i,j)-x(i)*x(k)/ &
				lambda(i,k)/4./aa(i,k)*mass(k)/(mass(i)+mass(k))* &
				(6.*cc(i,k)-5.)

				lambda11(i,j)=lambda11(i,j)+x(i)*x(k)/lambda(i,k)/ &
				2./aa(i,k)*(15./2.*(mass(i))**2+25./4. &
				*(mass(k))**2-3.*(mass(k))**2*bb(i,k)+ &
				4.*mass(i)*mass(k)*aa(i,k))/(mass(i)+mass(k))**2

				h00(i,j)=h00(i,j)+2.*x(i)*x(k)/eta(i,k)* &
				mass(i)*mass(k)/(mass(i)+mass(k))**2* &
				(5./3./aa(i,k)+mass(k)/mass(i))

				beta11(i,j)=beta11(i,j)+x(i)*x(k)/eta(i,k)* &
				mass(k)/(mass(i)+mass(k))**2*(5.*kb*t*mass(i)/aa(i,k)+ &
				4.*t*mass(k)/pi*(phi(i)+phi(k)))

				end if
			end do
		else
			lambda00(i,j)=-x(i)*x(j)/lambda(i,j)/2./aa(i,j)

			lambda01(i,j)=x(i)*x(j)/lambda(i,j)/4./aa(i,j)* &
			mass(i)/(mass(i)+mass(j))*(6.*cc(i,j)-5.)

			lambda11(i,j)=-x(i)*x(j)/lambda(i,j)/2./aa(i,j)* &
			mij/(mass(i)+mass(j))*(55./4.-3.*bb(i,j)-4.*aa(i,j))

			h00(i,j)=-2.*x(i)*x(j)/eta(i,j)*mij/(mass(i)+mass(j))* &
			(5./3./aa(i,j)-1)

			beta11(i,j)=x(i)*x(j)/eta(i,j)* &
			mass(i)*mass(j)/(mass(i)+mass(j))**2*(-5.*kb*t/aa(i,j)+ &
			4.*t/pi*(phi(i)+phi(j)))

		end if

		lambda_int(i)=lambda_int(i)+x(j)*mij*omega11(i,j)

		END DO
	END DO

  	DO i=1,5
		DO j=1,3
		if(i==j) then
			beta01(i,j)=-4.*t/pi*x(i)**2/eta(i,i)*phi(i)
			DO k=1,5
				if(k.NE.i) then
				beta01(i,j)=beta01(i,j)-4.*t/pi*x(i)*x(k)/eta(i,k)* &
				mass(k)/(mass(i)+mass(k))*phi(i)
				end if
			end do
		else
			beta01(i,j)=-4.*t/pi*x(i)*x(j)/eta(i,j)* &
			mass(j)/(mass(i)+mass(j))*phi(j)
		end if

		beta0011(j)=beta0011(j)+4.*t/pi*x(i)*phi(j)*x(j)/eta(i,j)

		END DO
	END DO

  END SUBROUTINE BRACKET

END MODULE BRACKET_INTEGRALS
