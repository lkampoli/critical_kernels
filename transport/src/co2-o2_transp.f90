
! Main program for the calculation of transport coefficients in a CO2/O2/CO/O/C mixture.
!
! Six-temperature model (T, T1, T2, T3, TVO2, TVCO).
!
! The program uses modules:
!
! - Constant.f90
! - Specific_heat.f90
! - Omega_integrals.f90
! - Bracket_integrals.f90
! - invers.f90
!
! Input variables: X(i), T, T1, T2, T3, TVO2, TVCO, press (not in this order)
!
! Output:
! - heat conductivity coefficient
! - vibrational heat conductivity coefficients
! - CO2, symmetric, bending, asymmetric mode; O2 and CO vibrational heat conductivity
! - shear viscosity coefficient
! - bulk viscosity coefficient
! - diffusion and thermal diffusion coefficients

PROGRAM transport

USE CONSTANT
USE SPECIFIC_HEAT
USE OMEGA_INTEGRALS
USE BRACKET_INTEGRALS

IMPLICIT NONE

INTEGER I,J,K,DELTA

! Matrices for the linear transport systems defining:
! - heat conductivity and thermal diffusion (LTH);
! - bulk viscosity (BVISC);
! - diffusion (LDIFF);
! - shear viscisity (HVISC).

REAL(dp), DIMENSION(10,10) :: LTH

REAL(dp), DIMENSION(8,8) :: BVISC

REAL(dp), DIMENSION(5,5) ::  LDIFF, HVISC, b1

! Vectors of right hand terms

REAL(dp), DIMENSION(10,1) :: b

REAL(dp), DIMENSION(5,1) :: b2

REAL(dp), DIMENSION(8,1) :: b3

REAL(dp) :: CU, CUT

! Neural network stuff ...
integer :: iT, iT1, iT2, iT3, iTVO2, iTVCO, iTsteps, iTstart, iTincr
integer :: iP, iPsteps, iPstart, nT, nX
integer :: ix1, ix2, ix3, ix4, ix5, iXsteps, iXstart, iXincr
integer :: nlines, dummy
real(dp):: T_tr, T_v1, T_v2, T_v3, T_O2, T_CO, TTT(6)
real(dp):: P, PPP, incrT, incrX, dT
real(dp):: x_CO2, x_CO, x_O2, x_O, x_C, XXX(5)
real(dp):: input(12) ! Assuming 12 input features, as written above

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!real(dp), dimension(10) :: Ttr,Tv1CO2,Tv2CO2,Tv3CO2,TvibO2,TvibCO
!real(dp), dimension(5) :: Ttr,Tv1CO2,Tv2CO2,Tv3CO2,TvibO2,TvibCO
real(dp), dimension(4) :: Ttr,Tv1CO2,Tv2CO2,Tv3CO2,TvibO2,TvibCO
!real(dp), dimension(1) :: Ttr,Tv1CO2,Tv2CO2,Tv3CO2,TvibO2,TvibCO

!Ttr    = (/2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000. /)
!Tv1CO2 = (/2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000. /)
!Tv2CO2 = (/2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000. /)
!Tv3CO2 = (/2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000. /)
!TvibO2 = (/2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000. /)
!TvibCO = (/2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000. /)
!
!Ttr    = (/ 2000., 4000., 6000., 8000., 10000. /)
!Tv1CO2 = (/ 2000., 4000., 6000., 8000., 10000. /)
!Tv2CO2 = (/ 2000., 4000., 6000., 8000., 10000. /)
!Tv3CO2 = (/ 2000., 4000., 6000., 8000., 10000. /)
!TvibO2 = (/ 2000., 4000., 6000., 8000., 10000. /)
!TvibCO = (/ 2000., 4000., 6000., 8000., 10000. /)

Ttr    = (/ 2500., 5000., 7500., 10000. /)
Tv1CO2 = (/ 2500., 5000., 7500., 10000. /)
Tv2CO2 = (/ 2500., 5000., 7500., 10000. /)
Tv3CO2 = (/ 2500., 5000., 7500., 10000. /)
TvibO2 = (/ 2500., 5000., 7500., 10000. /)
TvibCO = (/ 2500., 5000., 7500., 10000. /)

!Ttr    = 2500.
!Tv1CO2 = 2500.
!Tv2CO2 = 2500.
!Tv3CO2 = 2500.
!TvibO2 = 2500.
!TvibCO = 2500.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CALL ENERGY

ee1=en3(1,0,0)
ee2=en3(0,1,0)
ee3=en3(0,0,1)

! Output file which will be used for training the NN
!open (2, file = 'output_random.dat', status = 'unknown', action='write', position='append')
!open (2, file = 'output.dat', status = 'unknown', action='write', position='append')
!open (2, file = 'shear_viscosity.dat', status = 'unknown', action='write', position='append')
!TODO: uncomment back
!open (2, file = 'shear_viscosity_lite.dat', status = 'unknown', action='write', position='append')

! Input parameters: species molar fractions, temperatures, pressure

!do k=1,1000000

!! Molar fractions
!call execute_command_line ("./run_dirichlet.sh")
!open (5, file = 'molar_fractions.out', status = 'old')
!read(5,*) (XXX(i), i=1,5)
!close(5)
!
!x(1) = XXX(1)
!x(2) = XXX(2)
!x(3) = XXX(3)
!x(4) = XXX(4)
!x(5) = XXX(5)

!x(1)=0.2
!x(2)=0.2
!x(3)=0.2
!x(4)=0.2
!x(5)=0.2

!! Temperatures
!call execute_command_line ("./run_randomT.sh")
!open (3, file = 'random_temperatures.out', status = 'old')
!read(3,*) (TTT(i), i=1,6)
!close(3)
!T    = TTT(1)
!T1   = TTT(2)
!T2   = TTT(3)
!T3   = TTT(4)
!TVO2 = TTT(5)
!TVCO = TTT(6)

!T=500*k
!T12=5000
!T1=5000
!T2=5000
!T3=5000
!TVO2=5000
!TVCO=5000

!! Pressure
!call execute_command_line ("./run_randomP.sh")
!open (4, file = 'random_pressure.out', status = 'old')
!read(4,*) PPP
!close(4)
!press = PPP

press=101300.

! Check print
!write(*,"(12F15.4)") T, T1, T2, T3, TVO2, TVCO, press, x

!////////////////// LOOP ///////////////////////
!TODO: change nT and nX to run longer ...
nT    = 4 !1 !5 !10 !40  ! T = 500:500:20000 K
incrT = 2000. !500.
nX    = 2 !1 !10  ! X = 0.0:0.1:1.0
incrX = 0.2 !0.1

do iP = 1, 1
   press = 101325.
   do iT = 1, nT
      T = Ttr(iT) !incrT*iT 
      do iT1 = 1, nT
         !dT = T - incrT*iT1
         !if (dt > 0. .and. dt < 5000.) then
         T1 = Tv1CO2(iT1) !incrT*iT1
         do iT2 = 1, nT
            !dT = T - incrT*iT2
            !if (dt > 0. .and. dt < 5000.) then
            T2 = Tv2CO2(iT2) !incrT*iT2 
            do iT3 = 1, nT
               !dT = T - incrT*iT3
               !if (dt > 0. .and. dt < 5000.) then
               T3 = Tv3CO2(iT3) !incrT*iT3
               do iTVO2 = 1, nT
                  !dT = T - incrT*iTVO2
                  !if (dt > 0. .and. dt < 5000.) then
                  TVO2 = TvibO2(iTVO2) !incrT*iTVO2 
                  do iTVCO = 1, nT
                     !dT = T - incrT*iTVCO
                     !if (dt > 0. .and. dt < 5000.) then
                     TVCO = TvibCO(iTVCO) !incrT*iTVCO 
                     do ix1 = 1, nX ! x_CO2
                        x(1) = incrX*ix1
                        do ix2 = 1, nX ! x_CO
                           x(2) = incrX*ix2
                           do ix3 = 1, nX ! x_O2
                              x(3) = incrX*ix3
                              do ix4 = 1, nX ! x_O
                                 x(4) = incrX*ix4
                                 do ix5 = 1, nX ! x_C
                                    x(5) = incrX*ix5

ntot=press/kb/t

rho=0
do i=1,5
	rho=rho+x(i)*mass(i)*ntot
end do

! Calculation of vibrational energy, partition functions and specific heats
CALL PART_FUNC_CO2(T1,T2,T3)
CALL PART_FUNC_O2(TVO2)
CALL PART_FUNC_CO(TVCO)

CALL S_HEAT_CO2
CALL S_HEAT_O2
CALL S_HEAT_CO

! Calculation of bracket integrals
CALL OMEGA
CALL BRACKET

! Definition of matrix LTH for calculation of
! thermal conductivity and thermal diffuaion coefficients
! The system has a form	(1)
! LTH times a = b, a is the vector of unknowns
DO i=1,5
	DO j=1,5
		LTH(i,j)=Lambda00(i,j)
	END DO
END DO

DO i=1,5
	DO j=6,10
		LTH(i,j)=Lambda01(i,j-5)
	END DO
END DO

DO i=6,10
	DO j=1,5
		LTH(i,j)=LTH(j,i)
	END DO
END DO

DO i=6,10
	DO j=6,10
		LTH(i,j)=Lambda11(i-5,j-5)
	END DO
END DO

DO j=1,5
	LTH(1,j)=x(j)*mass(j)*ntot/rho
END DO

DO j=6,10
	LTH(1,j)=0.
END DO
! End of matrix LTH definition

! Definition of matrix LDIFF for calculation of
! diffuaion coefficients
! The system has a form	(2)
! LDIFF times D = B1, D a is the matrix of unknowns
DO i=1,5
	DO j=1,5
	  LDIFF(i,j)=LTH(i,j)
	END DO
END DO
! End of matrix LDIFF definition

! Definition of matrix HVISC for calculation of
! shear viscocity coefficient
! The system has a form	(3)
! HVISC times h = b2, h a is the vector of unknowns
DO i=1,5
	DO j=1,5
	  HVISC(i,j)=H00(i,j)
	END DO
END DO
! End of matrix HVISC definition

! Definition of matrix BVISC for calculation of
! bulk viscosity coefficients
! The system has a form	(4)
! BVISC times f = b3, f is the vector of unknowns
DO i=1,5
	DO j=1,5
		BVISC(i,j)=beta11(i,j)
	END DO
END DO

DO i=1,5
	DO j=6,8
		BVISC(i,j)=beta01(i,j-5)
	END DO
END DO

DO i=6,8
	DO j=1,5
		BVISC(i,j)=BVISC(j,i)
	END DO
END DO

DO i=6,8
	DO j=6,8
		BVISC(i,j)=0
	END DO
END DO

BVISC(6,6)=beta0011(1)
BVISC(7,7)=beta0011(2)
BVISC(8,8)=beta0011(3)

DO j=1,5
	BVISC(1,j)=x(j)*3./2.*kb*ntot/rho
END DO

DO j=6,6
	BVISC(1,j)=x(1)*kb/mass(1)
END DO

DO j=7,7
	BVISC(1,j)=x(2)*kb/mass(2)!*(1+c_v_o2)
END DO

DO j=8,8
	BVISC(1,j)=x(3)*kb/mass(3)!*(1+c_v_co)
END DO
! End of matrix BVISC definition

! Definition of vector b (right hand side of system (1))
DO i=1,5
	b(i,1)=0.
END DO
DO i=6,10
	b(i,1)=4./5./kb*x(i-5)
END DO
! End of vector b definition

! Definition of matrix b1 (right hand side of system (2))
DO i=1,5
	DO j=1,5
	if(i==j) then
		delta=1
	else
		delta=0
	end if
	B1(i,j)=8./25./kb*(delta-mass(i)*x(i)*ntot/rho);
	END DO
END DO

DO j=1,5
	B1(1,j)=0
END DO
! End of matrix b1 definition

! Definition of vector b2 (right hand side of system (3))
DO i=1,5
	b2(i,1)=2./kb/t*x(i)
END DO
! End of vector b2 definition

! Definition of vector b3 (right hand side of system (4))

!cu=kb*ntot/rho*(3./2.+x(1)+x(2)*(1+c_v_o2)+x(3)*(1+c_v_co))
!cut=kb*ntot/rho*(x(1)+x(2)*(1+c_v_o2)+x(3)*(1+c_v_co))

cu=kb*ntot/rho*(3./2.+x(1)+x(2)+x(3))
cut=kb*ntot/rho*(x(1)+x(2)+x(3))

DO i=1,5
	b3(i,1)=-x(i)*cut/cu
END DO

b3(6,1)=x(1)*ntot/rho*kb/cu

b3(7,1)=x(2)*ntot/rho*kb/cu !*(1+c_v_O2)

b3(8,1)=x(3)*ntot/rho*kb/cu !*(1+c_v_co)

b3(1,1)=0
! End of vector b3 definition

! Linear system solution using the Gauss method
! The solutions a, d, h, f are written to b, b1, b2, b3, respectively
CALL gaussj(LTH,10,10,b,1,1)
CALL gaussj(Ldiff,5,5,b1,5,5)
CALL gaussj(HVISC,5,5,b2,1,1)
CALL gaussj(BVISC,8,8,b3,1,1)

! Thermal diffusion coefficients THDIF(i)
DO i=1,5
	thdif(i)=-1./2./ntot*b(i,1)
END DO

! Thermal conductivity coefficient associated to translational energy, LTR
LTR=0
DO i=6,10
	LTR=LTR+5./4.*kb*x(i-5)*b(i,1)
END DO

! Thermal conductivity coefficients associated to internal energies
lrot_co2=3.*kb*t*x(1)/16./lambda_int(1)*kb
lrot_o2=3.*kb*t*x(2)/16./lambda_int(2)*kb
lrot_co=3.*kb*t*x(3)/16./lambda_int(3)*kb
lvibr_o2=3.*kb*t*x(2)/16./lambda_int(2)*kb*(c_v_o2)
lvibr_co=3.*kb*t*x(3)/16./lambda_int(3)*kb*(c_v_co)
lvibr_1=3.*kb*t*x(1)/16./lambda_int(1)*kb*(c_v_t1)
lvibr_2=3.*kb*t*x(1)/16./lambda_int(1)*kb*(c_v_t2)
lvibr_3=3.*kb*t*x(1)/16./lambda_int(1)*kb*(c_v_t3)

! Total thermal conductivity coefficient at the translational temperature gradient
ltot=ltr+lrot_co2+lrot_o2+lrot_co

! Diffusion coefficients DIFF(i,j)
DO i=1,5
	DO j=1,5
		diff(i,j)=1./2./ntot*b1(i,j)
	END DO
END DO

! Shear viscosity coefficient VISC
visc=0
DO i=1,5
	visc=visc+kb*t/2.*b2(i,1)*x(i)
END DO

! Bulk viscosity coefficient BULK_VISC
bulk_visc=0
DO i=1,5
	bulk_visc=bulk_visc-kb*t*b3(i,1)*x(i)
END DO

!! Output
!open(6,file='transport.txt',status='unknown')
!
!WRITE (6, *) 'Transport properties of CO2/O2/CO/O/C mixture, six-temperature model'
!WRITE (6, *) 'Harmonic oscillators'
!
!WRITE (6, *)
!
!WRITE (6, *) 'INPUT DATA:'
!
!WRITE (6, *)
!
!WRITE (6, *) 'Temperature, K        ',t
!WRITE (6, *) 'Temperature T12, K    ',t12
!WRITE (6, *) 'Temperature T1, K     ',t1
!WRITE (6, *) 'Temperature T2, K     ',t2
!WRITE (6, *) 'Temperature T3, K     ',t3
!WRITE (6, *) 'Temperature TVO2, K   ',tVO2
!WRITE (6, *) 'Temperature TVCO, K   ',tVCO
!WRITE (6, *) 'Pressure, Pa          ',press
!WRITE (6, *) 'CO2 molar fraction     ',x(1)
!WRITE (6, *) 'O2 molar fraction      ',x(2)
!WRITE (6, *) 'CO molar fraction      ',x(3)
!WRITE (6, *) 'O molar fraction       ',x(4)
!WRITE (6, *) 'C molar fraction       ',x(5)
!
!WRITE (6, *)
!
!WRITE (6, *)
!WRITE (6, *)
!
!WRITE (6, *) 'TRANSPORT COEFFICIENTS:'
!WRITE (6, *)
!
!WRITE (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', visc
!WRITE (6, '(1x, A45, E13.5)') 'Bulk viscosity coefficient, Pa.s              ', bulk_visc
!WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', ltot
!WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_O2, W/m/K     ', lvibr_O2
!WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_CO, W/m/K     ', lvibr_CO
!WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_12, W/m/K     ', lvibr_12
!WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_1, W/m/K      ', lvibr_1
!WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_2, W/m/K      ', lvibr_2
!WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_3, W/m/K      ', lvibr_3
!WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of CO2, m^2/s         ', THDIF(1)
!WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O2, m^2/s          ', THDIF(2)
!WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of CO, m^2/s          ', THDIF(3)
!WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O, m^2/s           ', THDIF(4)
!WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of C, m^2/s           ', THDIF(5)
!
!WRITE (6, *)
!WRITE (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
!WRITE (6, *)
!
!do i=1,5
!!WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
!WRITE (6, '(5e15.6)') (DIFF(i,j), j=1,5)
!end do

!WRITE(6,FMT='(F10.2,18E20.10)') T,ltot,lvibr_12,lvibr_3,lvibr_O2,lvibr_CO,visc,bulk_visc,kb*ntot/rho

! output file
!write(*,*) T, T1, T2, T3, TVO2, TVCO, press, x

!TODO: uncomment back
write(*,*) T, T1, T2, T3, TVO2, TVCO, x ! no pressure 

!write(2,"(50f15.7)") T, T1, T2, T3, TVO2, TVCO, press, x, &
!                     visc, bulk_visc, ltot, lvibr_O2, lvibr_CO, lvibr_1, lvibr_2, lvibr_3, &
!                     THDIF(1), THDIF(2), THDIF(3), THDIF(4), THDIF(5), &
!                     ((DIFF(i,j), i=1,5),j=1,5)
!write(2,"(46f15.7)") T, T1, T2, T3, TVO2, TVCO, press, x, &
!                     visc, bulk_visc, ltot, &
!                     THDIF(1), THDIF(2), THDIF(3), THDIF(4), THDIF(5), &
!                     ((DIFF(i,j), i=1,5),j=1,5)

!write(2,"(6e12.3, f9.1, 5f7.2, 1f15.7)") T, T1, T2, T3, TVO2, TVCO, press, x, visc
!TODO: uncomment back
!write(2,"(6e12.3, 5f7.2, 1f15.7)") T, T1, T2, T3, TVO2, TVCO, x, visc ! without pressure, scientific format for T's
!, bulk_visc, ltot, &
!                     THDIF(1), THDIF(2), THDIF(3), THDIF(4), THDIF(5), &
!                     ((DIFF(i,j), i=1,5),j=1,5)

!end do loop
enddo
enddo
enddo
enddo
enddo
!         endif
        enddo
!       endif
      enddo
!     endif
    enddo
!   endif
  enddo
! endif
enddo
enddo
enddo

close(2)

write(*,*) " The End "
END
