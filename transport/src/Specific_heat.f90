! In this module, vibrational energy, non-equilibrium vibrational partition
! functions, vibrational specific heat capacities, and mean vibrational energy
! en_int of molecular species are calculated. It uses the module constant.f90
! containing main constants and variables definition
! Input variables: T, Tv_CO2_1, Tv_CO2_2, Tv_CO2_3, Tv_CO, Tv_O2

MODULE Specific_heat

USE CONSTANT

IMPLICIT NONE

real(dp) :: zv_1, zv_2, zv_3, zv_co2, zv_o2, zv_co, &
            c_v_t1, c_v_t2, c_v_t3, c_v_co2, c_v_o2, c_v_co

real(dp), dimension(5) :: en_int

REAL(dp) :: ee1,ee2,ee3,ee4

INTEGER I1,I2,I3,IL,G

CONTAINS

! Calculation of vibrational energy levels for CO2, O2, CO
  SUBROUTINE ENERGY

    do i1=1,5
       en_int(i1)=0
    end do

    do i1=0,l1
      do i2=0,l2
        do i3=0,l3
          !en3(i1,i2,i3)=w1*i1+w2*i2+w3*i3!+wx11*i1*i1+wx12*i1*i2+wx13*i1*i3+!wx22*i2*i2+wx23*i2*i3+wx33*i3*i3+wxll*il*il
          en3(i1,i2,i3)=w1*i1+w2*i2+w3*i3+wx11*i1*i1+wx12*i1*i2+wx13*i1*i3+wx22*i2*i2+wx23*i2*i3+wx33*i3*i3+wxll*il*il
        end do
      end do
    end do
  
    do i1=0,l_o2
      !en_o2(i1)=we_o2*(i1)!-wexe_o2*i1**2
      en_o2(i1)=we_o2*(i1)-wexe_o2*i1**2
    end do
  
    do i1=0,l_co
      !en_co(i1)=we_co*(i1)!-wexe_co*i1**2
      en_co(i1)=we_co*(i1)-wexe_co*i1**2
    end do

  END SUBROUTINE ENERGY

! Calculation of non-equilibrium CO2 partition function Z_CO2(T1,T2,T3)
  SUBROUTINE PART_FUNC_CO2(T1,T2,T3)

    REAL(dp) :: EE, T1, T2, T3

    zv_1=0
    zv_2=0
    zv_3=0

    zv_co2=0

    do i1=0,L1
      ee=en3(i1,0,0)
      if(ee<d_co2*kb) then
        zv_1=zv_1+exp(-(i1*ee1)/kb/t1)
      end if
    end do

    do i2=0,L2
      g=i2+1
      ee=en3(0,i2,0)
      if(ee<d_co2*kb) then
        zv_2=zv_2+g*exp(-(i2*ee2)/kb/t2)
      end if
    end do

    do i3=0,L3
      ee=en3(0,0,i3)
      if(ee<d_co2*kb) then
        zv_3=zv_3+exp(-i3*ee3/kb/t3)
      end if
    end do

    zv_co2 = zv_1 * zv_2 * zv_3

  END SUBROUTINE PART_FUNC_CO2

! Calculation of non-equilibrium O2 partition function Z_O2(TVO2)
  SUBROUTINE PART_FUNC_O2(TVO2)

    INTEGER I1
  	REAL(dp) :: EE,TVO2

  	zv_O2=0

    DO i1=0,L_O2
      ee=en_O2(i1)
	  zv_O2=zv_O2+exp(-ee/kb/TVO2)
    end do

  END SUBROUTINE PART_FUNC_O2

! Calculation of non-equilibrium CO partition function Z_CO(TVCO)
  SUBROUTINE PART_FUNC_CO(TVCO)

  INTEGER I1
  REAL(dp) :: EE, TVCO

  zv_CO=0

  DO i1=0,L_CO
    ee=en_CO(i1)
	zv_CO=zv_CO+exp(-ee/kb/TVCO)
  end do

  END SUBROUTINE PART_FUNC_CO

! Calculation of non-equilibrium CO2 vibrational specific heats CVT1, CVT2, CVT3
  SUBROUTINE S_HEAT_CO2

    integer i1,i2,i3
    real(dp) :: ee
    real(dp) :: ppp,s1,s2,s3,ss1,ss2,ss3

    s1=0;
    s2=0;
    s3=0;
    ss1=0;
    ss2=0;
    ss3=0

    do i1=0,L1
      ee=en3(i1,0,0)
      if(ee<d_co2*kb) then
        ppp=exp(-(i1*ee1)/kb/t1)/zv_1
        s1=s1+(i1*ee1)/kb/t1*ppp
        ss1=ss1+((i1*ee1)/kb/t1)**2*ppp
        en_int(1)=en_int(1)+ee*ppp
      end if
    end do

    do i2=0,L2
      g=i2+1
      ee=en3(0,i2,0)
      if(ee<d_co2*kb) then
        ppp=g*exp(-(i2*ee2)/kb/t2)/zv_2
        s2=s2+(i2*ee2)/kb/t2*ppp
        ss2=ss2+((i2*ee2)/kb/t2)**2*ppp
        en_int(1)=en_int(1)+ee*ppp
      end if
    end do

    do i3=0,L3
      ee=en3(0,0,i3)
      if(ee<d_co2*kb) then
        ppp=exp(-i3*ee3/kb/t3)/zv_3
        s3=s3+i3*ee3/kb/t3*ppp
        ss3=ss3+i3*ee3/kb/t3*i3*ee3/kb/t3*ppp
        en_int(1)=en_int(1)+ee*ppp
      end if
    end do

    c_v_t1=(ss1-s1*s1)
    c_v_t2=(ss2-s2*s2)
    c_v_t3=(ss3-s3*s3)

    c_v_co2 = c_v_t1 * c_v_t2 * c_v_t3

  END SUBROUTINE S_HEAT_CO2

! Calculation of non-equilibrium O2 vibrational specific heat CV_O2
  SUBROUTINE s_heat_O2

	integer I1
	real(dp) :: ee
	real(dp) :: ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_O2
       ee=en_O2(i1)
       ppp=exp(-ee/kb/TVO2)/zv_O2
       s=s+ee/kb/TVO2*ppp;
       s0=s0+ee*ee/kb/TVO2/kb/TVO2*ppp;
	   en_int(2)=en_int(2)+ee*ppp
    END DO
    c_v_o2=(s0-s*s)

  END SUBROUTINE s_heat_O2

! Calculation of non-equilibrium CO vibrational specific heat CV_CO
  SUBROUTINE s_heat_CO

	integer I1
	real(dp) :: ee
	real(dp) :: ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_CO
       ee=en_CO(i1)
       ppp=exp(-ee/kb/TVCO)/zv_CO
       s=s+ee/kb/TVCO*ppp
       s0=s0+ee*ee/kb/TVCO/kb/TVCO*ppp
       en_int(3)=en_int(3)+ee*ppp
    END DO
    c_v_co=(s0-s*s)

  END SUBROUTINE s_heat_CO

END MODULE Specific_heat
