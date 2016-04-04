!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   goldbeter :     The Goldbeter model
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION Z, Y, beta, v_2, v_3

      Z=U(1)
      Y=U(2)

      beta=PAR(1)

      
      v_2 = (65.d0 * Z**2) / ( 1.d0**2 + Z**2 );
      v_3 = (500.d0 * Y**2 * Z**4) / ( (2.d0**2 + Y**2) * (0.9d0**4 + Z**4) );

      F(1) = 1.d0 + (7.3d0/1.d0)*beta - v_2/1.d0 + v_3/1.d0 + 1.d0*Y - (10.d0/1.d0)*Z;
      F(2) = v_2/1.d0 - v_3/1.d0 - Y;

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.d0

      U(1)=0.1d0
      U(2)=0.63655561d0

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
