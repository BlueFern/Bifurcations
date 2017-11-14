!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   koenigberger :     
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

!      DOUBLE PRECISION Z,Y,V,N,I
      DOUBLE PRECISION Z,Y,V,N,I,mbeta

      Z=U(1)
	  Y=U(2)
	  V=U(3)
	  N=U(4)
	  I=U(5)
	  
      mbeta=PAR(1)

	  F(1) = 0.23d0 * I**2 / (1.d0**2 + I**2) -0.00129d0 * (V- 100.d0)/(1 + exp(-(V+24.d0)/8.5d0)) &
	 + 0.00316d0*Z*(V+40.d0)/(Z+0.5d0) - 2.025d0*(Z**2)/(Z**2 + 1.d0**2) &
	 + 55.d0*(Y**2/(Y**2 + 2.d0**2))*(Z**4/(0.9d0**4 + Z**4)) &
	 - 0.24d0*Z*(1+(V+100.d0)/250.d0) +0.025d0*Y
	  
	  F(2) = + 2.025d0*(Z**2)/(Z**2 + 1.d0**2) - 55.d0*(Y**2/(Y**2 &
	 + 2.d0**2))*(Z**4/(0.9d0**4 + Z**4)) - 0.025d0*Y
	  
	  F(3) = 1970.d0*(-0.0432d0 - 0.00134d0*(V+25.d0) -2*(0.00129d0 * (V &
	 - 100.d0)/(1 + exp(-(V+24.d0)/8.5d0))) - 0.00316d0*Z*(V+40.d0)/(Z+0.5d0) - 0.00446d0*N*(V+94.d0))
	  
	  F(4) = 45.d0*((Z+0.d0)**2/((Z+0.d0)**2 + 0.13d0*exp(-(V+27.d0)/12.d0)) - N)
	  
	  F(5) = mbeta*(0.05d0 - 0.17d0) + 0.17d0 + 0.d0*Z**2/(0.3d0**2 + Z**2) - 0.1d0*I
	      


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.d0

      U(1)=1.d0
      U(2)=1.d0
      U(3)=1.d0
      U(4)=1.d0
      U(5)=1.d0
      
      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
