!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   ermentrout :     The Gonzalez Fernandez & Ermentrout model
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION Ca, V, N, v1, vca, vL, vK, v2, v4, v5, v6, ca3, ca4, gca, gL, gK, phin
      DOUBLE PRECISION kca, kd, bt, alphgca, C, unbufca, minf, v3, ninf

      Ca = U(1)
      V = U(2)
      N = U(3)

      v1 = PAR(1)

      vca = 80.d0;  
      vL = (-70.d0/vca);
      vK = (-90.d0/vca);
      v2 = (25.d0/vca);
      v4 = (14.5d0/vca);
      v5 = (8.d0/vca);
      v6 = (-15.d0/vca);
      ca3 = 400.d0;
      ca4 = (150.d0/ca3);
      gca = (1.57E-13); 
      gL = ((7.854E-14)/gca);
      gK = ((3.1416E-13)/gca);
      phin = 2.664d0;
      kca = (135.67537d0/phin);
      kd = 1000.d0;
      bt = 100000.d0;
      alphgca = 1255.6232d0;
      C = (1.9635E-14);

      unbufca = ((kd+(Ca*ca3))**2)/((kd+(Ca*ca3))**2+kd*bt);
      minf = (0.5d0*(1+tanh((V-v1)/v2)));
      v3 = (-(v5/2.d0)*tanh((Ca-1.d0)/ca4)+v6);
      ninf = 0.5d0*(1+tanh((V-v3)/v4));

      F(1) = unbufca*(((-alphgca*vca*minf)/(ca3*phin))*(V-1.d0)-kca*Ca);
      F(2) = -(gca/(C*phin))*(gL*(V-vL)+gK*N*(V-vK)+minf*(V-1.d0));
      F(3) = cosh((V-v3)/(2.d0*v4))*(ninf-N);

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=-0.5d0

      U(1)=1.916d0
      U(2)=-0.2646d0
      U(3)=0.4239d0

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
