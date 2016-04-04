!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   ms :     The Meyer & Stryer model
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION Cai, IP3, Cas, R, K_1, K_2, K_3, c_5, c_7, c_1, c_2, c_3, c_4, c_6, J_1, J_2

      Cai=U(1)
      IP3=U(2)
      Cas=U(3)

      R=PAR(1)

      K_1 = 0.1d0;         
      K_2 = 0.15d0;       
      K_3 = 1.d0;             

      c_5 = 2.d0;           
      c_7 = 0.6d0;          

      c_1 = 6.64d0/c_5;
      c_2 = 5.d0/(c_5*c_7);             
      c_3 = (c_7/c_5)*3.13E-5;  
      c_4 = 1.d0/(c_5*c_7);                         
      c_6 = 0.5d0/(c_5*c_7);  

      J_1 = c_1*Cas*(c_7**3*IP3**3/(K_1 + c_7*IP3)**3);         
      J_2 = c_2*((c_7*Cai)**2)/(K_2 + c_7*Cai)**2 - c_3*Cas**2; 

      F(1) = J_1 - J_2 + c_6*(1 - Cai**3.3d0);        
      F(2) = c_4*R*c_7*Cai/(c_7*Cai + K_3) - IP3;
      F(3) = J_2 - J_1;

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.d0

      U(1)=1.d0
      U(2)=0.d0
      U(3)=532.9074

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
