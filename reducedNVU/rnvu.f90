!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   rnvu :     The Reduced NVU model with initial JPLC=0
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

!      DOUBLE PRECISION Ca_i,AMp,AM,R,Mp,s_i,v_i,w_i,I_i,Ca_j,s_j,v_j,I_j,J_PLC,K_p
      DOUBLE PRECISION R,AMp,AM,Ca_i,Mp,s_i,v_i,w_i,I_i,Ca_j,s_j,v_j,I_j,J_PLC,K_p

      Ca_i=U(5)
      Mp=U(2)
      AMp=U(3)
      AM=U(4)
      R=U(1)
      s_i=U(6)
      v_i=U(7)
      w_i=U(8)
      I_i=U(9)
      Ca_j=U(10)
      s_j=U(11)
      v_j=U(12)
      I_j=U(13)

      J_PLC=PAR(1)
      K_p=PAR(2)

      F(5) = (0.23d0 * I_i**2 / (1.d0**2 + I_i**2)) - (2.025d0 * Ca_i**2 / (1.0d0**2 + Ca_i**2)) - &
      (0.24d0 * Ca_i * (1.d0 + (v_i + 100.d0) / 250.d0)) + (0.025d0 * s_i) - & 
      ((1.29E-3) * (v_i - 100.d0) /  (1.d0 + exp(-(v_i - (-24.d0)) / 8.5d0))) + &
      (55.d0 * s_i**2 / (2.0d0**2 + s_i**2) * Ca_i**4 / (0.9d0**4 + Ca_i**4)) + &
      ((3.16E-3) * Ca_i / (Ca_i + 0.5d0) * (v_i - (-30.d0))) + & 
      0.1d0*((6.1E-3) / (1 + exp(-(7.4E-3)*(30.d0*R/(0.1d0 * R) - 500.d0))) * (v_i - (-18.d0))) &
      + (-0.05d0 * (Ca_i - Ca_j)) 

      F(2) = 0.1d0 * AMp + (17.d0 * Ca_i**3) * (1.d0 - AM - AMp - Mp) - (0.5d0 + 0.4d0)*Mp    
                                                                                                                                                                                                                                                                                                                                                                                                    
      F(3) = 0.4d0*Mp + (17.d0*Ca_i**3)*AM - (0.1d0 + 0.5d0)*AMp   
                                                                                                                                                                                                                                                                                                                                                                                                                     
      F(4) = 0.5d0* AMp - (0.1d0 + (17.d0* Ca_i**3)) * AM      
                                                                                                                                                                                                                                                                                                                                                                                                                                
      F(1) = (20E-6) / (1E4) * ( R * (4000.d0) / (0.1d0 * R) - ((66E3) + (AMp + AM)* ((233E3) - (66E3))) * (R - ((20E-6) &
      + (AMp + AM)* ((0.6d0) - 1.d0) * (20E-6))) / ((20E-6) + (AMp + AM)* ((0.6d0) - 1.d0) * (20E-6)))   
                 
      F(6) = (2.025d0 * Ca_i**2 / (1.0d0**2 + Ca_i**2)) - &
      (55.d0 * s_i**2 / (2.0d0**2 + s_i**2) * Ca_i**4 / (0.9d0**4 + Ca_i**4)) - (0.025d0 * s_i)         
                                                                                                                                                                                                                                                                                                                                                
      F(7) = 1970.d0 * (- ((4.32E-2)) - ( (1.34E-3) * (v_i - (-25.d0))) - 2.d0*((1.29E-3) * (v_i - &
      100.d0) / (1.d0 + exp(-(v_i - (-24.d0)) / 8.5d0))) - ((3.16E-3) * Ca_i / (Ca_i + 0.5d0) * (v_i - (-30.d0))) &
      - ((4.46E-3) * w_i * (v_i - (-94.d0))) - ((6.1E-3) / (1.d0 + exp(-(7.4E-3)*(30.d0*R/(0.1d0 * R) &
      - 500.d0))) * (v_i - (-18.d0))) - ((7.5E2) * (exp((-7.4E-2) * v_i + (4.2E-4) * K_p - 12.6d0)) / 1970.d0 * (v_i &
      - ((4.5E-3) * K_p - 112.d0)))) + (-0.5d0 * (v_i - v_j))

      F(8) = 45.d0 * (((Ca_i + 0.d0)**2 / ((Ca_i + 0.d0)**2 + 0.13d0 * exp(-(v_i - (-27.d0)) / 12.d0))) - w_i)        
                                                                                                                                                                                                                                                                                                                                                                                 
      F(9) = (-0.05d0 * (I_i - I_j)) - (0.1d0 * I_i)      
                                                                                                                                                                                                                                                                                                                                                                                                                                
      F(10) = (0.23d0 * I_j**2 / (1.d0**2 + I_j**2)) - (0.5d0 * Ca_j**2 / (1.d0**2 + Ca_j**2)) + &
      (5.d0 * s_j**2 / (2.d0**2 + s_j**2) * Ca_j**4 / (0.9d0**4 + Ca_j**4)) - (0.24d0 * Ca_j) + &
      (0.025d0 * s_j) + ((6.6E-4) * (50.d0 - v_j) * 0.5d0 * (1 + tanh((log10(Ca_j) - (-0.18d0)) / 0.37d0))) &
      + 0.029d0 + ((6.1E-3) / (1.d0 + exp(-(7.4E-3)*(30.d0*R/(0.1d0 * R) - 500.d0))) * (v_j - (-18))) &
      - (-0.05d0 * (Ca_i - Ca_j))     
                                                                                                
      F(11) = (0.5d0 * Ca_j**2 / (1.d0**2 + Ca_j**2)) - (5.d0 * s_j**2 / (2.d0**2 &
      + s_j**2) * Ca_j**4 / (0.9d0**4 + Ca_j**4)) - (0.025d0 * s_j)     
                                                                                                                                                                                                                                                                                                                                                           
      F(12) = -1/25.8d0 * ((6927.d0 * (v_j - (-80.d0)) * ((0.2d0 * (1.d0 + tanh( ((log10(Ca_j) &
      - (-0.4d0)) * (v_j - (-80.8d0)) - 53.3d0) / ((1.32E-3) * (v_j + 53.3d0*(log10(Ca_j) - (-0.4d0)) &
      - (-80.8d0))**2 + 0.3d0)))) + (0.3d0 * (1 + tanh((log10(Ca_j) - (-0.28d0)) / 0.389d0))))) &
      + (955.d0 * (v_j - (-31.1d0)))) - (-0.5d0 * (v_i - v_j))    
                                                                                                                                                                      
      F(13) = J_PLC - (0.1d0 * I_j) - (-0.05d0 * (I_i - I_j))                 


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.1d0
      PAR(2)=9200.d0

      U(5)=1.d0
      U(2)=1.d0
      U(3)=1.d0
      U(4)=1.d0
      U(1)=1.d0
      U(6)=1.d0
      U(7)=1.d0
      U(8)=1.d0
      U(9)=1.d0
      U(10)=1.d0
      U(11)=1.d0
      U(12)=1.d0
      U(13)=1.d0

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
