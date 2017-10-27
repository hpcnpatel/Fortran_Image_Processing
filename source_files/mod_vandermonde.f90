!*************************************************************************!
! PROGRAM IS TO TO BUILD VANDERMONDE MATRIX FOR QR DECOMPOSITION POLYNOMIAL.!
! AX=B IS GENERAL SYSTEM OF EQUATION AND IN ORDER TO SOLVE THIS SYSTEM OF !
! EQUATIONS, WE CONSTRUCT A(WHICH IS A VANDARMONDE MATRIX) BY CONSTRUCTING!
! ARRAYS NAMED ARR IN THE FOLLOWING PROGRAM.                               !
!*************************************************************************!
MODULE VANDERMONDE_MODULE

    USE FACTORIAL
    USE CONSTANTS_MODULE
  
IMPLICIT NONE

CONTAINS
!*****************************************************************************!
! GENERALISED LAPLACE EQUATION.......
!
!*****************************************************************************!
! PP   KK   JJ
! ---  ---  ---  (  (II-1)     (JJ-II)     (KK-JJ))
! \    \    \    ( X^       * Y^         *Z^      )
! /    /    /    (--------------------------------)
! ---  ---  ---  (   (II-1)!* (JJ-II)!* (KK-JJ)!  )
! KK=1 JJ=1 II=1
!*****************************************************************************!
SUBROUTINE VANDERMONDE_3D(V_ARR,LEN)

INTEGER                                    ::X,Y,Z
INTEGER                                    ::I,J,K,EL
INTEGER                                    ::NN,KK,MM
INTEGER                                    ::POLY_ORDER
INTEGER                                    ::TOT
INTEGER                                    ::ELE_IN_ROW
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)    ::V_ARR
INTEGER,DIMENSION(3),INTENT(INOUT)         ::LEN
INTEGER,DIMENSION(3)                       ::CTR
!=====================
!=====================
POLY_ORDER=2
ELE_IN_ROW=((POLY_ORDER+1) * (POLY_ORDER+2) * (POLY_ORDER+3))/6
!=====================
CTR=(LEN-1)/2
!=====================
! (N+1)*(N+2)*(N+3)/6=NUMBER OF ELEMENTS IN N TH ORDER POLYNOMIAL
!=====================
TOT=0

DO K=-CTR(3),CTR(3)
DO J=-CTR(2),CTR(2)
DO I=-CTR(1),CTR(1)
    !DO J=-CTR+ABS(K),CTR-ABS(K)
    !DO I=-CTR+ABS(K+J),CTR-ABS(K+J)
    !!DO J=-CTR+ABS(K),CTR-ABS(K)
    !!DO I=-CTR+ABS(K)+ABS(J),CTR-ABS(K)-ABS(J)
    TOT=TOT+1
END DO
END DO
END DO

        IF(ALLOCATED(V_ARR))DEALLOCATE(V_ARR)
        ALLOCATE(V_ARR(TOT,ELE_IN_ROW))

TOT=0
DO Z=-CTR(3),CTR(3)
DO Y=-CTR(2),CTR(2)
DO X=-CTR(1),CTR(1)
    !DO Y=-CTR+ABS(Z),CTR-ABS(Z)
    !DO X=-CTR+ABS(Z+Y),CTR-ABS(Z+Y)
    !DO Y=-CTR+ABS(Z),CTR-ABS(Z)
    !DO X=-CTR+ABS(Z)+ABS(Y),CTR-ABS(Z)-ABS(Y)
    TOT=TOT+1
    EL=0
    DO KK=1,POLY_ORDER+1
    DO MM=1,KK
    DO NN=1,MM
        EL=EL+1
        V_ARR(TOT,EL)=((X**(NN-1))*(Y**(MM-NN))*(Z**(KK-MM)))/&
                        &(FACT(NN-1)*FACT(MM-NN)*FACT(KK-MM))
    ENDDO
    ENDDO
    ENDDO
!    WRITE(*,'(*(F5.2))')V_ARR(TOT,1:EL)
END DO
END DO
END DO

END SUBROUTINE VANDERMONDE_3D
!****************************************************************************!
SUBROUTINE VANDERMONDE_2D(V_ARR,LEN)

INTEGER                                    ::X,Y
INTEGER                                    ::I,J
INTEGER                                    ::NN,MM
INTEGER                                    ::POLY_ORDER
INTEGER                                    ::TOT
INTEGER                                    ::ELE_IN_ROW
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)    ::V_ARR
INTEGER,DIMENSION(3),INTENT(INOUT)         ::LEN
INTEGER,DIMENSION(3)                       ::CTR
!=====================
!=====================
POLY_ORDER=2
ELE_IN_ROW=((POLY_ORDER+1) * (POLY_ORDER+2) )/2
!=====================
CTR=(LEN-1)/2
!=====================
! (N+1)*(N+2)/2=NUMBER OF ELEMENTS IN N TH ORDER POLYNOMIAL
!=====================
TOT=0

DO J=-CTR(2),CTR(2)
DO I=-CTR(1),CTR(1)
    !DO J=-CTR+ABS(K),CTR-ABS(K)
    !DO I=-CTR+ABS(K+J),CTR-ABS(K+J)
    !!DO J=-CTR+ABS(K),CTR-ABS(K)
    !!DO I=-CTR+ABS(K)+ABS(J),CTR-ABS(K)-ABS(J)
    TOT=TOT+1
END DO
END DO

        IF(ALLOCATED(V_ARR))DEALLOCATE(V_ARR)
        ALLOCATE(V_ARR(TOT,ELE_IN_ROW))

TOT=0
DO Y=-CTR(2),CTR(2)
DO X=-CTR(1),CTR(1)
    !DO Y=-CTR+ABS(Z),CTR-ABS(Z)
    !DO X=-CTR+ABS(Z+Y),CTR-ABS(Z+Y)
    !DO Y=-CTR+ABS(Z),CTR-ABS(Z)
    !DO X=-CTR+ABS(Z)+ABS(Y),CTR-ABS(Z)-ABS(Y)
    TOT=TOT+1
    !!V_arr(TOT,1:ELE_IN_ROW)=GENERAL_POLYNOMIAL(X,Y,POLY_ORDER,ELE_IN_ROW)
    V_ARR(TOT,1:ELE_IN_ROW)=(/(((((X**(NN-1))*(Y**(MM-NN)))/(FACT(NN-1)*FACT(MM-NN)))&
    &,NN=1,MM),MM=1,POLY_ORDER+1)/)
!    write(*,*)X,Y,real(V_arr(TOT,1:ELE_IN_ROW),4)
END DO
END DO

END SUBROUTINE VANDERMONDE_2D
!****************************************************************************!
SUBROUTINE VANDERMONDE_1D(V_ARR,LEN)

INTEGER                                    ::X
INTEGER                                    ::I
INTEGER                                    ::NN
INTEGER                                    ::POLY_ORDER
INTEGER                                    ::TOT
INTEGER                                    ::ELE_IN_ROW
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)    ::V_ARR
INTEGER,DIMENSION(3),INTENT(INOUT)         ::LEN
INTEGER,DIMENSION(3)                       ::CTR
!=====================
!=====================
POLY_ORDER=2
ELE_IN_ROW=POLY_ORDER+1
!=====================
CTR=(LEN-1)/2
!=====================
! (N+1)=NUMBER OF ELEMEMTNS IN N TH ORDER POLYNOMIAL
!=====================
TOT=0

DO I=-CTR(1),CTR(1)
    TOT=TOT+1
END DO

        IF(ALLOCATED(V_ARR))DEALLOCATE(V_ARR)
        ALLOCATE(V_ARR(TOT,ELE_IN_ROW))
TOT=0
DO X=-CTR(1),CTR(1)
    TOT=TOT+1
    V_ARR(TOT,1:ELE_IN_ROW)=(/((((X**(NN-1)))/(FACT(NN-1)))&
    &,NN=1,POLY_ORDER+1)/)
!    WRITE(*,'(*(F5.2))')V_ARR(TOT,:)
END DO

END SUBROUTINE VANDERMONDE_1D
!****************************************************************************!
FUNCTION GENERAL_POLYNOMIAL(X,Y,POLY_ORDER,ELE_IN_ROW)RESULT(arr)

INTEGER                                                 ::tot,MM,NN
INTEGER,INTENT(IN)                                      ::Y,POLY_ORDER,ELE_IN_ROW,X
REAL(KIND=rk8),DIMENSION(ELE_IN_ROW)                    ::arr

tot=0
DO MM=1,POLY_ORDER+1
DO NN=1,MM

tot=tot+1
arr(tot)=(X**(NN-1)*Y**(MM-NN))/FACT(NN-1)*FACT(MM-NN)

ENDDO
ENDDO

END FUNCTION
!****************************************************************************!
SUBROUTINE VANDERMONDE_2D_other(V_ARR,LEN)
!*****************************************************************************!
! EQUATION.......
!
!*****************************************************************************!
!NEW_EQU = (a0+a1X+a2Y)(X**2+Y**2)**0 + 
!          (b0+b1X+b2Y)(X**2+Y**2)**1 +
!          (c0+c1X+c2Y)(X**2+Y**2)**2 +
!          (d0+d1X+d2Y)(X**2+Y**2)**3 +
!          ____________ __________
!          (First term) (Second term)**(N)

!POLYNOMIAL ORDER of NEW_EQU:: First term + second term = highest POLY order
!                                       1 +    N*2      = hugest POLY order

!First term is linear and second term is addition of terms.
!N is a power term of second term which will provide the polynomial order.
!
! NUMBER OF ELEMENTS IN N TH ORDER POLYNOMIAL
! TOT no of ELE/VARIABLES in a equ is = 3 (for first term) * 2**N (for second term)
! 
!*****************************************************************************!
INTEGER                                    ::X,Y
INTEGER                                    ::I,J,II,JJ
INTEGER                                    ::LL,KK
INTEGER                                    ::N,NN,MM
INTEGER                                    ::POLY_ORDER
INTEGER                                    ::TOT
INTEGER                                    ::ELE_IN_ROW
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)    ::V_ARR
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)      ::T_ARR
INTEGER,DIMENSION(3),INTENT(INOUT)         ::LEN
INTEGER,DIMENSION(3)                       ::CTR
INTEGER,DIMENSION(:,:),allocatable           ::aa
!=====================
! NUMBER OF ELEMENTS IN N TH ORDER POLYNOMIAL
! 3 is for first term * 2**power is for second term
!=====================
CTR=(LEN-1)/2

TOT=0

DO J=-CTR(2),CTR(2)
DO I=-CTR(1),CTR(1)
    TOT=TOT+1
END DO
END DO


!=====================
N=2
jj=N
!=====================
POLY_ORDER=1+N*2

ELE_IN_ROW=0
DO j=N,1,-1
ELE_IN_ROW=j*3+ELE_IN_ROW
ENDDO
        IF(ALLOCATED(V_ARR))DEALLOCATE(V_ARR)
        ALLOCATE(V_ARR(tot,1:ELE_IN_ROW))


TOT=0
DO Y=-CTR(2),CTR(2)
DO X=-CTR(1),CTR(1)
!X=0;Y=2

    TOT=TOT+1

LL=0
DO N=1,jj

        !ELE_IN_ROW=3*(2**N)
        ELE_IN_ROW=N*3

        IF(ALLOCATED(T_ARR))DEALLOCATE(T_ARR)
        ALLOCATE(T_ARR(1:ELE_IN_ROW))
        T_arr=0.0
        IF(ALLOCATED(AA))DEALLOCATE(AA)
        ALLOCATE(aa(0:N-1,0:N-1))
        aa=0
!=========================================================
!Pascal's triage formation and its coefficients are 
!stored in as columns in j. 
!        DO i=0,N-1
!        DO j=0,i
!            aa(i,0)=1
!            if(j .ne. 0)then
!                aa(i,j)=aa(i-1,j)+aa(i-1,j-1)
!            endif
!        ENDDO
!        ENDDO
!=========================================================
        DO i=0,N-1
        DO j=0,i
            aa(i,0)=1
            IF(j .ne. 0)then
                aa(i,j)=aa(i-1,j)+aa(i-1,j-1)
            ENDIF
            IF(i==N-1)then    
                MM=(((j+1)*3)-(2))
                NN=MM+2
              T_arr(MM:NN)=(/(1)*(aa(i,j))*(X**(2*(i-j)))*(Y**(2*j)),&
                                            & (X)*(aa(i,j))*(X**(2*(i-j)))*(Y**(2*j)),&
                                            & (Y)*(aa(i,j))*(X**(2*(i-j)))*(Y**(2*j))/)
            ENDIF
        ENDDO
        ENDDO
!=========================================================
        KK=LL+1
        LL=KK+ELE_IN_ROW-1

!if(X==1 .and. Y==0)write(*,*)X,Y,'T_arr, to be added===========',int(T_arr(1:ELE_IN_ROW),2)
    V_ARR(TOT,KK:LL)=T_arr(1:ELE_IN_ROW)
END DO
write(*,*)X,Y,'|',LL,'|',int(v_arr(tot,1:LL),2)
END DO
END DO

END SUBROUTINE

END MODULE VANDERMONDE_MODULE
