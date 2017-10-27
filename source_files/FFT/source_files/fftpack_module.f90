module fftpack_module
use constants_module
implicit none
      ! real             replaced by real(kind=rk)
      ! double precision replaced by real(kind=dk)
      ! all constants replaced by 1.234_rk or 1.234_dk notation
      ! expressions involving atan(1.) replaced by atan(1._dk)
      ! float(.) replaced by real(.,kind=rk)
      ! all variables have been declared now; forced by implicit none

      public


      integer,parameter   :: debug=0

integer,parameter        :: CFFT1_sign_type=1
integer,parameter        :: RFFT1_sign_type=2
integer,parameter        :: COST1_sign_type=3
integer,parameter        :: SINT1_sign_type=4
integer,parameter        :: COSQ1_sign_type=5
integer,parameter        :: SINQ1_sign_type=6

type fft_real_sign_type
  integer                                    :: switch=0  !=0 not yet defined,1 for MCFTI1,2 for RFFT1I
  integer                                    :: sign_type
  integer                                    :: N
  integer                                    :: NF
  integer      ,allocatable,dimension(:)     :: FAC
  real(kind=dk),allocatable,dimension(:)     :: WA
  real(kind=dk),allocatable,dimension(:)     :: WSAVE
  real(kind=dk),allocatable,dimension(:)     :: WORK
  real(kind=dk),allocatable,dimension(:)     :: WORK_a
end type fft_real_sign_type

complex(kind=dk),parameter :: root3=exp(-2._dk*PI_long/3._dk*(0._dk,1._dk))
complex(kind=dk),parameter :: root51=exp(-2._dk*PI_long/5._dk*(0._dk,1._dk))
complex(kind=dk),parameter :: root52=exp(-4._dk*PI_long/5._dk*(0._dk,1._dk))

interface factor
   module procedure factor_integer
   module procedure factor_real
end interface


contains
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F2KB (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CC(IN1,L1,IDO,2),CH(IN2,L1,2,IDO),WA(IDO,1,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
!
      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      DO 101 K=1,L1
         CHOLD1 = CC(1,K,1,1)+CC(1,K,1,2)
         CC(1,K,1,2) = CC(1,K,1,1)-CC(1,K,1,2)
         CC(1,K,1,1) = CHOLD1
         CHOLD2 = CC(2,K,1,1)+CC(2,K,1,2)
         CC(2,K,1,2) = CC(2,K,1,1)-CC(2,K,1,2)
         CC(2,K,1,1) = CHOLD2
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         CH(1,K,1,1) = CC(1,K,1,1)+CC(1,K,1,2)
         CH(1,K,2,1) = CC(1,K,1,1)-CC(1,K,1,2)
         CH(2,K,1,1) = CC(2,K,1,1)+CC(2,K,1,2)
         CH(2,K,2,1) = CC(2,K,1,1)-CC(2,K,1,2)
  103 CONTINUE
      IF(IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
         DO 104 K=1,L1
            CH(1,K,1,I) = CC(1,K,I,1)+CC(1,K,I,2)
            TR2 = CC(1,K,I,1)-CC(1,K,I,2)
            CH(2,K,1,I) = CC(2,K,I,1)+CC(2,K,I,2)
            TI2 = CC(2,K,I,1)-CC(2,K,I,2)
            CH(2,K,2,I) = WA(I,1,1)*TI2+WA(I,1,2)*TR2
            CH(1,K,2,I) = WA(I,1,1)*TR2-WA(I,1,2)*TI2
  104    CONTINUE
  105 CONTINUE
      END subroutine C1F2KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F2KF (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: SN
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CC(IN1,L1,IDO,2),CH(IN2,L1,2,IDO),WA(IDO,1,2)
!
      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(2*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         CHOLD1 = SN*(CC(1,K,1,1)+CC(1,K,1,2))
         CC(1,K,1,2) = SN*(CC(1,K,1,1)-CC(1,K,1,2))
         CC(1,K,1,1) = CHOLD1
         CHOLD2 = SN*(CC(2,K,1,1)+CC(2,K,1,2))
         CC(2,K,1,2) = SN*(CC(2,K,1,1)-CC(2,K,1,2))
         CC(2,K,1,1) = CHOLD2

         ! complex
       !!  CHOLD1 = SN*(CC(K,1,1)+CC(K,1,2))
       !!  CC(K,1,2) = SN*(CC(K,1,1)-CC(K,1,2))
       !!  CC(K,1,1) = CHOLD1
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         CH(1,K,1,1) = SN*(CC(1,K,1,1)+CC(1,K,1,2))
         CH(1,K,2,1) = SN*(CC(1,K,1,1)-CC(1,K,1,2))
         CH(2,K,1,1) = SN*(CC(2,K,1,1)+CC(2,K,1,2))
         CH(2,K,2,1) = SN*(CC(2,K,1,1)-CC(2,K,1,2))

         ! complex
       !!  CH(K,1,1) = SN*(CC(K,1,1)+CC(K,1,2))
       !!  CH(K,2,1) = SN*(CC(K,1,1)-CC(K,1,2))
  107 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         CH(1,K,1,1) = CC(1,K,1,1)+CC(1,K,1,2)
         CH(1,K,2,1) = CC(1,K,1,1)-CC(1,K,1,2)
         CH(2,K,1,1) = CC(2,K,1,1)+CC(2,K,1,2)
         CH(2,K,2,1) = CC(2,K,1,1)-CC(2,K,1,2)

         ! complex
         !!CH(K,1,1) = CC(K,1,1)+CC(K,1,2)
         !!CH(K,2,1) = CC(K,1,1)-CC(K,1,2)
  103 CONTINUE
      DO 105 I=2,IDO
         DO 104 K=1,L1
            CH(1,K,1,I) = CC(1,K,I,1)+CC(1,K,I,2)
            TR2 = CC(1,K,I,1)-CC(1,K,I,2)
            CH(2,K,1,I) = CC(2,K,I,1)+CC(2,K,I,2)
            TI2 = CC(2,K,I,1)-CC(2,K,I,2)
            CH(2,K,2,I) = WA(I,1,1)*TI2-WA(I,1,2)*TR2
            CH(1,K,2,I) = WA(I,1,1)*TR2+WA(I,1,2)*TI2

         ! complex
         !!   CH(K,1,I) = CC(K,I,1)+CC(K,I,2)
         !!   TT = CC(K,I,1)-CC(K,I,2)
         !!   CH(1,K,2,I) = WA(I,1,1)*TR2-(-WA(I,1,2))*TI2
         !!   CH(2,K,2,I) = WA(I,1,1)*TI2+(-WA(I,1,2))*TR2
         !!   re x*y^  =  re x * re y + im x im y
         !!   im x*y^  = -re x * im y + im x re y
         !!   CH(K,2,I) = TT*conjg(WA(I,1))
  104    CONTINUE
  105 CONTINUE
      END subroutine C1F2KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F3KB (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(IN1,L1,IDO,3),CH(IN2,L1,3,IDO),WA(IDO,2,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TAUR,TAUI
      !  (TAUR,TAUI)=exp(+2*PI/3*i)
      DATA TAUR,TAUI /-.5_rk,.866025403784439_rk/
!
      if(debug>0) write(*,*) 'C1F3KB'  !test_print
      TAUR=real(root3,kind=rk);TAUI=-aimag(root3)
      !!write(*,'(*(a,e11.3))') 'TAUR ',TAUR,' TAUI ',TAUI

      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      write(*,*) 'C1F3KB 101'  !test_print
      DO 101 K=1,L1
         TR2 = CC(1,K,1,2)+CC(1,K,1,3)
         CR2 = CC(1,K,1,1)+TAUR*TR2
         CC(1,K,1,1) = CC(1,K,1,1)+TR2
         TI2 = CC(2,K,1,2)+CC(2,K,1,3)
         CI2 = CC(2,K,1,1)+TAUR*TI2
         CC(2,K,1,1) = CC(2,K,1,1)+TI2
         CR3 = TAUI*(CC(1,K,1,2)-CC(1,K,1,3))
         CI3 = TAUI*(CC(2,K,1,2)-CC(2,K,1,3))
         CC(1,K,1,2) = CR2-CI3
         CC(1,K,1,3) = CR2+CI3
         CC(2,K,1,2) = CI2+CR3
         CC(2,K,1,3) = CI2-CR3
     !!   write(*,'(2i4,3(2e10.3))') k,i,CH(1,K,1,1),CH(2,K,1,1),CH(1,K,1,2),CH(2,K,1,2),CH(1,K,1,3),CH(2,K,1,3)
  101 CONTINUE
      RETURN
  102 continue
      !test! write(*,*) 'C1F3KB 103'  !test_print
      DO 103 K=1,L1
         TR2 = CC(1,K,1,2)+CC(1,K,1,3)
         CR2 = CC(1,K,1,1)+TAUR*TR2
         CH(1,K,1,1) = CC(1,K,1,1)+TR2
         TI2 = CC(2,K,1,2)+CC(2,K,1,3)
         CI2 = CC(2,K,1,1)+TAUR*TI2
         CH(2,K,1,1) = CC(2,K,1,1)+TI2
         CR3 = TAUI*(CC(1,K,1,2)-CC(1,K,1,3))
         CI3 = TAUI*(CC(2,K,1,2)-CC(2,K,1,3))
         CH(1,K,2,1) = CR2-CI3
         CH(1,K,3,1) = CR2+CI3
         CH(2,K,2,1) = CI2+CR3
         CH(2,K,3,1) = CI2-CR3
  !!      write(*,*)
  !!      write(*,'(2i4,3(2e10.3))') k,i,CC(1,K,1,1),CC(2,K,1,1),CC(1,K,1,2),CC(2,K,1,2),CC(1,K,1,3),CC(2,K,1,3)
  !!      write(*,'(2i4,3(2e10.3))') k,i,TR2,TI2,CR2,CI2,CR3,CI3 
  !!      write(*,'(2i4,3(2e10.3))') k,i,CH(1,K,1,1),CH(2,K,1,1),CH(1,K,2,1),CH(2,K,2,1),CH(1,K,3,1),CH(2,K,3,1)
  103 CONTINUE
      IF (IDO .EQ. 1) RETURN
      !test! write(*,*) 'C1F3KB 105'  !test_print
      DO 105 I=2,IDO
        DO 104 K=1,L1
            TR2 = CC(1,K,I,2)+CC(1,K,I,3)
            CR2 = CC(1,K,I,1)+TAUR*TR2
            CH(1,K,1,I) = CC(1,K,I,1)+TR2
            TI2 = CC(2,K,I,2)+CC(2,K,I,3)
            CI2 = CC(2,K,I,1)+TAUR*TI2
            CH(2,K,1,I) = CC(2,K,I,1)+TI2
            CR3 = TAUI*(CC(1,K,I,2)-CC(1,K,I,3))
            CI3 = TAUI*(CC(2,K,I,2)-CC(2,K,I,3))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(1,K,2,I) = WA(I,1,1)*DR2-WA(I,1,2)*DI2
            CH(2,K,2,I) = WA(I,1,1)*DI2+WA(I,1,2)*DR2
            CH(1,K,3,I) = WA(I,2,1)*DR3-WA(I,2,2)*DI3
            CH(2,K,3,I) = WA(I,2,1)*DI3+WA(I,2,2)*DR3
    !!    write(*,*)
    !!    write(*,'(2i4,3(2e11.4))') k,i,TR2,TI2,CR2,CI2,CR3,CI3
    !!    write(*,'(2i4,3(2e11.4))') k,i,DR2,DI2,DR3,DI3,DR2-DR3,DI2-DI3
    !!    write(*,'(2i4,3(2e11.4))') k,i,CH(1,K,1,I),CH(2,K,1,I),CH(1,K,2,I),CH(2,K,2,I),CH(1,K,3,I),CH(2,K,3,I)
  104    CONTINUE
  105 CONTINUE
  if( debug > 0 ) write(*,*) 'C1F3KB end'
!!      DO I=1,IDO
!!        DO K=1,L1
!!        write(*,'(2i4,3(2e10.3))') k,i,CH(1,K,1,I),CH(2,K,1,I),CH(1,K,2,I),CH(2,K,2,I),CH(1,K,3,I),CH(2,K,3,I)
!!        enddo
!!      enddo
      END subroutine C1F3KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F3KF (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(IN1,L1,IDO,3),CH(IN2,L1,3,IDO),WA(IDO,2,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TAUR,TAUI
      real(kind=rk)     :: SN
      !  (TAUR,TAUI)=exp(-2*PI/3*i)
      DATA TAUR,TAUI /-.5_rk,-.866025403784439_rk/
!
      if(debug>0) write(*,*) 'C1F3KF' !test_print
      TAUR=real(root3,kind=rk);TAUI=aimag(root3)
      !!write(*,'(*(a,e11.3))') 'TAUR ',TAUR,' TAUI ',TAUI

      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(3*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         TR2 = CC(1,K,1,2)+CC(1,K,1,3)
         CR2 = CC(1,K,1,1)+TAUR*TR2
         CC(1,K,1,1) = SN*(CC(1,K,1,1)+TR2)
         TI2 = CC(2,K,1,2)+CC(2,K,1,3)
         CI2 = CC(2,K,1,1)+TAUR*TI2
         CC(2,K,1,1) = SN*(CC(2,K,1,1)+TI2)
         CR3 = TAUI*(CC(1,K,1,2)-CC(1,K,1,3))
         CI3 = TAUI*(CC(2,K,1,2)-CC(2,K,1,3))
         CC(1,K,1,2) = SN*(CR2-CI3)
         CC(1,K,1,3) = SN*(CR2+CI3)
         CC(2,K,1,2) = SN*(CI2+CR3)
         CC(2,K,1,3) = SN*(CI2-CR3)
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         TR2 = CC(1,K,1,2)+CC(1,K,1,3)
         CR2 = CC(1,K,1,1)+TAUR*TR2
         CH(1,K,1,1) = SN*(CC(1,K,1,1)+TR2)
         TI2 = CC(2,K,1,2)+CC(2,K,1,3)
         CI2 = CC(2,K,1,1)+TAUR*TI2
         CH(2,K,1,1) = SN*(CC(2,K,1,1)+TI2)
         CR3 = TAUI*(CC(1,K,1,2)-CC(1,K,1,3))
         CI3 = TAUI*(CC(2,K,1,2)-CC(2,K,1,3))
         CH(1,K,2,1) = SN*(CR2-CI3)
         CH(1,K,3,1) = SN*(CR2+CI3)
         CH(2,K,2,1) = SN*(CI2+CR3)
         CH(2,K,3,1) = SN*(CI2-CR3)
  107 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         TR2 = CC(1,K,1,2)+CC(1,K,1,3)
         CR2 = CC(1,K,1,1)+TAUR*TR2
         CH(1,K,1,1) = CC(1,K,1,1)+TR2
         TI2 = CC(2,K,1,2)+CC(2,K,1,3)
         CI2 = CC(2,K,1,1)+TAUR*TI2
         CH(2,K,1,1) = CC(2,K,1,1)+TI2
         CR3 = TAUI*(CC(1,K,1,2)-CC(1,K,1,3))
         CI3 = TAUI*(CC(2,K,1,2)-CC(2,K,1,3))
         CH(1,K,2,1) = CR2-CI3
         CH(1,K,3,1) = CR2+CI3
         CH(2,K,2,1) = CI2+CR3
         CH(2,K,3,1) = CI2-CR3
  103 CONTINUE
      DO 105 I=2,IDO
        DO 104 K=1,L1
            TR2 = CC(1,K,I,2)+CC(1,K,I,3)
            CR2 = CC(1,K,I,1)+TAUR*TR2
            CH(1,K,1,I) = CC(1,K,I,1)+TR2
            TI2 = CC(2,K,I,2)+CC(2,K,I,3)
            CI2 = CC(2,K,I,1)+TAUR*TI2
            CH(2,K,1,I) = CC(2,K,I,1)+TI2
            CR3 = TAUI*(CC(1,K,I,2)-CC(1,K,I,3))
            CI3 = TAUI*(CC(2,K,I,2)-CC(2,K,I,3))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(2,K,2,I) = WA(I,1,1)*DI2-WA(I,1,2)*DR2
            CH(1,K,2,I) = WA(I,1,1)*DR2+WA(I,1,2)*DI2
            CH(2,K,3,I) = WA(I,2,1)*DI3-WA(I,2,2)*DR3
            CH(1,K,3,I) = WA(I,2,1)*DR3+WA(I,2,2)*DI3
  104    CONTINUE
  105 CONTINUE
      END subroutine C1F3KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F4KB (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      REAL(kind=rk)     :: CC(IN1,L1,IDO,4),CH(IN2,L1,4,IDO),WA(IDO,3,2)
!
! FFTPACK 5.1 auxiliary routine
!
      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      print*,'in 101 von C1F4KB'
      DO 101 K=1,L1
      !!stop 'in C1F4KB'
         TI1 = CC(2,K,1,1)-CC(2,K,1,3)
         TI2 = CC(2,K,1,1)+CC(2,K,1,3)
         TI3 = CC(2,K,1,2)+CC(2,K,1,4)
         TR1 = CC(1,K,1,1)-CC(1,K,1,3)
         TR2 = CC(1,K,1,1)+CC(1,K,1,3)
         TR3 = CC(1,K,1,2)+CC(1,K,1,4)

         TR4 = CC(2,K,1,4)-CC(2,K,1,2)
         TI4 = CC(1,K,1,2)-CC(1,K,1,4)

         !!TR4 = CC(1,K,1,4)-CC(1,K,1,2)   ! auf jeden Fall falsch
         !!TI4 = CC(2,K,1,4)-CC(2,K,1,2)   ! auf jeden Fall falsch
         !!TR4 = CC(1,K,1,2)-CC(1,K,1,4)   ! auch falsch
         !!TI4 = CC(2,K,1,2)-CC(2,K,1,4)   ! auch falsch

         CC(1,K,1,1) = TR2+TR3
         CC(1,K,1,3) = TR2-TR3
         CC(2,K,1,1) = TI2+TI3
         CC(2,K,1,3) = TI2-TI3
         CC(1,K,1,2) = TR1+TR4
         CC(1,K,1,4) = TR1-TR4
         CC(2,K,1,2) = TI1+TI4
         CC(2,K,1,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         TI1 = CC(2,K,1,1)-CC(2,K,1,3)
         TI2 = CC(2,K,1,1)+CC(2,K,1,3)
         TR4 = CC(2,K,1,4)-CC(2,K,1,2)
         TI3 = CC(2,K,1,2)+CC(2,K,1,4)
         TR1 = CC(1,K,1,1)-CC(1,K,1,3)
         TR2 = CC(1,K,1,1)+CC(1,K,1,3)
         TI4 = CC(1,K,1,2)-CC(1,K,1,4)
         TR3 = CC(1,K,1,2)+CC(1,K,1,4)
         CH(1,K,1,1) = TR2+TR3
         CH(1,K,3,1) = TR2-TR3
         CH(2,K,1,1) = TI2+TI3
         CH(2,K,3,1) = TI2-TI3
         CH(1,K,2,1) = TR1+TR4
         CH(1,K,4,1) = TR1-TR4
         CH(2,K,2,1) = TI1+TI4
         CH(2,K,4,1) = TI1-TI4
  103 CONTINUE
      IF(IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
         DO 104 K=1,L1
            TI1 = CC(2,K,I,1)-CC(2,K,I,3)
            TI2 = CC(2,K,I,1)+CC(2,K,I,3)
            TI3 = CC(2,K,I,2)+CC(2,K,I,4)
            TR4 = CC(2,K,I,4)-CC(2,K,I,2)
            TR1 = CC(1,K,I,1)-CC(1,K,I,3)
            TR2 = CC(1,K,I,1)+CC(1,K,I,3)
            TI4 = CC(1,K,I,2)-CC(1,K,I,4)
            TR3 = CC(1,K,I,2)+CC(1,K,I,4)
            CH(1,K,1,I) = TR2+TR3
            CR3 = TR2-TR3
            CH(2,K,1,I) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(1,K,2,I) = WA(I,1,1)*CR2-WA(I,1,2)*CI2
            CH(2,K,2,I) = WA(I,1,1)*CI2+WA(I,1,2)*CR2
            CH(1,K,3,I) = WA(I,2,1)*CR3-WA(I,2,2)*CI3
            CH(2,K,3,I) = WA(I,2,1)*CI3+WA(I,2,2)*CR3
            CH(1,K,4,I) = WA(I,3,1)*CR4-WA(I,3,2)*CI4
            CH(2,K,4,I) = WA(I,3,1)*CI4+WA(I,3,2)*CR4
  104    CONTINUE
  105 CONTINUE
      END subroutine C1F4KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F4KF (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: SN
      REAL(kind=rk)     :: CC(IN1,L1,IDO,4),CH(IN2,L1,4,IDO),WA(IDO,3,2)
!
! FFTPACK 5.1 auxiliary routine
!
      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(4*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         TI1 = CC(2,K,1,1)-CC(2,K,1,3)
         TI2 = CC(2,K,1,1)+CC(2,K,1,3)
         TR4 = CC(2,K,1,2)-CC(2,K,1,4)
         TI3 = CC(2,K,1,2)+CC(2,K,1,4)
         TR1 = CC(1,K,1,1)-CC(1,K,1,3)
         TR2 = CC(1,K,1,1)+CC(1,K,1,3)
         TI4 = CC(1,K,1,4)-CC(1,K,1,2)
         TR3 = CC(1,K,1,2)+CC(1,K,1,4)
         CC(1,K,1,1) = SN*(TR2+TR3)
         CC(1,K,1,3) = SN*(TR2-TR3)
         CC(2,K,1,1) = SN*(TI2+TI3)
         CC(2,K,1,3) = SN*(TI2-TI3)
         CC(1,K,1,2) = SN*(TR1+TR4)
         CC(1,K,1,4) = SN*(TR1-TR4)
         CC(2,K,1,2) = SN*(TI1+TI4)
         CC(2,K,1,4) = SN*(TI1-TI4)
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         TI1 = CC(2,K,1,1)-CC(2,K,1,3)
         TI2 = CC(2,K,1,1)+CC(2,K,1,3)
         TR4 = CC(2,K,1,2)-CC(2,K,1,4)
         TI3 = CC(2,K,1,2)+CC(2,K,1,4)
         TR1 = CC(1,K,1,1)-CC(1,K,1,3)
         TR2 = CC(1,K,1,1)+CC(1,K,1,3)
         TI4 = CC(1,K,1,4)-CC(1,K,1,2)
         TR3 = CC(1,K,1,2)+CC(1,K,1,4)
         CH(1,K,1,1) = SN*(TR2+TR3)
         CH(1,K,3,1) = SN*(TR2-TR3)
         CH(2,K,1,1) = SN*(TI2+TI3)
         CH(2,K,3,1) = SN*(TI2-TI3)
         CH(1,K,2,1) = SN*(TR1+TR4)
         CH(1,K,4,1) = SN*(TR1-TR4)
         CH(2,K,2,1) = SN*(TI1+TI4)
         CH(2,K,4,1) = SN*(TI1-TI4)
  107 CONTINUE
      RETURN
  102 continue
!test!print*,'C1F4KF 103'  !test_print
      DO 103 K=1,L1
         TI1 = CC(2,K,1,1)-CC(2,K,1,3)
         TI2 = CC(2,K,1,1)+CC(2,K,1,3)
         TR4 = CC(2,K,1,2)-CC(2,K,1,4)
         TI3 = CC(2,K,1,2)+CC(2,K,1,4)
         TR1 = CC(1,K,1,1)-CC(1,K,1,3)
         TR2 = CC(1,K,1,1)+CC(1,K,1,3)
         TI4 = CC(1,K,1,4)-CC(1,K,1,2)
         TR3 = CC(1,K,1,2)+CC(1,K,1,4)
         CH(1,K,1,1) = TR2+TR3
         CH(1,K,3,1) = TR2-TR3
         CH(2,K,1,1) = TI2+TI3
         CH(2,K,3,1) = TI2-TI3
         CH(1,K,2,1) = TR1+TR4
         CH(1,K,4,1) = TR1-TR4
         CH(2,K,2,1) = TI1+TI4
         CH(2,K,4,1) = TI1-TI4
         !!print*,k,CH(1,K,1,1),CH(2,K,1,1),CH(1,K,2,1),CH(2,K,2,1)
  103 CONTINUE
!test!print*,'C1F4KF 105'  !test_print
      DO 105 I=2,IDO
         DO 104 K=1,L1
!!         write(*,*)
!!         write(*,'(i4,i4,4(1x,2e10.3))') k,i,CC(1,K,I,1),CC(2,K,I,1)  &
!!                                           &,CC(1,K,I,2),CC(2,K,I,2)  &
!!                                           &,CC(1,K,I,3),CC(2,K,I,3)  &
!!                                           &,CC(1,K,I,4),CC(2,K,I,4)
            TI1 = CC(2,K,I,1)-CC(2,K,I,3)
            TI2 = CC(2,K,I,1)+CC(2,K,I,3)
            TI3 = CC(2,K,I,2)+CC(2,K,I,4)
            TR4 = CC(2,K,I,2)-CC(2,K,I,4)
            TR1 = CC(1,K,I,1)-CC(1,K,I,3)
            TR2 = CC(1,K,I,1)+CC(1,K,I,3)
            TI4 = CC(1,K,I,4)-CC(1,K,I,2)
            TR3 = CC(1,K,I,2)+CC(1,K,I,4)
            CH(1,K,1,I) = TR2+TR3
            CR3 = TR2-TR3
            CH(2,K,1,I) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(1,K,2,I) = WA(I,1,1)*CR2+WA(I,1,2)*CI2
            CH(2,K,2,I) = WA(I,1,1)*CI2-WA(I,1,2)*CR2
            CH(1,K,3,I) = WA(I,2,1)*CR3+WA(I,2,2)*CI3
            CH(2,K,3,I) = WA(I,2,1)*CI3-WA(I,2,2)*CR3
            CH(1,K,4,I) = WA(I,3,1)*CR4+WA(I,3,2)*CI4
            CH(2,K,4,I) = WA(I,3,1)*CI4-WA(I,3,2)*CR4
!!         write(*,'(i4,i4,4(1x,2e10.3))') k,i,TR1,TI1  &
!!                                           &,TR2,TI2  &
!!                                           &,TR3,TI3  &
!!                                           &,TR4,TI4  
!!         write(*,'(i4,i4,4(1x,2e10.3))') k,i,CH(1,K,1,I),CH(2,K,1,I) &
!!                                           &,CH(1,K,2,I),CH(2,K,2,I) &
!!                                           &,CH(1,K,3,I),CH(2,K,3,I) &
!!                                           &,CH(1,K,4,I),CH(2,K,4,I)
  104    CONTINUE
  105 CONTINUE
!!      DO I=1,IDO
!!         DO K=1,L1
!!         write(*,'(i4,i4,4(1x,2e10.3))') k,i,CH(1,K,1,I),CH(2,K,1,I) &
!!                                           &,CH(1,K,2,I),CH(2,K,2,I) &
!!                                           &,CH(1,K,3,I),CH(2,K,3,I) &
!!                                           &,CH(1,K,4,I),CH(2,K,4,I)
!!         enddo
!!     enddo
      END subroutine C1F4KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F5KB (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(IN1,L1,IDO,5),CH(IN2,L1,5,IDO),WA(IDO,4,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      real(kind=rk)     :: CHOLD1,CHOLD2
      !   (TR11,TI11)=exp(+2*PI/5*i), (TR12,TI12)=exp(+4*PI/5*i)
      DATA TR11,TI11,TR12,TI12 /.3090169943749474_rk,.9510565162951536_rk,-.8090169943749474_rk,.5877852522924731_rk/
!
! FFTPACK 5.1 auxiliary routine
!
!!      write(*,*) 'C1F5KB'
      TR11=real(root51,kind=rk);TI11=-aimag(root51)
      TR12=real(root52,kind=rk);TI12=-aimag(root52)

      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,K,1,2)-CC(2,K,1,5)
         TI2 = CC(2,K,1,2)+CC(2,K,1,5)
         TI4 = CC(2,K,1,3)-CC(2,K,1,4)
         TI3 = CC(2,K,1,3)+CC(2,K,1,4)
         TR5 = CC(1,K,1,2)-CC(1,K,1,5)
         TR2 = CC(1,K,1,2)+CC(1,K,1,5)
         TR4 = CC(1,K,1,3)-CC(1,K,1,4)
         TR3 = CC(1,K,1,3)+CC(1,K,1,4)
         CHOLD1 = CC(1,K,1,1)+TR2+TR3
         CHOLD2 = CC(2,K,1,1)+TI2+TI3
         CR2 = CC(1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,K,1,1)+TR12*TI2+TR11*TI3
         CC(1,K,1,1) = CHOLD1
         CC(2,K,1,1) = CHOLD2
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CC(1,K,1,2) = CR2-CI5
         CC(1,K,1,5) = CR2+CI5
         CC(2,K,1,2) = CI2+CR5
         CC(2,K,1,3) = CI3+CR4
         CC(1,K,1,3) = CR3-CI4
         CC(1,K,1,4) = CR3+CI4
         CC(2,K,1,4) = CI3-CR4
         CC(2,K,1,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         TI5 = CC(2,K,1,2)-CC(2,K,1,5)
         TI2 = CC(2,K,1,2)+CC(2,K,1,5)
         TI4 = CC(2,K,1,3)-CC(2,K,1,4)
         TI3 = CC(2,K,1,3)+CC(2,K,1,4)
         TR5 = CC(1,K,1,2)-CC(1,K,1,5)
         TR2 = CC(1,K,1,2)+CC(1,K,1,5)
         TR4 = CC(1,K,1,3)-CC(1,K,1,4)
         TR3 = CC(1,K,1,3)+CC(1,K,1,4)
         CH(1,K,1,1) = CC(1,K,1,1)+TR2+TR3
         CH(2,K,1,1) = CC(2,K,1,1)+TI2+TI3
         CR2 = CC(1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,K,1,1)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2,1) = CR2-CI5
         CH(1,K,5,1) = CR2+CI5
         CH(2,K,2,1) = CI2+CR5
         CH(2,K,3,1) = CI3+CR4
         CH(1,K,3,1) = CR3-CI4
         CH(1,K,4,1) = CR3+CI4
         CH(2,K,4,1) = CI3-CR4
         CH(2,K,5,1) = CI2-CR5
  103 CONTINUE
      IF(IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
         DO 104 K=1,L1
            TI5 = CC(2,K,I,2)-CC(2,K,I,5)
            TI2 = CC(2,K,I,2)+CC(2,K,I,5)
            TI4 = CC(2,K,I,3)-CC(2,K,I,4)
            TI3 = CC(2,K,I,3)+CC(2,K,I,4)
            TR5 = CC(1,K,I,2)-CC(1,K,I,5)
            TR2 = CC(1,K,I,2)+CC(1,K,I,5)
            TR4 = CC(1,K,I,3)-CC(1,K,I,4)
            TR3 = CC(1,K,I,3)+CC(1,K,I,4)
            CH(1,K,1,I) = CC(1,K,I,1)+TR2+TR3
            CH(2,K,1,I) = CC(2,K,I,1)+TI2+TI3
            CR2 = CC(1,K,I,1)+TR11*TR2+TR12*TR3
            CI2 = CC(2,K,I,1)+TR11*TI2+TR12*TI3
            CR3 = CC(1,K,I,1)+TR12*TR2+TR11*TR3
            CI3 = CC(2,K,I,1)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(1,K,2,I) = WA(I,1,1)*DR2-WA(I,1,2)*DI2
            CH(2,K,2,I) = WA(I,1,1)*DI2+WA(I,1,2)*DR2
            CH(1,K,3,I) = WA(I,2,1)*DR3-WA(I,2,2)*DI3
            CH(2,K,3,I) = WA(I,2,1)*DI3+WA(I,2,2)*DR3
            CH(1,K,4,I) = WA(I,3,1)*DR4-WA(I,3,2)*DI4
            CH(2,K,4,I) = WA(I,3,1)*DI4+WA(I,3,2)*DR4
            CH(1,K,5,I) = WA(I,4,1)*DR5-WA(I,4,2)*DI5
            CH(2,K,5,I) = WA(I,4,1)*DI5+WA(I,4,2)*DR5
  104    CONTINUE
  105 CONTINUE
      END subroutine C1F5KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1F5KF (IDO,L1,NA,CC,IN1,CH,IN2,WA)
      integer           :: IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(IN1,L1,IDO,5),CH(IN2,L1,5,IDO),WA(IDO,4,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      real(kind=rk)     :: CHOLD1,CHOLD2
      real(kind=rk)     :: SN
      !   (TR11,TI11)=exp(-2*PI/5*i), (TR12,TI12)=exp(-4*PI/5*i)
      DATA TR11,TI11,TR12,TI12 /.3090169943749474_rk,-.9510565162951536_rk,-.8090169943749474_rk,-.5877852522924731_rk/
!
! FFTPACK 5.1 auxiliary routine
!
!!      write(*,*) 'C1F5KB IDO,L1,NA,IN1,IN2=',IDO,L1,NA,IN1,IN2
     !! write(*,*) 'T1=',cmplx(TR11,TI11),' T2=',cmplx(TR12,TI12)
     !! write(*,*) 'root51=',root51,'root52=',root52

      TR11=real(root51,kind=rk);TI11=aimag(root51)
      TR12=real(root52,kind=rk);TI12=aimag(root52)

      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(5*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
!test!print*,'C1F5KF 101'  !test_print
      DO 101 K=1,L1
         TI5 = CC(2,K,1,2)-CC(2,K,1,5)
         TI2 = CC(2,K,1,2)+CC(2,K,1,5)
         TI4 = CC(2,K,1,3)-CC(2,K,1,4)
         TI3 = CC(2,K,1,3)+CC(2,K,1,4)
         TR5 = CC(1,K,1,2)-CC(1,K,1,5)
         TR2 = CC(1,K,1,2)+CC(1,K,1,5)
         TR4 = CC(1,K,1,3)-CC(1,K,1,4)
         TR3 = CC(1,K,1,3)+CC(1,K,1,4)
         CHOLD1 = SN*(CC(1,K,1,1)+TR2+TR3)
         CHOLD2 = SN*(CC(2,K,1,1)+TI2+TI3)
         CR2 = CC(1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,K,1,1)+TR12*TI2+TR11*TI3
         CC(1,K,1,1) = CHOLD1
         CC(2,K,1,1) = CHOLD2
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CC(1,K,1,2) = SN*(CR2-CI5)
         CC(1,K,1,5) = SN*(CR2+CI5)
         CC(2,K,1,2) = SN*(CI2+CR5)
         CC(2,K,1,3) = SN*(CI3+CR4)
         CC(1,K,1,3) = SN*(CR3-CI4)
         CC(1,K,1,4) = SN*(CR3+CI4)
         CC(2,K,1,4) = SN*(CI3-CR4)
         CC(2,K,1,5) = SN*(CI2-CR5)
  101 CONTINUE
      RETURN
  106 continue
!test!print*,'C1F5KF 107'  !test_print
      DO 107 K=1,L1
         TI5 = CC(2,K,1,2)-CC(2,K,1,5)
         TI2 = CC(2,K,1,2)+CC(2,K,1,5)
         TI4 = CC(2,K,1,3)-CC(2,K,1,4)
         TI3 = CC(2,K,1,3)+CC(2,K,1,4)
         TR5 = CC(1,K,1,2)-CC(1,K,1,5)
         TR2 = CC(1,K,1,2)+CC(1,K,1,5)
         TR4 = CC(1,K,1,3)-CC(1,K,1,4)
         TR3 = CC(1,K,1,3)+CC(1,K,1,4)
         CH(1,K,1,1) = SN*(CC(1,K,1,1)+TR2+TR3)
         CH(2,K,1,1) = SN*(CC(2,K,1,1)+TI2+TI3)
         CR2 = CC(1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,K,1,1)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2,1) = SN*(CR2-CI5)
         CH(1,K,5,1) = SN*(CR2+CI5)
         CH(2,K,2,1) = SN*(CI2+CR5)
         CH(2,K,3,1) = SN*(CI3+CR4)
         CH(1,K,3,1) = SN*(CR3-CI4)
         CH(1,K,4,1) = SN*(CR3+CI4)
         CH(2,K,4,1) = SN*(CI3-CR4)
         CH(2,K,5,1) = SN*(CI2-CR5)
  107 CONTINUE
      RETURN
  102 continue
!test!print*,'C1F5KF 103'  !test_print
      DO 103 K=1,L1
         TI5 = CC(2,K,1,2)-CC(2,K,1,5)
         TI2 = CC(2,K,1,2)+CC(2,K,1,5)
         TI4 = CC(2,K,1,3)-CC(2,K,1,4)
         TI3 = CC(2,K,1,3)+CC(2,K,1,4)
         TR5 = CC(1,K,1,2)-CC(1,K,1,5)
         TR2 = CC(1,K,1,2)+CC(1,K,1,5)
         TR4 = CC(1,K,1,3)-CC(1,K,1,4)
         TR3 = CC(1,K,1,3)+CC(1,K,1,4)
         CH(1,K,1,1) = CC(1,K,1,1)+TR2+TR3
         CH(2,K,1,1) = CC(2,K,1,1)+TI2+TI3
         CR2 = CC(1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,K,1,1)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2,1) = CR2-CI5
         CH(1,K,5,1) = CR2+CI5
         CH(2,K,2,1) = CI2+CR5
         CH(2,K,3,1) = CI3+CR4
         CH(1,K,3,1) = CR3-CI4
         CH(1,K,4,1) = CR3+CI4
         CH(2,K,4,1) = CI3-CR4
         CH(2,K,5,1) = CI2-CR5
  103 CONTINUE
!test!print*,'C1F5KF 105'  !test_print
      DO 105 I=2,IDO
         DO 104 K=1,L1
            TI5 = CC(2,K,I,2)-CC(2,K,I,5)
            TI2 = CC(2,K,I,2)+CC(2,K,I,5)
            TI4 = CC(2,K,I,3)-CC(2,K,I,4)
            TI3 = CC(2,K,I,3)+CC(2,K,I,4)
            TR5 = CC(1,K,I,2)-CC(1,K,I,5)
            TR2 = CC(1,K,I,2)+CC(1,K,I,5)
            TR4 = CC(1,K,I,3)-CC(1,K,I,4)
            TR3 = CC(1,K,I,3)+CC(1,K,I,4)
            CH(1,K,1,I) = CC(1,K,I,1)+TR2+TR3
            CH(2,K,1,I) = CC(2,K,I,1)+TI2+TI3
            CR2 = CC(1,K,I,1)+TR11*TR2+TR12*TR3
            CI2 = CC(2,K,I,1)+TR11*TI2+TR12*TI3
            CR3 = CC(1,K,I,1)+TR12*TR2+TR11*TR3
            CI3 = CC(2,K,I,1)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(1,K,2,I) = WA(I,1,1)*DR2+WA(I,1,2)*DI2
            CH(2,K,2,I) = WA(I,1,1)*DI2-WA(I,1,2)*DR2
            CH(1,K,3,I) = WA(I,2,1)*DR3+WA(I,2,2)*DI3
            CH(2,K,3,I) = WA(I,2,1)*DI3-WA(I,2,2)*DR3
            CH(1,K,4,I) = WA(I,3,1)*DR4+WA(I,3,2)*DI4
            CH(2,K,4,I) = WA(I,3,1)*DI4-WA(I,3,2)*DR4
            CH(1,K,5,I) = WA(I,4,1)*DR5+WA(I,4,2)*DI5
            CH(2,K,5,I) = WA(I,4,1)*DI5-WA(I,4,2)*DR5
  104    CONTINUE
  105 CONTINUE
      END subroutine C1F5KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1FGKB (IDO,IP,L1,LID,NA,CC,CC1,IN1,CH,CH1,IN2,WA)
      integer           :: LID,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDLJ
      real(kind=rk)     :: WAR,WAI    
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CH(IN2,L1,IDO,IP),CC(IN1,L1,IP,IDO),CC1(IN1,LID,IP),CH1(IN2,LID,IP),WA(IDO,IP-1,2)
!
! FFTPACK 5.1 auxiliary routine
!
      IPP2 = IP+2
      IPPH = (IP+1)/2
      DO 110 KI=1,LID
         CH1(1,KI,1) = CC1(1,KI,1)
         CH1(2,KI,1) = CC1(2,KI,1)
  110 CONTINUE
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 112 KI=1,LID
            CH1(1,KI,J) =  CC1(1,KI,J)+CC1(1,KI,JC)
            CH1(1,KI,JC) = CC1(1,KI,J)-CC1(1,KI,JC)
            CH1(2,KI,J) =  CC1(2,KI,J)+CC1(2,KI,JC)
            CH1(2,KI,JC) = CC1(2,KI,J)-CC1(2,KI,JC)
  112    CONTINUE
  111 CONTINUE
      DO 118 J=2,IPPH
         DO 117 KI=1,LID
            CC1(1,KI,1) = CC1(1,KI,1)+CH1(1,KI,J)
            CC1(2,KI,1) = CC1(2,KI,1)+CH1(2,KI,J)
  117    CONTINUE
  118 CONTINUE
      DO 116 L=2,IPPH
         LC = IPP2-L
         DO 113 KI=1,LID
            CC1(1,KI,L) = CH1(1,KI,1)+WA(1,L-1,1)*CH1(1,KI,2)
            CC1(1,KI,LC) = WA(1,L-1,2)*CH1(1,KI,IP)
            CC1(2,KI,L) = CH1(2,KI,1)+WA(1,L-1,1)*CH1(2,KI,2)
            CC1(2,KI,LC) = WA(1,L-1,2)*CH1(2,KI,IP)
  113    CONTINUE
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = MOD((L-1)*(J-1),IP)
            WAR = WA(1,IDLJ,1)
            WAI = WA(1,IDLJ,2)
            DO 114 KI=1,LID
               CC1(1,KI,L) = CC1(1,KI,L)+WAR*CH1(1,KI,J)
               CC1(1,KI,LC) = CC1(1,KI,LC)+WAI*CH1(1,KI,JC)
               CC1(2,KI,L) = CC1(2,KI,L)+WAR*CH1(2,KI,J)
               CC1(2,KI,LC) = CC1(2,KI,LC)+WAI*CH1(2,KI,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      IF(IDO .GT. 1 .OR. NA .EQ. 1) GO TO 136
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 KI=1,LID
            CHOLD1 = CC1(1,KI,J)-CC1(2,KI,JC)
            CHOLD2 = CC1(1,KI,J)+CC1(2,KI,JC)
            CC1(1,KI,J) = CHOLD1
            CC1(2,KI,JC) = CC1(2,KI,J)-CC1(1,KI,JC)
            CC1(2,KI,J) = CC1(2,KI,J)+CC1(1,KI,JC)
            CC1(1,KI,JC) = CHOLD2
  119    CONTINUE
  120 CONTINUE
      RETURN
  136 DO 137 KI=1,LID
         CH1(1,KI,1) = CC1(1,KI,1)
         CH1(2,KI,1) = CC1(2,KI,1)
  137 CONTINUE
      DO 135 J=2,IPPH
         JC = IPP2-J
         DO 134 KI=1,LID
            CH1(1,KI,J) = CC1(1,KI,J)-CC1(2,KI,JC)
            CH1(1,KI,JC) = CC1(1,KI,J)+CC1(2,KI,JC)
            CH1(2,KI,JC) = CC1(2,KI,J)-CC1(1,KI,JC)
            CH1(2,KI,J) = CC1(2,KI,J)+CC1(1,KI,JC)
  134    CONTINUE
  135 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 131 I=1,IDO
         DO 130 K=1,L1
            CC(1,K,1,I) = CH(1,K,I,1)
            CC(2,K,1,I) = CH(2,K,I,1)
  130    CONTINUE
  131 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            CC(1,K,J,1) = CH(1,K,1,J)
            CC(2,K,J,1) = CH(2,K,1,J)
  122    CONTINUE
  123 CONTINUE
      DO 126 J=2,IP
         DO 125 I=2,IDO
            DO 124 K=1,L1
               CC(1,K,J,I) = WA(I,J-1,1)*CH(1,K,I,J)-WA(I,J-1,2)*CH(2,K,I,J)
               CC(2,K,J,I) = WA(I,J-1,1)*CH(2,K,I,J)+WA(I,J-1,2)*CH(1,K,I,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
!!     write(*,*) 'C1FGKB'
!!      DO J=1,IP
!!      DO I=1,IDO
!!        DO K=1,L1
!!        write(*,'(3i4,3(2e10.3))') j,i,k,CC(1,K,J,I),CC(2,K,J,I)
!!        enddo
!!      enddo
!!      enddo
      END subroutine C1FGKB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE C1FGKF (IDO,IP,L1,LID,NA,CC,CC1,IN1,CH,CH1,IN2,WA)
      integer           :: LID,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDLJ
      real(kind=rk)     :: WAR,WAI    
      real(kind=rk)     :: CHOLD1,CHOLD2
      real(kind=rk)     :: SN
      REAL(kind=rk)     :: CH(IN2,L1,IDO,IP),CC(IN1,L1,IP,IDO),CC1(IN1,LID,IP),CH1(IN2,LID,IP),WA(IDO,IP-1,2)
!
! FFTPACK 5.1 auxiliary routine
!
      IPP2 = IP+2
      IPPH = (IP+1)/2
      DO 110 KI=1,LID
         CH1(1,KI,1) = CC1(1,KI,1)
         CH1(2,KI,1) = CC1(2,KI,1)
  110 CONTINUE
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 112 KI=1,LID
            CH1(1,KI,J) =  CC1(1,KI,J)+CC1(1,KI,JC)
            CH1(1,KI,JC) = CC1(1,KI,J)-CC1(1,KI,JC)
            CH1(2,KI,J) =  CC1(2,KI,J)+CC1(2,KI,JC)
            CH1(2,KI,JC) = CC1(2,KI,J)-CC1(2,KI,JC)
  112    CONTINUE
  111 CONTINUE
      DO 118 J=2,IPPH
         DO 117 KI=1,LID
            CC1(1,KI,1) = CC1(1,KI,1)+CH1(1,KI,J)
            CC1(2,KI,1) = CC1(2,KI,1)+CH1(2,KI,J)
  117    CONTINUE
  118 CONTINUE
      DO 116 L=2,IPPH
         LC = IPP2-L
         DO 113 KI=1,LID
            CC1(1,KI,L) = CH1(1,KI,1)+WA(1,L-1,1)*CH1(1,KI,2)
            CC1(1,KI,LC) = -WA(1,L-1,2)*CH1(1,KI,IP)
            CC1(2,KI,L) = CH1(2,KI,1)+WA(1,L-1,1)*CH1(2,KI,2)
            CC1(2,KI,LC) = -WA(1,L-1,2)*CH1(2,KI,IP)
  113    CONTINUE
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = MOD((L-1)*(J-1),IP)
            WAR = WA(1,IDLJ,1)
            WAI = -WA(1,IDLJ,2)
            DO 114 KI=1,LID
               CC1(1,KI,L) = CC1(1,KI,L)+WAR*CH1(1,KI,J)
               CC1(1,KI,LC) = CC1(1,KI,LC)+WAI*CH1(1,KI,JC)
               CC1(2,KI,L) = CC1(2,KI,L)+WAR*CH1(2,KI,J)
               CC1(2,KI,LC) = CC1(2,KI,LC)+WAI*CH1(2,KI,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      IF (IDO .GT. 1) GO TO 136
      SN = 1._rk/REAL(IP*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 146
      DO 149 KI=1,LID
         CC1(1,KI,1) = SN*CC1(1,KI,1)
         CC1(2,KI,1) = SN*CC1(2,KI,1)
  149 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 KI=1,LID
            CHOLD1 = SN*(CC1(1,KI,J)-CC1(2,KI,JC))
            CHOLD2 = SN*(CC1(1,KI,J)+CC1(2,KI,JC))
            CC1(1,KI,J) = CHOLD1
            CC1(2,KI,JC) = SN*(CC1(2,KI,J)-CC1(1,KI,JC))
            CC1(2,KI,J) = SN*(CC1(2,KI,J)+CC1(1,KI,JC))
            CC1(1,KI,JC) = CHOLD2
  119    CONTINUE
  120 CONTINUE
      RETURN
  146 DO 147 KI=1,LID
         CH1(1,KI,1) = SN*CC1(1,KI,1)
         CH1(2,KI,1) = SN*CC1(2,KI,1)
  147 CONTINUE
      DO 145 J=2,IPPH
         JC = IPP2-J
         DO 144 KI=1,LID
            CH1(1,KI,J) = SN*(CC1(1,KI,J)-CC1(2,KI,JC))
            CH1(2,KI,J) = SN*(CC1(2,KI,J)+CC1(1,KI,JC))
            CH1(1,KI,JC) = SN*(CC1(1,KI,J)+CC1(2,KI,JC))
            CH1(2,KI,JC) = SN*(CC1(2,KI,J)-CC1(1,KI,JC))
  144    CONTINUE
  145 CONTINUE
      RETURN
  136 DO 137 KI=1,LID
         CH1(1,KI,1) = CC1(1,KI,1)
         CH1(2,KI,1) = CC1(2,KI,1)
  137 CONTINUE
      DO 135 J=2,IPPH
         JC = IPP2-J
         DO 134 KI=1,LID
            CH1(1,KI,J) = CC1(1,KI,J)-CC1(2,KI,JC)
            CH1(2,KI,J) = CC1(2,KI,J)+CC1(1,KI,JC)
            CH1(1,KI,JC) = CC1(1,KI,J)+CC1(2,KI,JC)
            CH1(2,KI,JC) = CC1(2,KI,J)-CC1(1,KI,JC)
  134    CONTINUE
  135 CONTINUE
      DO 131 I=1,IDO
         DO 130 K=1,L1
            CC(1,K,1,I) = CH(1,K,I,1)
            CC(2,K,1,I) = CH(2,K,I,1)
  130    CONTINUE
  131 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            CC(1,K,J,1) = CH(1,K,1,J)
            CC(2,K,J,1) = CH(2,K,1,J)
  122    CONTINUE
  123 CONTINUE
      DO 126 J=2,IP
         DO 125 I=2,IDO
            DO 124 K=1,L1
               CC(1,K,J,I) = WA(I,J-1,1)*CH(1,K,I,J)+WA(I,J-1,2)*CH(2,K,I,J)
               CC(2,K,J,I) = WA(I,J-1,1)*CH(2,K,I,J)-WA(I,J-1,2)*CH(1,K,I,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      END subroutine C1FGKF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE C1FM1B (N,INC,C,CH,WA,FNF,FAC,fft_sign)
      SUBROUTINE C1FM1B (N,INC,C,CH,fft_sign)
      integer           :: N,INC
      integer           :: INC2,LID,NBR
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: FNF
      !Kuester COMPLEX(kind=rk)   :: C(*)
      REAL(kind=rk)     :: C(*)
      REAL(kind=rk)     :: CH(*)
      !!REAL(kind=rk)     :: WA(*)
      !!REAL(kind=rk)     :: FAC(*)
      type(fft_real_sign_type) :: fft_sign
!
! FFTPACK 5.1 auxiliary routine
!
      if(N/=fft_sign%N) then
              write(*,*) 'C1FM1B: error N/=fft_sign%N' 
              write(*,*) 'C1FM1B: N=',N          
              write(*,*) 'C1FM1B: fft_sign%N=',fft_sign%N          
      endif

!test!      write(*,*) 'C1FM1B'   !test_print
      INC2 = INC+INC
      !!NF = FNF
      NF = fft_sign%NF
      NA = 0
      L1 = 1
      IW = 1
      DO 125 K1=1,NF
         !!IP = FAC(K1)
         IP = fft_sign%FAC(K1)
         L2 = IP*L1
         IDO = N/L2
         LID = L1*IDO
         NBR = 1+NA+2*MIN(IP-2,4)
        !! print*,k1,ip,L2,IDO,LID,NBR
         GO TO (52,62,53,63,54,64,55,65,56,66),NBR
   52    CALL C1F2KB (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   62    CALL C1F2KB (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   53    CALL C1F3KB (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   63    CALL C1F3KB (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   54    CALL C1F4KB (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   64    CALL C1F4KB (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   55    CALL C1F5KB (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   65    CALL C1F5KB (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   56    CALL C1FGKB (IDO,IP,L1,LID,NA,C,C,INC2,CH,CH,2,fft_sign%WA(IW))
         GO TO 120
   66    CALL C1FGKB (IDO,IP,L1,LID,NA,CH,CH,2,C,C,INC2,fft_sign%WA(IW))
  120    L1 = L2
         IW = IW+(IP-1)*(IDO+IDO)
         IF(IP .LE. 5) NA = 1-NA
  125 CONTINUE
      END subroutine C1FM1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE C1FM1F (N,INC,C,CH,WA,FNF,FAC,fft_sign)
      SUBROUTINE C1FM1F (N,INC,C,CH,fft_sign)
      integer           :: N,INC
      integer           :: INC2,LID,NBR
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: FNF
      !Kuester COMPLEX(kind=rk)   :: C(*)
      REAL(kind=rk)     :: C(*)
      REAL(kind=rk)     :: CH(*)
      !!REAL(kind=rk)     :: WA(*)
      !!REAL(kind=rk)     :: FAC(*)
      type(fft_real_sign_type) :: fft_sign
!
! FFTPACK 5.1 auxiliary routine
!
      if(N/=fft_sign%N) then
              write(*,*) 'C1FM1F: error N/=fft_sign%N' 
      endif

      INC2 = INC+INC
      !!NF = FNF
      NF = fft_sign%NF
      NA = 0
      L1 = 1
      IW = 1
      DO 125 K1=1,fft_sign%NF
         !!IP = FAC(K1)
         IP = fft_sign%FAC(K1)
         L2 = IP*L1
         IDO = N/L2
         LID = L1*IDO
         NBR = 1+NA+2*MIN(IP-2,4)
         GO TO (52,62,53,63,54,64,55,65,56,66),NBR
   52    CALL C1F2KF (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   62    CALL C1F2KF (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   53    CALL C1F3KF (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   63    CALL C1F3KF (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   54    CALL C1F4KF (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   64    CALL C1F4KF (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   55    CALL C1F5KF (IDO,L1,NA,C,INC2,CH,2,fft_sign%WA(IW))
         GO TO 120
   65    CALL C1F5KF (IDO,L1,NA,CH,2,C,INC2,fft_sign%WA(IW))
         GO TO 120
   56    CALL C1FGKF (IDO,IP,L1,LID,NA,C,C,INC2,CH,CH,2,fft_sign%WA(IW))
         GO TO 120
   66    CALL C1FGKF (IDO,IP,L1,LID,NA,CH,CH,2,C,C,INC2,fft_sign%WA(IW))
  120    L1 = L2
         IW = IW+(IP-1)*(IDO+IDO)
         IF(IP .LE. 5) NA = 1-NA
  125 CONTINUE
      END subroutine C1FM1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFT1B (N, INC, C, LENC, WSAVE, LENSAV,WORK, LENWRK,fft_sign, IER)
      !!SUBROUTINE CFFT1B (N, INC, C, LENC, WORK, LENWRK,fft_sign, IER)
      SUBROUTINE CFFT1B (N, INC, C, LENC, fft_sign, IER)
      integer                     :: N, INC, LENC, LENSAV
      !!integer                     :: LENWRK
      integer                     :: IER
      integer                     :: IW1
      !Kuester COMPLEX(kind=rk)   :: C(LENC)
      REAL(kind=rk)               :: C(2*LENC)
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      type(fft_real_sign_type)    :: fft_sign
!
      IER = 0
!
      IF (LENC .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('CFFT1B ', 4)
    !!  ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4) THEN
    !!    IER = 2
    !!    CALL XERFFT ('CFFT1B ', 6)
     !!  ELSEIF (LENWRK .LT. 2*N) THEN
      !!   IER = 3
       !!  CALL XERFFT ('CFFT1B ', 8)
      ENDIF

!
      IF (N .EQ. 1) RETURN
!
      !!IW1 = N+N+1
      !!CALL C1FM1B (N,INC,C,WORK,WSAVE,WSAVE(IW1),WSAVE(IW1+1),fft_sign)
      CALL C1FM1B (N,INC,C,fft_sign%WORK,fft_sign)
      END subroutine CFFT1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFT1F (N, INC, C, LENC, WSAVE, LENSAV,WORK, LENWRK,fft_sign, IER)
      !!SUBROUTINE CFFT1F (N, INC, C, LENC, WORK, LENWRK,fft_sign, IER)
      SUBROUTINE CFFT1F (N, INC, C, LENC, fft_sign, IER)
      integer     :: N, INC, LENC
      !!integer    ::  LENWRK
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IW1
      !Kuester COMPLEX(kind=rk)   :: C(LENC)
      REAL(kind=rk)     :: C(2*LENC)
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
      IF (LENC .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('CFFT1F ', 4)
    !!  ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4) THEN
    !!    IER = 2
    !!    CALL XERFFT ('CFFT1F ', 6)
     !! ELSEIF (LENWRK .LT. 2*N) THEN
     !!   IER = 3
     !!   CALL XERFFT ('CFFT1F ', 8)
      ENDIF

!
      IF (N .EQ. 1) RETURN
!
      IW1 = N+N+1
      !!CALL C1FM1F (N,INC,C,WORK,WSAVE,WSAVE(IW1),WSAVE(IW1+1),fft_sign)
      CALL C1FM1F (N,INC,C,fft_sign%WORK,fft_sign)
      END subroutine CFFT1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFT1I (N, WSAVE, LENSAV,fft_sign, IER)
      SUBROUTINE CFFT1I (N, fft_sign, IER)
      integer     :: N
      integer     ::  IER
      !!integer     :: LENSAV
      integer     :: LENWRK
      integer           :: IW1
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
  !!    IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4) THEN
  !!      IER = 2
  !!      CALL XERFFT ('CFFTMI ', 3)
  !!    ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!IW1 = N+N+1
      !!CALL MCFTI1 (N,WSAVE,WSAVE(IW1),WSAVE(IW1+1),fft_sign)
      CALL MCFTI1 (N,fft_sign)
      fft_sign%sign_type=CFFT1_sign_type

      LENWRK=2*N
      if(allocated(fft_sign%work) ) deallocate(fft_sign%work)
      allocate(fft_sign%work(LENWRK))
!
      END subroutine CFFT1I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFT2B (LDIM, L, M, C, WSAVE, LENSAV,WORK, LENWRK, IER)
      SUBROUTINE CFFT2B (LDIM, L, M, C,WORK, LENWRK,fft_sign, IER)
      integer     :: L, M, LDIM
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      !Kuester COMPLEX(kind=rk)   :: C(LDIM,M)
      integer           :: IER1
      integer           :: IW
      real(kind=rk)     :: C(2*LDIM,M)
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      integer           :: L_len
      integer           :: M_len
      type(fft_real_sign_type) :: fft_sign(:)
!
! Initialize error return
!
      IER = 0
      L_len= 2*L + INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4
      M_len= 2*M + INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
!
      IF (L .GT. LDIM) THEN
        IER = 5
        CALL XERFFT ('CFFT2B', -2)
        GO TO 100
  !!    ELSEIF (LENSAV .LT. L_len + M_len ) THEN
  !!      IER = 2
  !!      CALL XERFFT ('CFFT2B', 6)
  !!      GO TO 100
      ELSEIF (LENWRK .LT. 2*L*M) THEN
        IER = 3
        CALL XERFFT ('CFFT2B', 8)
        GO TO 100
      ENDIF
!
! Transform X lines of C array
      !!IW = 2*L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk)) + 3

      IW = L_len -1
      !!CALL CFFTMB(L,    1, M, LDIM, C, (L-1) + LDIM*(M-1) +1,WSAVE(IW), M_len ,WORK, 2*L*M, IER1)
      CALL CFFTMB(L,    1, M, LDIM, C, (L-1) + LDIM*(M-1) +1, WORK, 2*L*M,fft_sign(2), IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('CFFT2B',-5)
        GO TO 100
      ENDIF
!
! Transform Y lines of C array
      IW = 1
      !!CALL CFFTMB(M, LDIM, L,    1, C, (M-1)*LDIM + L       ,WSAVE(IW), L_len ,WORK, 2*M*L, IER1)
      CALL CFFTMB(M, LDIM, L,    1, C, (M-1)*LDIM + L       , WORK, 2*M*L,fft_sign(1), IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('CFFT2B',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine CFFT2B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFT2F (LDIM, L, M, C, WSAVE, LENSAV,WORK, LENWRK, IER)
      SUBROUTINE CFFT2F (LDIM, L, M, C, WORK, LENWRK,fft_sign, IER)
      integer     :: L, M, LDIM
      integer     ::  LENWRK, IER
      integer     :: LENSAV
      integer           :: IW
      !Kuester COMPLEX(kind=rk)   :: C(LDIM,M)
      integer           :: IER1
      real(kind=rk)     :: C(2*LDIM,M)
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      integer           :: L_len
      integer           :: M_len
      type(fft_real_sign_type) :: fft_sign(:)
!
! Initialize error return
!
      IER = 0
      L_len=2*L + INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4
      M_len=2*M + INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
!
      IF (L .GT. LDIM) THEN
        IER = 5
        CALL XERFFT ('CFFT2F', -2)
        GO TO 100
   !!   ELSEIF (LENSAV .LT. L_len + M_len ) THEN
   !!     IER = 2
   !!     CALL XERFFT ('CFFT2F', 6)
   !!     GO TO 100
      ELSEIF (LENWRK .LT. 2*L*M) THEN
        IER = 3
        CALL XERFFT ('CFFT2F', 8)
        GO TO 100
      ENDIF
!
!!     write(*,*) 'Transform X lines of C array'
! Transform X lines of C array
      !!IW = 2*L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk)) + 3
      IW= L_len-1
!!      write(*,'(*(a,i0))') 'IW ',IW,' L ',L,' M ',M,' LDIM ',LDIM
 !!     CALL CFFTMF(L,    1, M, LDIM, C, (L-1) + LDIM*(M-1) +1,WSAVE(IW), M_len ,WORK, 2*L*M, IER1)
 !!   CALL CFFTMB(L,    1, M, LDIM, C, (L-1) + LDIM*(M-1) +1,WSAVE(IW), M_len ,WORK, 2*L*M, IER1) ! sieht genauso aus in CFFT2B
      CALL CFFTMF(L,    1, M, LDIM, C, (L-1) + LDIM*(M-1) +1, WORK, 2*L*M,fft_sign(2), IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('CFFT2F',-5)
        GO TO 100
      ENDIF
!
!!     write(*,*) 'Transform Y lines of C array'
! Transform Y lines of C array
      IW = 1
 !!     CALL CFFTMF(M, LDIM, L,    1, C, (M-1)*LDIM + L       ,WSAVE(IW), L_len ,WORK, 2*M*L, IER1)
 !!   CALL CFFTMB(M, LDIM, L,    1, C, (M-1)*LDIM + L       ,WSAVE(IW), L_len ,WORK, 2*M*L, IER1) ! sieht genauso aus in CFFT2B
      CALL CFFTMF(M, LDIM, L,    1, C, (M-1)*LDIM + L       , WORK, 2*M*L,fft_sign(1), IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('CFFT2F',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine CFFT2F
!************************************************************************************************************************
      SUBROUTINE CFFT3F ( L1, L2, L3, C, LENC, fft_sign, IER)
      integer                     :: L1, L2, L3
      integer                     :: LDIM
      !Kuester COMPLEX(kind=rk)   :: C(LDIM,M)
      integer                     :: IER,IER1
      !!real(kind=rk)     :: C(2*LDIM,M)
      real(kind=rk)               :: C(*)
      integer                     :: LOT,JUMP,NN,INC,LENC
      integer                     :: ii
      integer                     :: offset
      integer                     :: total_dim
      integer                     :: LENWRK
      real(kind=rk),allocatable,dimension(:) :: WORK
      type(fft_real_sign_type) :: fft_sign(:)
!
! Initialize error return
!
      IER = 0
!
!
      total_dim=L1*L2*L3
      !!LENC=size(C)
      LENWRK=2*total_dim
      allocate(WORK(LENWRK))

         write(*,*) 'Transform along dimension 1 of C array'
         NN=L1;JUMP=L1;INC=1;LOT=L2*L3 
         offset=1
         call   CFFTMF(LOT, JUMP, NN, INC, C(2*offset), LENC, WORK, LENWRK, fft_sign(1), IER)
         IF (IER1 .NE. 0) THEN
           IER = 20
           CALL XERFFT ('CFFT3F',-5)
           GO TO 100
         ENDIF

         NN=L2;JUMP=1;INC=1;LOT=L1  
      do ii=1,L3
         offset=1+L1*L2*(ii-1)
         write(*,*) 'Transform along dimension 2 of C array ',ii
         call   CFFTMF(LOT, JUMP, NN, INC, C(2*offset), LENC, WORK, LENWRK, fft_sign(2), IER)
         IF (IER1 .NE. 0) THEN
           IER = 20
           CALL XERFFT ('CFFT3F',-5)
           GO TO 100
         ENDIF
      enddo

         NN=L3;JUMP=1;INC=L1*L2;LOT=L1*L2
         offset=1
         call   CFFTMF(LOT, JUMP, NN, INC, C(2*offset), LENC, WORK, LENWRK, fft_sign(3), IER)
         write(*,*) 'Transform along dimension 3 of C array '
         IF (IER1 .NE. 0) THEN
           IER = 20
           CALL XERFFT ('CFFT3F',-5)
           GO TO 100
         ENDIF

  100 CONTINUE
      END subroutine CFFT3F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CFFT2I (L, M, fft_sign, IER)
      integer           :: LENSAV
      INTEGER           :: L, M, IER
      integer           :: IER1
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      integer           :: L_len,M_len
      type(fft_real_sign_type) :: fft_sign(:)
!
! Initialize error return
!
      IER = 0

      !!L_len=2*L + INT(LOG(REAL(L,kind=rk))/LOG(2._rk)) + 4
      !!M_len=2*M + INT(LOG(REAL(M,kind=rk))/LOG(2._rk)) + 4
!
      !!IF (LENSAV .LT. 2*L + INT(LOG(REAL(L,kind=rk))/LOG(2._rk)) +2*M + INT(LOG(REAL(M,kind=rk))/LOG(2._rk)) +8) THEN
      !!IF (LENSAV .LT. L_len + M_len ) THEN
      !!  IER = 2
      !!  CALL XERFFT ('CFFT2I', 4)
      !!  GO TO 100
      !!ENDIF
!
      !!CALL CFFTMI_orig (L, WSAVE(1), 2*L + INT(LOG(REAL(L,kind=rk))/LOG(2._rk)) + 4,IER1)
      !!CALL CFFTMI_orig (L, WSAVE(1), L_len ,IER1)
      !!IF (IER1 .NE. 0) THEN
      !!  IER = 20
      !!  CALL XERFFT ('CFFT2I',-5)
      !!  GO TO 100
      !!ENDIF

      !!CALL CFFTMI_orig (M, WSAVE(2*L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk)) + 3),2*M + INT(LOG(REAL(M,kind=rk))/LOG(2._rk)) + 4, IER1)
      !!CALL CFFTMI_orig (M, WSAVE(L_len-1),M_len, IER1)
      !!IF (IER1 .NE. 0) THEN
      !!  IER = 20
      !!  CALL XERFFT ('CFFT2I',-5)
      !!ENDIF

      call CFFTMI (L, fft_sign(1))
      call CFFTMI (M, fft_sign(2))
!
  100 CONTINUE
      END subroutine CFFT2I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFTMB (LOT, JUMP, N, INC, C, LENC, WSAVE, LENSAV,WORK, LENWRK, IER)
      SUBROUTINE CFFTMB (LOT, JUMP, N, INC, C, LENC, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENC
      integer     ::  LENWRK, IER
      integer     :: LENSAV
      integer           :: IW1
      !Kuester COMPLEX(kind=rk)   :: C(LENC)
      real(kind=rk)     :: C(2*LENC)
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      integer           :: N_len
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL XERCON
!
      IER = 0
      !!N_len=2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
!
      IF (LENC .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('CFFTMB ', 6)
     !! ELSEIF (LENSAV .LT. N_len) THEN
     !!   IER = 2
     !!   CALL XERFFT ('CFFTMB ', 8)
      ELSEIF (LENWRK .LT. 2*LOT*N) THEN
        IER = 3
        CALL XERFFT ('CFFTMB ', 10)
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('CFFTMB ', -1)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!IW1 = N+N+1
      !!CALL CMFM1B (LOT,JUMP,N,INC,C,WORK,WSAVE,WSAVE(IW1),WSAVE(IW1+1))
      CALL CMFM1B (LOT,JUMP,N,INC,C,WORK,fft_sign)
      END subroutine CFFTMB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CFFTMF (LOT, JUMP, N, INC, C, LENC, WSAVE, LENSAV,WORK, LENWRK, IER)
      SUBROUTINE CFFTMF (LOT, JUMP, N, INC, C, LENC, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENC
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IW1
      !Kuester COMPLEX(kind=rk)   :: C(LENC)
      real(kind=rk)     :: C(2*LENC)
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      integer           :: N_len
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL  XERCON
!
      IER = 0
      !!N_len=2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
!
      IF (LENC .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('CFFTMF ', 6)
      !!ELSEIF (LENSAV .LT. N_len ) THEN
      !!  IER = 2
      !!  CALL XERFFT ('CFFTMF ', 8)
      ELSEIF (LENWRK .LT. 2*LOT*N) THEN
        IER = 3
        CALL XERFFT ('CFFTMF ', 10)
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('CFFTMF ', -1)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!IW1 = N+N+1
      !!CALL CMFM1F (LOT,JUMP,N,INC,C,WORK,WSAVE,WSAVE(IW1),WSAVE(IW1+1))
      CALL CMFM1F (LOT,JUMP,N,INC,C,WORK,fft_sign)
      END subroutine CFFTMF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CFFTMI_orig (N, WSAVE, LENSAV,IER)
      integer     :: N
      integer    ::  IER
      integer   :: LENSAV
      integer           :: IW1
      REAL(kind=rk)     :: WSAVE(LENSAV)
      integer           :: N_len
!
      IER = 0
      N_len=2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
!
      IF (LENSAV .LT. N_len ) THEN
        IER = 2
        CALL XERFFT ('CFFTMI ', 3)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      IW1 = N+N+1
      CALL MCFTI1_orig (N,WSAVE,WSAVE(IW1),WSAVE(IW1+1))
      END subroutine CFFTMI_orig
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CFFTMI (N, fft_sign)
      integer     :: N
      integer     ::  IER
      integer     :: LENSAV
      integer           :: IW1
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
      !!stop' mit WSAVE, LENSAV,fft_sign falsch'
!
     !! IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4) THEN
     !!   IER = 2
     !!   CALL XERFFT ('CFFTMI ', 3)
     !! ENDIF
!
      !!IF (N .EQ. 1) RETURN
!
      !!IW1 = N+N+1
      !!CALL MCFTI1 (N,WSAVE,WSAVE(IW1),WSAVE(IW1+1),fft_sign)
      CALL MCFTI1 (N,fft_sign)
      END subroutine CFFTMI
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF2KB (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,2),CH(2,IN2,L1,2,IDO),WA(IDO,1,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
!
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         CHOLD1 = CC(1,M1,K,1,1)+CC(1,M1,K,1,2)
         CC(1,M1,K,1,2) = CC(1,M1,K,1,1)-CC(1,M1,K,1,2)
         CC(1,M1,K,1,1) = CHOLD1
         CHOLD2 = CC(2,M1,K,1,1)+CC(2,M1,K,1,2)
         CC(2,M1,K,1,2) = CC(2,M1,K,1,1)-CC(2,M1,K,1,2)
         CC(2,M1,K,1,1) = CHOLD2
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(1,M2,K,1,1) = CC(1,M1,K,1,1)+CC(1,M1,K,1,2)
         CH(1,M2,K,2,1) = CC(1,M1,K,1,1)-CC(1,M1,K,1,2)
         CH(2,M2,K,1,1) = CC(2,M1,K,1,1)+CC(2,M1,K,1,2)
         CH(2,M2,K,2,1) = CC(2,M1,K,1,1)-CC(2,M1,K,1,2)
  103 CONTINUE
      IF(IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
         DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            CH(1,M2,K,1,I) = CC(1,M1,K,I,1)+CC(1,M1,K,I,2)
            TR2 = CC(1,M1,K,I,1)-CC(1,M1,K,I,2)
            CH(2,M2,K,1,I) = CC(2,M1,K,I,1)+CC(2,M1,K,I,2)
            TI2 = CC(2,M1,K,I,1)-CC(2,M1,K,I,2)
            CH(2,M2,K,2,I) = WA(I,1,1)*TI2+WA(I,1,2)*TR2
            CH(1,M2,K,2,I) = WA(I,1,1)*TR2-WA(I,1,2)*TI2
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF2KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF2KF (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: SN
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,2),CH(2,IN2,L1,2,IDO),WA(IDO,1,2)
!
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(2*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         CHOLD1 = SN*(CC(1,M1,K,1,1)+CC(1,M1,K,1,2))
         CC(1,M1,K,1,2) = SN*(CC(1,M1,K,1,1)-CC(1,M1,K,1,2))
         CC(1,M1,K,1,1) = CHOLD1
         CHOLD2 = SN*(CC(2,M1,K,1,1)+CC(2,M1,K,1,2))
         CC(2,M1,K,1,2) = SN*(CC(2,M1,K,1,1)-CC(2,M1,K,1,2))
         CC(2,M1,K,1,1) = CHOLD2
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         M2 = M2S
         DO 107 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(1,M2,K,1,1) = SN*(CC(1,M1,K,1,1)+CC(1,M1,K,1,2))
         CH(1,M2,K,2,1) = SN*(CC(1,M1,K,1,1)-CC(1,M1,K,1,2))
         CH(2,M2,K,1,1) = SN*(CC(2,M1,K,1,1)+CC(2,M1,K,1,2))
         CH(2,M2,K,2,1) = SN*(CC(2,M1,K,1,1)-CC(2,M1,K,1,2))
  107 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(1,M2,K,1,1) = CC(1,M1,K,1,1)+CC(1,M1,K,1,2)
         CH(1,M2,K,2,1) = CC(1,M1,K,1,1)-CC(1,M1,K,1,2)
         CH(2,M2,K,1,1) = CC(2,M1,K,1,1)+CC(2,M1,K,1,2)
         CH(2,M2,K,2,1) = CC(2,M1,K,1,1)-CC(2,M1,K,1,2)
  103 CONTINUE
      DO 105 I=2,IDO
         DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            CH(1,M2,K,1,I) = CC(1,M1,K,I,1)+CC(1,M1,K,I,2)
            TR2 = CC(1,M1,K,I,1)-CC(1,M1,K,I,2)
            CH(2,M2,K,1,I) = CC(2,M1,K,I,1)+CC(2,M1,K,I,2)
            TI2 = CC(2,M1,K,I,1)-CC(2,M1,K,I,2)
            CH(2,M2,K,2,I) = WA(I,1,1)*TI2-WA(I,1,2)*TR2
            CH(1,M2,K,2,I) = WA(I,1,1)*TR2+WA(I,1,2)*TI2
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF2KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF3KB (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,3),CH(2,IN2,L1,3,IDO),WA(IDO,2,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TAUR,TAUI
      !  (TAUR,TAUI)=exp(+2*PI/3*i)
      DATA TAUR,TAUI /-.5_rk,.866025403784439_rk/
!
!!      write(*,*) 'CMF3KB'
      TAUR=real(root3,kind=rk);TAUI=-aimag(root3)

      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,3)
         CR2 = CC(1,M1,K,1,1)+TAUR*TR2
         CC(1,M1,K,1,1) = CC(1,M1,K,1,1)+TR2
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,3)
         CI2 = CC(2,M1,K,1,1)+TAUR*TI2
         CC(2,M1,K,1,1) = CC(2,M1,K,1,1)+TI2
         CR3 = TAUI*(CC(1,M1,K,1,2)-CC(1,M1,K,1,3))
         CI3 = TAUI*(CC(2,M1,K,1,2)-CC(2,M1,K,1,3))
         CC(1,M1,K,1,2) = CR2-CI3
         CC(1,M1,K,1,3) = CR2+CI3
         CC(2,M1,K,1,2) = CI2+CR3
         CC(2,M1,K,1,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,3)
         CR2 = CC(1,M1,K,1,1)+TAUR*TR2
         CH(1,M2,K,1,1) = CC(1,M1,K,1,1)+TR2
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,3)
         CI2 = CC(2,M1,K,1,1)+TAUR*TI2
         CH(2,M2,K,1,1) = CC(2,M1,K,1,1)+TI2
         CR3 = TAUI*(CC(1,M1,K,1,2)-CC(1,M1,K,1,3))
         CI3 = TAUI*(CC(2,M1,K,1,2)-CC(2,M1,K,1,3))
         CH(1,M2,K,2,1) = CR2-CI3
         CH(1,M2,K,3,1) = CR2+CI3
         CH(2,M2,K,2,1) = CI2+CR3
         CH(2,M2,K,3,1) = CI2-CR3
  103 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
        DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            TR2 = CC(1,M1,K,I,2)+CC(1,M1,K,I,3)
            CR2 = CC(1,M1,K,I,1)+TAUR*TR2
            CH(1,M2,K,1,I) = CC(1,M1,K,I,1)+TR2
            TI2 = CC(2,M1,K,I,2)+CC(2,M1,K,I,3)
            CI2 = CC(2,M1,K,I,1)+TAUR*TI2
            CH(2,M2,K,1,I) = CC(2,M1,K,I,1)+TI2
            CR3 = TAUI*(CC(1,M1,K,I,2)-CC(1,M1,K,I,3))
            CI3 = TAUI*(CC(2,M1,K,I,2)-CC(2,M1,K,I,3))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(2,M2,K,2,I) = WA(I,1,1)*DI2+WA(I,1,2)*DR2
            CH(1,M2,K,2,I) = WA(I,1,1)*DR2-WA(I,1,2)*DI2
            CH(2,M2,K,3,I) = WA(I,2,1)*DI3+WA(I,2,2)*DR3
            CH(1,M2,K,3,I) = WA(I,2,1)*DR3-WA(I,2,2)*DI3
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF3KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF3KF (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,3),CH(2,IN2,L1,3,IDO),WA(IDO,2,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TAUR,TAUI
      real(kind=rk)     :: SN
      !  (TAUR,TAUI)=exp(-2*PI/3*i)
      DATA TAUR,TAUI /-.5_rk,-.866025403784439_rk/
!
!!      write(*,*) 'CMF3KF'
      TAUR=real(root3,kind=rk);TAUI=aimag(root3)

      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(3*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,3)
         CR2 = CC(1,M1,K,1,1)+TAUR*TR2
         CC(1,M1,K,1,1) = SN*(CC(1,M1,K,1,1)+TR2)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,3)
         CI2 = CC(2,M1,K,1,1)+TAUR*TI2
         CC(2,M1,K,1,1) = SN*(CC(2,M1,K,1,1)+TI2)
         CR3 = TAUI*(CC(1,M1,K,1,2)-CC(1,M1,K,1,3))
         CI3 = TAUI*(CC(2,M1,K,1,2)-CC(2,M1,K,1,3))
         CC(1,M1,K,1,2) = SN*(CR2-CI3)
         CC(1,M1,K,1,3) = SN*(CR2+CI3)
         CC(2,M1,K,1,2) = SN*(CI2+CR3)
         CC(2,M1,K,1,3) = SN*(CI2-CR3)
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         M2 = M2S
         DO 107 M1=1,M1D,IM1
         M2 = M2+IM2
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,3)
         CR2 = CC(1,M1,K,1,1)+TAUR*TR2
         CH(1,M2,K,1,1) = SN*(CC(1,M1,K,1,1)+TR2)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,3)
         CI2 = CC(2,M1,K,1,1)+TAUR*TI2
         CH(2,M2,K,1,1) = SN*(CC(2,M1,K,1,1)+TI2)
         CR3 = TAUI*(CC(1,M1,K,1,2)-CC(1,M1,K,1,3))
         CI3 = TAUI*(CC(2,M1,K,1,2)-CC(2,M1,K,1,3))
         CH(1,M2,K,2,1) = SN*(CR2-CI3)
         CH(1,M2,K,3,1) = SN*(CR2+CI3)
         CH(2,M2,K,2,1) = SN*(CI2+CR3)
         CH(2,M2,K,3,1) = SN*(CI2-CR3)
  107 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,3)
         CR2 = CC(1,M1,K,1,1)+TAUR*TR2
         CH(1,M2,K,1,1) = CC(1,M1,K,1,1)+TR2
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,3)
         CI2 = CC(2,M1,K,1,1)+TAUR*TI2
         CH(2,M2,K,1,1) = CC(2,M1,K,1,1)+TI2
         CR3 = TAUI*(CC(1,M1,K,1,2)-CC(1,M1,K,1,3))
         CI3 = TAUI*(CC(2,M1,K,1,2)-CC(2,M1,K,1,3))
         CH(1,M2,K,2,1) = CR2-CI3
         CH(1,M2,K,3,1) = CR2+CI3
         CH(2,M2,K,2,1) = CI2+CR3
         CH(2,M2,K,3,1) = CI2-CR3
  103 CONTINUE
      DO 105 I=2,IDO
        DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            TR2 = CC(1,M1,K,I,2)+CC(1,M1,K,I,3)
            CR2 = CC(1,M1,K,I,1)+TAUR*TR2
            CH(1,M2,K,1,I) = CC(1,M1,K,I,1)+TR2
            TI2 = CC(2,M1,K,I,2)+CC(2,M1,K,I,3)
            CI2 = CC(2,M1,K,I,1)+TAUR*TI2
            CH(2,M2,K,1,I) = CC(2,M1,K,I,1)+TI2
            CR3 = TAUI*(CC(1,M1,K,I,2)-CC(1,M1,K,I,3))
            CI3 = TAUI*(CC(2,M1,K,I,2)-CC(2,M1,K,I,3))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(2,M2,K,2,I) = WA(I,1,1)*DI2-WA(I,1,2)*DR2
            CH(1,M2,K,2,I) = WA(I,1,1)*DR2+WA(I,1,2)*DI2
            CH(2,M2,K,3,I) = WA(I,2,1)*DI3-WA(I,2,2)*DR3
            CH(1,M2,K,3,I) = WA(I,2,1)*DR3+WA(I,2,2)*DI3
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF3KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF4KB (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,4),CH(2,IN2,L1,4,IDO),WA(IDO,3,2)
!
! FFTPACK 5.0 auxiliary routine
!
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         TI1 = CC(2,M1,K,1,1)-CC(2,M1,K,1,3)
         TI2 = CC(2,M1,K,1,1)+CC(2,M1,K,1,3)
         TR4 = CC(2,M1,K,1,4)-CC(2,M1,K,1,2)
         TI3 = CC(2,M1,K,1,2)+CC(2,M1,K,1,4)
         TR1 = CC(1,M1,K,1,1)-CC(1,M1,K,1,3)
         TR2 = CC(1,M1,K,1,1)+CC(1,M1,K,1,3)
         TI4 = CC(1,M1,K,1,2)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,2)+CC(1,M1,K,1,4)
         CC(1,M1,K,1,1) = TR2+TR3
         CC(1,M1,K,1,3) = TR2-TR3
         CC(2,M1,K,1,1) = TI2+TI3
         CC(2,M1,K,1,3) = TI2-TI3
         CC(1,M1,K,1,2) = TR1+TR4
         CC(1,M1,K,1,4) = TR1-TR4
         CC(2,M1,K,1,2) = TI1+TI4
         CC(2,M1,K,1,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         TI1 = CC(2,M1,K,1,1)-CC(2,M1,K,1,3)
         TI2 = CC(2,M1,K,1,1)+CC(2,M1,K,1,3)
         TR4 = CC(2,M1,K,1,4)-CC(2,M1,K,1,2)
         TI3 = CC(2,M1,K,1,2)+CC(2,M1,K,1,4)
         TR1 = CC(1,M1,K,1,1)-CC(1,M1,K,1,3)
         TR2 = CC(1,M1,K,1,1)+CC(1,M1,K,1,3)
         TI4 = CC(1,M1,K,1,2)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,2)+CC(1,M1,K,1,4)
         CH(1,M2,K,1,1) = TR2+TR3
         CH(1,M2,K,3,1) = TR2-TR3
         CH(2,M2,K,1,1) = TI2+TI3
         CH(2,M2,K,3,1) = TI2-TI3
         CH(1,M2,K,2,1) = TR1+TR4
         CH(1,M2,K,4,1) = TR1-TR4
         CH(2,M2,K,2,1) = TI1+TI4
         CH(2,M2,K,4,1) = TI1-TI4
  103 CONTINUE
      IF(IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
         DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            TI1 = CC(2,M1,K,I,1)-CC(2,M1,K,I,3)
            TI2 = CC(2,M1,K,I,1)+CC(2,M1,K,I,3)
            TI3 = CC(2,M1,K,I,2)+CC(2,M1,K,I,4)
            TR4 = CC(2,M1,K,I,4)-CC(2,M1,K,I,2)
            TR1 = CC(1,M1,K,I,1)-CC(1,M1,K,I,3)
            TR2 = CC(1,M1,K,I,1)+CC(1,M1,K,I,3)
            TI4 = CC(1,M1,K,I,2)-CC(1,M1,K,I,4)
            TR3 = CC(1,M1,K,I,2)+CC(1,M1,K,I,4)
            CH(1,M2,K,1,I) = TR2+TR3
            CR3 = TR2-TR3
            CH(2,M2,K,1,I) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(1,M2,K,2,I) = WA(I,1,1)*CR2-WA(I,1,2)*CI2
            CH(2,M2,K,2,I) = WA(I,1,1)*CI2+WA(I,1,2)*CR2
            CH(1,M2,K,3,I) = WA(I,2,1)*CR3-WA(I,2,2)*CI3
            CH(2,M2,K,3,I) = WA(I,2,1)*CI3+WA(I,2,2)*CR3
            CH(1,M2,K,4,I) = WA(I,3,1)*CR4-WA(I,3,2)*CI4
            CH(2,M2,K,4,I) = WA(I,3,1)*CI4+WA(I,3,2)*CR4
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF4KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
   !  54    CALL CMF4KF (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
      SUBROUTINE CMF4KF (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: SN
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,4),CH(2,IN2,L1,4,IDO),WA(IDO,3,2)
!
! FFTPACK 5.0 auxiliary routine
!
     write(*,*) 'CMF4KF: size(CC)=',size(CC)
!!      write(*,'(a,*(1x,i0))') 'CMF4KF(LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)=',LOT,IDO,L1,NA,IM1,IN1,IM2,IN2
      write(*,*) 'IN1 ',IN1,' L1 ',L1,' IDO ',IDO,' IN1*L1 ',IN1*L1
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      write(*,'(*(a,i0))') 'M1D ',M1D,' IM1 ',IM1,' M2S ',M2S,' NA ',NA
      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(4*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         TI1 = CC(2,M1,K,1,1)-CC(2,M1,K,1,3)
         TI2 = CC(2,M1,K,1,1)+CC(2,M1,K,1,3)
         TR4 = CC(2,M1,K,1,2)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,2)+CC(2,M1,K,1,4)
         TR1 = CC(1,M1,K,1,1)-CC(1,M1,K,1,3)
         TR2 = CC(1,M1,K,1,1)+CC(1,M1,K,1,3)
         TI4 = CC(1,M1,K,1,4)-CC(1,M1,K,1,2)
         TR3 = CC(1,M1,K,1,2)+CC(1,M1,K,1,4)
         CC(1,M1,K,1,1) = SN*(TR2+TR3)
         CC(1,M1,K,1,3) = SN*(TR2-TR3)
         CC(2,M1,K,1,1) = SN*(TI2+TI3)
         CC(2,M1,K,1,3) = SN*(TI2-TI3)
         CC(1,M1,K,1,2) = SN*(TR1+TR4)
         CC(1,M1,K,1,4) = SN*(TR1-TR4)
         CC(2,M1,K,1,2) = SN*(TI1+TI4)
         CC(2,M1,K,1,4) = SN*(TI1-TI4)
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         M2 = M2S
         DO 107 M1=1,M1D,IM1
         M2 = M2+IM2
         TI1 = CC(2,M1,K,1,1)-CC(2,M1,K,1,3)
         TI2 = CC(2,M1,K,1,1)+CC(2,M1,K,1,3)
         TR4 = CC(2,M1,K,1,2)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,2)+CC(2,M1,K,1,4)
         TR1 = CC(1,M1,K,1,1)-CC(1,M1,K,1,3)
         TR2 = CC(1,M1,K,1,1)+CC(1,M1,K,1,3)
         TI4 = CC(1,M1,K,1,4)-CC(1,M1,K,1,2)
         TR3 = CC(1,M1,K,1,2)+CC(1,M1,K,1,4)
         CH(1,M2,K,1,1) = SN*(TR2+TR3)
         CH(1,M2,K,3,1) = SN*(TR2-TR3)
         CH(2,M2,K,1,1) = SN*(TI2+TI3)
         CH(2,M2,K,3,1) = SN*(TI2-TI3)
         CH(1,M2,K,2,1) = SN*(TR1+TR4)
         CH(1,M2,K,4,1) = SN*(TR1-TR4)
         CH(2,M2,K,2,1) = SN*(TI1+TI4)
         CH(2,M2,K,4,1) = SN*(TI1-TI4)
  107 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         TI1 = CC(2,M1,K,1,1)-CC(2,M1,K,1,3)
         TI2 = CC(2,M1,K,1,1)+CC(2,M1,K,1,3)
         TR4 = CC(2,M1,K,1,2)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,2)+CC(2,M1,K,1,4)
         TR1 = CC(1,M1,K,1,1)-CC(1,M1,K,1,3)
         TR2 = CC(1,M1,K,1,1)+CC(1,M1,K,1,3)
         TI4 = CC(1,M1,K,1,4)-CC(1,M1,K,1,2)
         TR3 = CC(1,M1,K,1,2)+CC(1,M1,K,1,4)
         CH(1,M2,K,1,1) = TR2+TR3
         CH(1,M2,K,3,1) = TR2-TR3
         CH(2,M2,K,1,1) = TI2+TI3
         CH(2,M2,K,3,1) = TI2-TI3
         CH(1,M2,K,2,1) = TR1+TR4
         CH(1,M2,K,4,1) = TR1-TR4
         CH(2,M2,K,2,1) = TI1+TI4
         CH(2,M2,K,4,1) = TI1-TI4
  103 CONTINUE
      DO 105 I=2,IDO
         DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            TI1 = CC(2,M1,K,I,1)-CC(2,M1,K,I,3)
            TI2 = CC(2,M1,K,I,1)+CC(2,M1,K,I,3)
            TI3 = CC(2,M1,K,I,2)+CC(2,M1,K,I,4)
            TR4 = CC(2,M1,K,I,2)-CC(2,M1,K,I,4)
            TR1 = CC(1,M1,K,I,1)-CC(1,M1,K,I,3)
            TR2 = CC(1,M1,K,I,1)+CC(1,M1,K,I,3)
            TI4 = CC(1,M1,K,I,4)-CC(1,M1,K,I,2)
            TR3 = CC(1,M1,K,I,2)+CC(1,M1,K,I,4)
            CH(1,M2,K,1,I) = TR2+TR3
            CR3 = TR2-TR3
            CH(2,M2,K,1,I) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(1,M2,K,2,I) = WA(I,1,1)*CR2+WA(I,1,2)*CI2
            CH(2,M2,K,2,I) = WA(I,1,1)*CI2-WA(I,1,2)*CR2
            CH(1,M2,K,3,I) = WA(I,2,1)*CR3+WA(I,2,2)*CI3
            CH(2,M2,K,3,I) = WA(I,2,1)*CI3-WA(I,2,2)*CR3
            CH(1,M2,K,4,I) = WA(I,3,1)*CR4+WA(I,3,2)*CI4
            CH(2,M2,K,4,I) = WA(I,3,1)*CI4-WA(I,3,2)*CR4
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF4KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF5KB (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,5),CH(2,IN2,L1,5,IDO),WA(IDO,4,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      real(kind=rk)     :: CHOLD1,CHOLD2
      !   (TR11,TI11)=exp(+2*PI/5*i), (TR12,TI12)=exp(+4*PI/5*i)
      DATA TR11,TI11,TR12,TI12 /.3090169943749474_rk,.9510565162951536_rk,-.8090169943749474_rk,.5877852522924731_rk/
!
! FFTPACK 5.0 auxiliary routine
!
!!      write(*,*) 'CMF5KB'
      TR11=real(root51,kind=rk);TI11=-aimag(root51)
      TR12=real(root52,kind=rk);TI12=-aimag(root52)


      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1 .OR. NA .EQ. 1) GO TO 102
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         TI5 = CC(2,M1,K,1,2)-CC(2,M1,K,1,5)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,5)
         TI4 = CC(2,M1,K,1,3)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,3)+CC(2,M1,K,1,4)
         TR5 = CC(1,M1,K,1,2)-CC(1,M1,K,1,5)
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,5)
         TR4 = CC(1,M1,K,1,3)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,3)+CC(1,M1,K,1,4)
         CHOLD1 = CC(1,M1,K,1,1)+TR2+TR3
         CHOLD2 = CC(2,M1,K,1,1)+TI2+TI3
         CR2 = CC(1,M1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,M1,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,M1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,M1,K,1,1)+TR12*TI2+TR11*TI3
         CC(1,M1,K,1,1) = CHOLD1
         CC(2,M1,K,1,1) = CHOLD2
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CC(1,M1,K,1,2) = CR2-CI5
         CC(1,M1,K,1,5) = CR2+CI5
         CC(2,M1,K,1,2) = CI2+CR5
         CC(2,M1,K,1,3) = CI3+CR4
         CC(1,M1,K,1,3) = CR3-CI4
         CC(1,M1,K,1,4) = CR3+CI4
         CC(2,M1,K,1,4) = CI3-CR4
         CC(2,M1,K,1,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         TI5 = CC(2,M1,K,1,2)-CC(2,M1,K,1,5)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,5)
         TI4 = CC(2,M1,K,1,3)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,3)+CC(2,M1,K,1,4)
         TR5 = CC(1,M1,K,1,2)-CC(1,M1,K,1,5)
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,5)
         TR4 = CC(1,M1,K,1,3)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,3)+CC(1,M1,K,1,4)
         CH(1,M2,K,1,1) = CC(1,M1,K,1,1)+TR2+TR3
         CH(2,M2,K,1,1) = CC(2,M1,K,1,1)+TI2+TI3
         CR2 = CC(1,M1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,M1,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,M1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,M1,K,1,1)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,M2,K,2,1) = CR2-CI5
         CH(1,M2,K,5,1) = CR2+CI5
         CH(2,M2,K,2,1) = CI2+CR5
         CH(2,M2,K,3,1) = CI3+CR4
         CH(1,M2,K,3,1) = CR3-CI4
         CH(1,M2,K,4,1) = CR3+CI4
         CH(2,M2,K,4,1) = CI3-CR4
         CH(2,M2,K,5,1) = CI2-CR5
  103 CONTINUE
      IF(IDO .EQ. 1) RETURN
      DO 105 I=2,IDO
         DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            TI5 = CC(2,M1,K,I,2)-CC(2,M1,K,I,5)
            TI2 = CC(2,M1,K,I,2)+CC(2,M1,K,I,5)
            TI4 = CC(2,M1,K,I,3)-CC(2,M1,K,I,4)
            TI3 = CC(2,M1,K,I,3)+CC(2,M1,K,I,4)
            TR5 = CC(1,M1,K,I,2)-CC(1,M1,K,I,5)
            TR2 = CC(1,M1,K,I,2)+CC(1,M1,K,I,5)
            TR4 = CC(1,M1,K,I,3)-CC(1,M1,K,I,4)
            TR3 = CC(1,M1,K,I,3)+CC(1,M1,K,I,4)
            CH(1,M2,K,1,I) = CC(1,M1,K,I,1)+TR2+TR3
            CH(2,M2,K,1,I) = CC(2,M1,K,I,1)+TI2+TI3
            CR2 = CC(1,M1,K,I,1)+TR11*TR2+TR12*TR3
            CI2 = CC(2,M1,K,I,1)+TR11*TI2+TR12*TI3
            CR3 = CC(1,M1,K,I,1)+TR12*TR2+TR11*TR3
            CI3 = CC(2,M1,K,I,1)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(1,M2,K,2,I) = WA(I,1,1)*DR2-WA(I,1,2)*DI2
            CH(2,M2,K,2,I) = WA(I,1,1)*DI2+WA(I,1,2)*DR2
            CH(1,M2,K,3,I) = WA(I,2,1)*DR3-WA(I,2,2)*DI3
            CH(2,M2,K,3,I) = WA(I,2,1)*DI3+WA(I,2,2)*DR3
            CH(1,M2,K,4,I) = WA(I,3,1)*DR4-WA(I,3,2)*DI4
            CH(2,M2,K,4,I) = WA(I,3,1)*DI4+WA(I,3,2)*DR4
            CH(1,M2,K,5,I) = WA(I,4,1)*DR5-WA(I,4,2)*DI5
            CH(2,M2,K,5,I) = WA(I,4,1)*DI5+WA(I,4,2)*DR5
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF5KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMF5KF (LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
      integer           :: LOT,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=rk)     :: CC(2,IN1,L1,IDO,5),CH(2,IN2,L1,5,IDO),WA(IDO,4,2)
      real(kind=rk)     :: TR1,TI1,TR2,TI2,TR3,TI3,TR4,TI4,TR5,TI5
      real(kind=rk)     :: CR1,CI1,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5
      real(kind=rk)     :: DR1,DI1,DR2,DI2,DR3,DI3,DR4,DI4,DR5,DI5
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      real(kind=rk)     :: CHOLD1,CHOLD2
      real(kind=rk)     :: SN
      !   (TR11,TI11)=exp(-2*PI/5*i), (TR12,TI12)=exp(-4*PI/5*i)
      DATA TR11,TI11,TR12,TI12 /.3090169943749474_rk,-.9510565162951536_rk,-.8090169943749474_rk,-.5877852522924731_rk/
!
! FFTPACK 5.0 auxiliary routine
!
!!      write(*,*) 'CMF5KF'
      TR11=real(root51,kind=rk);TI11=aimag(root51)
      TR12=real(root52,kind=rk);TI12=aimag(root52)

      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IF (IDO .GT. 1) GO TO 102
      SN = 1._rk/REAL(5*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 106
      DO 101 K=1,L1
         DO 101 M1=1,M1D,IM1
         TI5 = CC(2,M1,K,1,2)-CC(2,M1,K,1,5)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,5)
         TI4 = CC(2,M1,K,1,3)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,3)+CC(2,M1,K,1,4)
         TR5 = CC(1,M1,K,1,2)-CC(1,M1,K,1,5)
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,5)
         TR4 = CC(1,M1,K,1,3)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,3)+CC(1,M1,K,1,4)
         CHOLD1 = SN*(CC(1,M1,K,1,1)+TR2+TR3)
         CHOLD2 = SN*(CC(2,M1,K,1,1)+TI2+TI3)
         CR2 = CC(1,M1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,M1,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,M1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,M1,K,1,1)+TR12*TI2+TR11*TI3
         CC(1,M1,K,1,1) = CHOLD1
         CC(2,M1,K,1,1) = CHOLD2
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CC(1,M1,K,1,2) = SN*(CR2-CI5)
         CC(1,M1,K,1,5) = SN*(CR2+CI5)
         CC(2,M1,K,1,2) = SN*(CI2+CR5)
         CC(2,M1,K,1,3) = SN*(CI3+CR4)
         CC(1,M1,K,1,3) = SN*(CR3-CI4)
         CC(1,M1,K,1,4) = SN*(CR3+CI4)
         CC(2,M1,K,1,4) = SN*(CI3-CR4)
         CC(2,M1,K,1,5) = SN*(CI2-CR5)
  101 CONTINUE
      RETURN
  106 DO 107 K=1,L1
         M2 = M2S
         DO 107 M1=1,M1D,IM1
         M2 = M2+IM2
         TI5 = CC(2,M1,K,1,2)-CC(2,M1,K,1,5)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,5)
         TI4 = CC(2,M1,K,1,3)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,3)+CC(2,M1,K,1,4)
         TR5 = CC(1,M1,K,1,2)-CC(1,M1,K,1,5)
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,5)
         TR4 = CC(1,M1,K,1,3)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,3)+CC(1,M1,K,1,4)
         CH(1,M2,K,1,1) = SN*(CC(1,M1,K,1,1)+TR2+TR3)
         CH(2,M2,K,1,1) = SN*(CC(2,M1,K,1,1)+TI2+TI3)
         CR2 = CC(1,M1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,M1,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,M1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,M1,K,1,1)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,M2,K,2,1) = SN*(CR2-CI5)
         CH(1,M2,K,5,1) = SN*(CR2+CI5)
         CH(2,M2,K,2,1) = SN*(CI2+CR5)
         CH(2,M2,K,3,1) = SN*(CI3+CR4)
         CH(1,M2,K,3,1) = SN*(CR3-CI4)
         CH(1,M2,K,4,1) = SN*(CR3+CI4)
         CH(2,M2,K,4,1) = SN*(CI3-CR4)
         CH(2,M2,K,5,1) = SN*(CI2-CR5)
  107 CONTINUE
      RETURN
  102 DO 103 K=1,L1
         M2 = M2S
         DO 103 M1=1,M1D,IM1
         M2 = M2+IM2
         TI5 = CC(2,M1,K,1,2)-CC(2,M1,K,1,5)
         TI2 = CC(2,M1,K,1,2)+CC(2,M1,K,1,5)
         TI4 = CC(2,M1,K,1,3)-CC(2,M1,K,1,4)
         TI3 = CC(2,M1,K,1,3)+CC(2,M1,K,1,4)
         TR5 = CC(1,M1,K,1,2)-CC(1,M1,K,1,5)
         TR2 = CC(1,M1,K,1,2)+CC(1,M1,K,1,5)
         TR4 = CC(1,M1,K,1,3)-CC(1,M1,K,1,4)
         TR3 = CC(1,M1,K,1,3)+CC(1,M1,K,1,4)
         CH(1,M2,K,1,1) = CC(1,M1,K,1,1)+TR2+TR3
         CH(2,M2,K,1,1) = CC(2,M1,K,1,1)+TI2+TI3
         CR2 = CC(1,M1,K,1,1)+TR11*TR2+TR12*TR3
         CI2 = CC(2,M1,K,1,1)+TR11*TI2+TR12*TI3
         CR3 = CC(1,M1,K,1,1)+TR12*TR2+TR11*TR3
         CI3 = CC(2,M1,K,1,1)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,M2,K,2,1) = CR2-CI5
         CH(1,M2,K,5,1) = CR2+CI5
         CH(2,M2,K,2,1) = CI2+CR5
         CH(2,M2,K,3,1) = CI3+CR4
         CH(1,M2,K,3,1) = CR3-CI4
         CH(1,M2,K,4,1) = CR3+CI4
         CH(2,M2,K,4,1) = CI3-CR4
         CH(2,M2,K,5,1) = CI2-CR5
  103 CONTINUE
      DO 105 I=2,IDO
         DO 104 K=1,L1
         M2 = M2S
         DO 104 M1=1,M1D,IM1
         M2 = M2+IM2
            TI5 = CC(2,M1,K,I,2)-CC(2,M1,K,I,5)
            TI2 = CC(2,M1,K,I,2)+CC(2,M1,K,I,5)
            TI4 = CC(2,M1,K,I,3)-CC(2,M1,K,I,4)
            TI3 = CC(2,M1,K,I,3)+CC(2,M1,K,I,4)
            TR5 = CC(1,M1,K,I,2)-CC(1,M1,K,I,5)
            TR2 = CC(1,M1,K,I,2)+CC(1,M1,K,I,5)
            TR4 = CC(1,M1,K,I,3)-CC(1,M1,K,I,4)
            TR3 = CC(1,M1,K,I,3)+CC(1,M1,K,I,4)
            CH(1,M2,K,1,I) = CC(1,M1,K,I,1)+TR2+TR3
            CH(2,M2,K,1,I) = CC(2,M1,K,I,1)+TI2+TI3
            CR2 = CC(1,M1,K,I,1)+TR11*TR2+TR12*TR3
            CI2 = CC(2,M1,K,I,1)+TR11*TI2+TR12*TI3
            CR3 = CC(1,M1,K,I,1)+TR12*TR2+TR11*TR3
            CI3 = CC(2,M1,K,I,1)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(1,M2,K,2,I) = WA(I,1,1)*DR2+WA(I,1,2)*DI2
            CH(2,M2,K,2,I) = WA(I,1,1)*DI2-WA(I,1,2)*DR2
            CH(1,M2,K,3,I) = WA(I,2,1)*DR3+WA(I,2,2)*DI3
            CH(2,M2,K,3,I) = WA(I,2,1)*DI3-WA(I,2,2)*DR3
            CH(1,M2,K,4,I) = WA(I,3,1)*DR4+WA(I,3,2)*DI4
            CH(2,M2,K,4,I) = WA(I,3,1)*DI4-WA(I,3,2)*DR4
            CH(1,M2,K,5,I) = WA(I,4,1)*DR5+WA(I,4,2)*DI5
            CH(2,M2,K,5,I) = WA(I,4,1)*DI5-WA(I,4,2)*DR5
  104    CONTINUE
  105 CONTINUE
      END subroutine CMF5KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMFGKB (LOT,IDO,IP,L1,LID,NA,CC,CC1,IM1,IN1,CH,CH1,IM2,IN2,WA)
      integer           :: LOT,LID,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDLJ
      real(kind=rk)     :: WAR,WAI    
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CH(2,IN2,L1,IDO,IP),CC(2,IN1,L1,IP,IDO),CC1(2,IN1,LID,IP),CH1(2,IN2,LID,IP),WA(IDO,IP-1,2)
!
! FFTPACK 5.0 auxiliary routine
!
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      DO 110 KI=1,LID
         M2 = M2S
         DO 110 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = CC1(2,M1,KI,1)
  110 CONTINUE
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 112 KI=1,LID
         M2 = M2S
         DO 112 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) =  CC1(1,M1,KI,J)+CC1(1,M1,KI,JC)
            CH1(1,M2,KI,JC) = CC1(1,M1,KI,J)-CC1(1,M1,KI,JC)
            CH1(2,M2,KI,J) =  CC1(2,M1,KI,J)+CC1(2,M1,KI,JC)
            CH1(2,M2,KI,JC) = CC1(2,M1,KI,J)-CC1(2,M1,KI,JC)
  112    CONTINUE
  111 CONTINUE
      DO 118 J=2,IPPH
         DO 117 KI=1,LID
         M2 = M2S
         DO 117 M1=1,M1D,IM1
         M2 = M2+IM2
            CC1(1,M1,KI,1) = CC1(1,M1,KI,1)+CH1(1,M2,KI,J)
            CC1(2,M1,KI,1) = CC1(2,M1,KI,1)+CH1(2,M2,KI,J)
  117    CONTINUE
  118 CONTINUE
      DO 116 L=2,IPPH
         LC = IPP2-L
         DO 113 KI=1,LID
         M2 = M2S
         DO 113 M1=1,M1D,IM1
         M2 = M2+IM2
            CC1(1,M1,KI,L) = CH1(1,M2,KI,1)+WA(1,L-1,1)*CH1(1,M2,KI,2)
            CC1(1,M1,KI,LC) = WA(1,L-1,2)*CH1(1,M2,KI,IP)
            CC1(2,M1,KI,L) = CH1(2,M2,KI,1)+WA(1,L-1,1)*CH1(2,M2,KI,2)
            CC1(2,M1,KI,LC) = WA(1,L-1,2)*CH1(2,M2,KI,IP)
  113    CONTINUE
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = MOD((L-1)*(J-1),IP)
            WAR = WA(1,IDLJ,1)
            WAI = WA(1,IDLJ,2)
            DO 114 KI=1,LID
               M2 = M2S
               DO 114 M1=1,M1D,IM1
               M2 = M2+IM2
               CC1(1,M1,KI,L) = CC1(1,M1,KI,L)+WAR*CH1(1,M2,KI,J)
               CC1(1,M1,KI,LC) = CC1(1,M1,KI,LC)+WAI*CH1(1,M2,KI,JC)
               CC1(2,M1,KI,L) = CC1(2,M1,KI,L)+WAR*CH1(2,M2,KI,J)
               CC1(2,M1,KI,LC) = CC1(2,M1,KI,LC)+WAI*CH1(2,M2,KI,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      IF(IDO .GT. 1 .OR. NA .EQ. 1) GO TO 136
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 KI=1,LID
         DO 119 M1=1,M1D,IM1
            CHOLD1 = CC1(1,M1,KI,J)-CC1(2,M1,KI,JC)
            CHOLD2 = CC1(1,M1,KI,J)+CC1(2,M1,KI,JC)
            CC1(1,M1,KI,J) = CHOLD1
            CC1(2,M1,KI,JC) = CC1(2,M1,KI,J)-CC1(1,M1,KI,JC)
            CC1(2,M1,KI,J) = CC1(2,M1,KI,J)+CC1(1,M1,KI,JC)
            CC1(1,M1,KI,JC) = CHOLD2
  119    CONTINUE
  120 CONTINUE
      RETURN
  136 DO 137 KI=1,LID
         M2 = M2S
         DO 137 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = CC1(2,M1,KI,1)
  137 CONTINUE
      DO 135 J=2,IPPH
         JC = IPP2-J
         DO 134 KI=1,LID
         M2 = M2S
         DO 134 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) = CC1(1,M1,KI,J)-CC1(2,M1,KI,JC)
            CH1(1,M2,KI,JC) = CC1(1,M1,KI,J)+CC1(2,M1,KI,JC)
            CH1(2,M2,KI,JC) = CC1(2,M1,KI,J)-CC1(1,M1,KI,JC)
            CH1(2,M2,KI,J) = CC1(2,M1,KI,J)+CC1(1,M1,KI,JC)
  134    CONTINUE
  135 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 131 I=1,IDO
         DO 130 K=1,L1
         M2 = M2S
         DO 130 M1=1,M1D,IM1
         M2 = M2+IM2
            CC(1,M1,K,1,I) = CH(1,M2,K,I,1)
            CC(2,M1,K,1,I) = CH(2,M2,K,I,1)
  130    CONTINUE
  131 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
         M2 = M2S
         DO 122 M1=1,M1D,IM1
         M2 = M2+IM2
            CC(1,M1,K,J,1) = CH(1,M2,K,1,J)
            CC(2,M1,K,J,1) = CH(2,M2,K,1,J)
  122    CONTINUE
  123 CONTINUE
      DO 126 J=2,IP
         DO 125 I=2,IDO
            DO 124 K=1,L1
               M2 = M2S
               DO 124 M1=1,M1D,IM1
               M2 = M2+IM2
               CC(1,M1,K,J,I) = WA(I,J-1,1)*CH(1,M2,K,I,J)-WA(I,J-1,2)*CH(2,M2,K,I,J)
               CC(2,M1,K,J,I) = WA(I,J-1,1)*CH(2,M2,K,I,J)+WA(I,J-1,2)*CH(1,M2,K,I,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      END subroutine CMFGKB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE CMFGKF (LOT,IDO,IP,L1,LID,NA,CC,CC1,IM1,IN1,CH,CH1,IM2,IN2,WA)
      integer           :: LOT,LID,IN1,IN2,IM1,IM2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDLJ
      real(kind=rk)     :: WAR,WAI    
      real(kind=rk)     :: SN
      real(kind=rk)     :: CHOLD1,CHOLD2
      REAL(kind=rk)     :: CH(2,IN2,L1,IDO,IP),CC(2,IN1,L1,IP,IDO),CC1(2,IN1,LID,IP),CH1(2,IN2,LID,IP),WA(IDO,IP-1,2)
!
! FFTPACK 5.0 auxiliary routine
!
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      DO 110 KI=1,LID
         M2 = M2S
         DO 110 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = CC1(2,M1,KI,1)
  110 CONTINUE
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 112 KI=1,LID
         M2 = M2S
         DO 112 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) =  CC1(1,M1,KI,J)+CC1(1,M1,KI,JC)
            CH1(1,M2,KI,JC) = CC1(1,M1,KI,J)-CC1(1,M1,KI,JC)
            CH1(2,M2,KI,J) =  CC1(2,M1,KI,J)+CC1(2,M1,KI,JC)
            CH1(2,M2,KI,JC) = CC1(2,M1,KI,J)-CC1(2,M1,KI,JC)
  112    CONTINUE
  111 CONTINUE
      DO 118 J=2,IPPH
         DO 117 KI=1,LID
         M2 = M2S
         DO 117 M1=1,M1D,IM1
         M2 = M2+IM2
            CC1(1,M1,KI,1) = CC1(1,M1,KI,1)+CH1(1,M2,KI,J)
            CC1(2,M1,KI,1) = CC1(2,M1,KI,1)+CH1(2,M2,KI,J)
  117    CONTINUE
  118 CONTINUE
      DO 116 L=2,IPPH
         LC = IPP2-L
         DO 113 KI=1,LID
         M2 = M2S
         DO 113 M1=1,M1D,IM1
         M2 = M2+IM2
            CC1(1,M1,KI,L) = CH1(1,M2,KI,1)+WA(1,L-1,1)*CH1(1,M2,KI,2)
            CC1(1,M1,KI,LC) = -WA(1,L-1,2)*CH1(1,M2,KI,IP)
            CC1(2,M1,KI,L) = CH1(2,M2,KI,1)+WA(1,L-1,1)*CH1(2,M2,KI,2)
            CC1(2,M1,KI,LC) = -WA(1,L-1,2)*CH1(2,M2,KI,IP)
  113    CONTINUE
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = MOD((L-1)*(J-1),IP)
            WAR = WA(1,IDLJ,1)
            WAI = -WA(1,IDLJ,2)
            DO 114 KI=1,LID
               M2 = M2S
               DO 114 M1=1,M1D,IM1
               M2 = M2+IM2
               CC1(1,M1,KI,L) = CC1(1,M1,KI,L)+WAR*CH1(1,M2,KI,J)
               CC1(1,M1,KI,LC) = CC1(1,M1,KI,LC)+WAI*CH1(1,M2,KI,JC)
               CC1(2,M1,KI,L) = CC1(2,M1,KI,L)+WAR*CH1(2,M2,KI,J)
               CC1(2,M1,KI,LC) = CC1(2,M1,KI,LC)+WAI*CH1(2,M2,KI,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      IF (IDO .GT. 1) GO TO 136
      SN = 1._rk/REAL(IP*L1,kind=rk)
      IF (NA .EQ. 1) GO TO 146
      DO 149 KI=1,LID
         M2 = M2S
         DO 149 M1=1,M1D,IM1
         M2 = M2+IM2
         CC1(1,M1,KI,1) = SN*CC1(1,M1,KI,1)
         CC1(2,M1,KI,1) = SN*CC1(2,M1,KI,1)
  149 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 KI=1,LID
         DO 119 M1=1,M1D,IM1
            CHOLD1 = SN*(CC1(1,M1,KI,J)-CC1(2,M1,KI,JC))
            CHOLD2 = SN*(CC1(1,M1,KI,J)+CC1(2,M1,KI,JC))
            CC1(1,M1,KI,J) = CHOLD1
            CC1(2,M1,KI,JC) = SN*(CC1(2,M1,KI,J)-CC1(1,M1,KI,JC))
            CC1(2,M1,KI,J) = SN*(CC1(2,M1,KI,J)+CC1(1,M1,KI,JC))
            CC1(1,M1,KI,JC) = CHOLD2
  119    CONTINUE
  120 CONTINUE
      RETURN
  146 DO 147 KI=1,LID
         M2 = M2S
         DO 147 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = SN*CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = SN*CC1(2,M1,KI,1)
  147 CONTINUE
      DO 145 J=2,IPPH
         JC = IPP2-J
         DO 144 KI=1,LID
         M2 = M2S
         DO 144 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) = SN*(CC1(1,M1,KI,J)-CC1(2,M1,KI,JC))
            CH1(2,M2,KI,J) = SN*(CC1(2,M1,KI,J)+CC1(1,M1,KI,JC))
            CH1(1,M2,KI,JC) = SN*(CC1(1,M1,KI,J)+CC1(2,M1,KI,JC))
            CH1(2,M2,KI,JC) = SN*(CC1(2,M1,KI,J)-CC1(1,M1,KI,JC))
  144    CONTINUE
  145 CONTINUE
      RETURN
  136 DO 137 KI=1,LID
         M2 = M2S
         DO 137 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = CC1(2,M1,KI,1)
  137 CONTINUE
      DO 135 J=2,IPPH
         JC = IPP2-J
         DO 134 KI=1,LID
         M2 = M2S
         DO 134 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) = CC1(1,M1,KI,J)-CC1(2,M1,KI,JC)
            CH1(2,M2,KI,J) = CC1(2,M1,KI,J)+CC1(1,M1,KI,JC)
            CH1(1,M2,KI,JC) = CC1(1,M1,KI,J)+CC1(2,M1,KI,JC)
            CH1(2,M2,KI,JC) = CC1(2,M1,KI,J)-CC1(1,M1,KI,JC)
  134    CONTINUE
  135 CONTINUE
      DO 131 I=1,IDO
         DO 130 K=1,L1
         M2 = M2S
         DO 130 M1=1,M1D,IM1
         M2 = M2+IM2
            CC(1,M1,K,1,I) = CH(1,M2,K,I,1)
            CC(2,M1,K,1,I) = CH(2,M2,K,I,1)
  130    CONTINUE
  131 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
         M2 = M2S
         DO 122 M1=1,M1D,IM1
         M2 = M2+IM2
            CC(1,M1,K,J,1) = CH(1,M2,K,1,J)
            CC(2,M1,K,J,1) = CH(2,M2,K,1,J)
  122    CONTINUE
  123 CONTINUE
      DO 126 J=2,IP
         DO 125 I=2,IDO
            DO 124 K=1,L1
               M2 = M2S
               DO 124 M1=1,M1D,IM1
               M2 = M2+IM2
               CC(1,M1,K,J,I) = WA(I,J-1,1)*CH(1,M2,K,I,J)+WA(I,J-1,2)*CH(2,M2,K,I,J)
               CC(2,M1,K,J,I) = WA(I,J-1,1)*CH(2,M2,K,I,J)-WA(I,J-1,2)*CH(1,M2,K,I,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      END subroutine CMFGKF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CMFM1B (LOT,JUMP,N,INC,C,CH,WA,FNF,FAC)
      SUBROUTINE CMFM1B (LOT,JUMP,N,INC,C,CH,fft_sign)
      integer           :: LOT,JUMP,N,INC
      integer           :: NF,NA,LID,NBR
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      !Kuester COMPLEX(kind=rk)   :: C(*)
      real(kind=rk)     :: fnf  
      real(kind=rk)     :: C(*)
      REAL(kind=rk)     :: CH(*)
      !!REAL(kind=rk)     :: WA(*), FAC(*)
      type(fft_real_sign_type) :: fft_sign
!
! FFTPACK 5.0 auxiliary routine
!
      !!NF = FNF
      NF = fft_sign%NF
      NA = 0
      L1 = 1
      IW = 1
      DO 125 K1=1,fft_sign%NF
         IP = fft_sign%FAC(K1)
         L2 = IP*L1
         IDO = N/L2
         LID = L1*IDO
         NBR = 1+NA+2*MIN(IP-2,4)
         GO TO (52,62,53,63,54,64,55,65,56,66),NBR
   52    CALL CMF2KB (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   62    CALL CMF2KB (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   53    CALL CMF3KB (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   63    CALL CMF3KB (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   54    CALL CMF4KB (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   64    CALL CMF4KB (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   55    CALL CMF5KB (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   65    CALL CMF5KB (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   56    CALL CMFGKB (LOT,IDO,IP,L1,LID,NA,C,C,JUMP,INC,CH,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   66    CALL CMFGKB (LOT,IDO,IP,L1,LID,NA,CH,CH,1,LOT,C,C,JUMP,INC,fft_sign%WA(IW))
  120    L1 = L2
         IW = IW+(IP-1)*(IDO+IDO)
         IF(IP .LE. 5) NA = 1-NA
  125 CONTINUE
      END subroutine CMFM1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE CMFM1F (LOT,JUMP,N,INC,C,CH,WA,FNF,FAC)
      SUBROUTINE CMFM1F (LOT,JUMP,N,INC,C,CH,fft_sign)   
      integer           :: LOT,JUMP,N,INC
      integer           :: NF,NA,LID,NBR
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      real(kind=rk)     :: fnf  
      !Kuester COMPLEX(kind=rk)   :: C(*)
      real(kind=rk)     :: C(*)
      REAL(kind=rk)     :: CH(*)
      !!REAL(kind=rk)     :: WA(*), FAC(*)
      type(fft_real_sign_type) :: fft_sign
!
! FFTPACK 5.0 auxiliary routine
!
      !!NF = FNF
      NF = fft_sign%NF
      NA = 0
      L1 = 1
      IW = 1
      DO 125 K1=1,fft_sign%NF
         IP = fft_sign%FAC(K1)
         L2 = IP*L1
         IDO = N/L2
         LID = L1*IDO
         NBR = 1+NA+2*MIN(IP-2,4)
   !!      write(*,'(*(a,i0))') 'CMFM1F: K1 =',K1,' IP =',IP,' L2 =',L2,' IDO =',IDO,' NPR=',NBR &
   !!                           & ,' LOT ',LOT,' JUMP ',JUMP,' INC ',INC 
   !!      if(INC < (LOT-1)*JUMP+1) then
!test!   !!              write(*,*) 'Warning, INC < (LOT-1)*JUMP+1'  !test_print
   !!      endif
         GO TO (52,62,53,63,54,64,55,65,56,66),NBR
   52    CALL CMF2KF (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   62    CALL CMF2KF (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   53    CALL CMF3KF (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   63    CALL CMF3KF (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   54    CALL CMF4KF (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   64    CALL CMF4KF (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   55    CALL CMF5KF (LOT,IDO,L1,NA,C,JUMP,INC,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   65    CALL CMF5KF (LOT,IDO,L1,NA,CH,1,LOT,C,JUMP,INC,fft_sign%WA(IW))
         GO TO 120
   56    CALL CMFGKF (LOT,IDO,IP,L1,LID,NA,C,C,JUMP,INC,CH,CH,1,LOT,fft_sign%WA(IW))
         GO TO 120
   66    CALL CMFGKF (LOT,IDO,IP,L1,LID,NA,CH,CH,1,LOT,C,C,JUMP,INC,fft_sign%WA(IW))
  120    L1 = L2
         IW = IW+(IP-1)*(IDO+IDO)
         IF(IP .LE. 5) NA = 1-NA
  125 CONTINUE
      END subroutine CMFM1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQ1B (N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE COSQ1B (N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE COSQ1B (N, INC, X, LENX, fft_sign, IER)
      integer     :: N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IER1
      real(kind=rk)     :: X1
      real(kind=rk)     :: SSQRT2
    !!  REAL(kind=rk)     :: WSAVE(LENSAV)
    !!  real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSQ1B', 6)
        GO TO 300
    !!  ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
    !!    IER = 2
    !!    CALL XERFFT ('COSQ1B', 8)
    !!    GO TO 300
    !!  ELSEIF (LENWRK .LT. N) THEN
    !!    IER = 3
    !!    CALL XERFFT ('COSQ1B', 10)
    !!    GO TO 300
      ENDIF
!
      IF (N-2) 300,102,103
 102  SSQRT2 = 1._rk/SQRT(2._rk)
      X1 = X(1,1)+X(1,2)
      X(1,2) = SSQRT2*(X(1,1)-X(1,2))
      X(1,1) = X1
      RETURN
  103 continue
      !!CALL COSQB1 (N,INC,X,WSAVE,WORK, fft_sign,IER1)
      CALL COSQB1 (N,INC,X,fft_sign%WORK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQ1B',-5)
      ENDIF
!
  300 continue
      END subroutine COSQ1B
!************************************************************************************************************************
!!      SUBROUTINE COSQ1B_work_alloc (N, INC, X, LENX, fft_sign, IER)
!!      integer     :: N, INC, LENX
!!      integer    ::  LENWRK, IER
!!      !!integer   :: LENSAV
!!      integer           :: IER1
!!      real(kind=rk)     :: X1
!!      real(kind=rk)     :: SSQRT2
!!      real(kind=rk),allocatable,dimension(:) :: WORK
!!      real(kind=rk)     :: X(INC,*)
!!      !!real(kind=rk)     :: X(*)
!!      type(fft_real_sign_type) :: fft_sign
!!!
!!      IER = 0
!!!
!!      IF (LENX .LT. INC*(N-1) + 1) THEN
!!        IER = 1
!!        CALL XERFFT ('COSQ1B', 6)
!!        GO TO 300
!!      ENDIF
!!!
!!      LENWRK=N
!!      allocate(WORK(LENWRK))
!!
!!      IF (N-2) 300,102,103
!! 102  SSQRT2 = 1._rk/SQRT(2._rk)
!!      X1 = X(1,1)+X(1,2)
!!      X(1,2) = SSQRT2*(X(1,1)-X(1,2))
!!      X(1,1) = X1
!!      RETURN
!!  103 continue
!!      !!CALL COSQB1 (N,INC,X,WSAVE,WORK, fft_sign,IER1)
!!      CALL COSQB1 (N,INC,X,WORK, fft_sign,IER1)
!!      IF (IER1 .NE. 0) THEN
!!        IER = 20
!!        CALL XERFFT ('COSQ1B',-5)
!!      ENDIF
!!!
!!  300 continue
!!      END subroutine COSQ1B_work_alloc
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQ1F (N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE COSQ1F (N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE COSQ1F (N, INC, X, LENX, fft_sign, IER)
      integer     :: N, INC, LENX
      !!integer    ::  LENWRK
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IER1
      real(kind=rk)     :: SSQRT2
      real(kind=rk)     :: TSQX
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSQ1F', 6)
        GO TO 300
   !!   ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
   !!     IER = 2
   !!     CALL XERFFT ('COSQ1F', 8)
   !!     GO TO 300
   !!   ELSEIF (LENWRK .LT. N) THEN
   !!     IER = 3
   !!     CALL XERFFT ('COSQ1F', 10)
   !!     GO TO 300
      ENDIF
!
      IF (N-2) 102,101,103
  101 SSQRT2 = 1._rk/SQRT(2._rk)
      TSQX = SSQRT2*X(1,2)
      X(1,2) = .5_rk*X(1,1)-TSQX
      X(1,1) = .5_rk*X(1,1)+TSQX
  102 RETURN
  103 continue
      !!CALL COSQF1 (N,INC,X,WSAVE,WORK,fft_sign, IER1)
      CALL COSQF1 (N,INC,X,fft_sign%WORK,fft_sign, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQ1F',-5)
      ENDIF
!
  300 CONTINUE
      END subroutine COSQ1F
!************************************************************************************************************************

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQ1I (N, WSAVE, LENSAV,fft_sign, IER)
      SUBROUTINE COSQ1I (N, fft_sign, IER)
      integer     :: N
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      real(kind=rk)     :: PIH,DT,FK
     !! REAL(kind=rk)     :: WSAVE(LENSAV)
      integer     :: LENWRK
      type(fft_real_sign_type) :: fft_sign
      integer                  :: size_WA
      integer                  :: NF
!
if( debug > 0) write(*,*) 'COSQ1I start '
      if( fft_sign%switch == 2 .and. fft_sign%N == N .and. allocated(fft_sign%FAC) .and. allocated(fft_sign%WA)) then
           goto 1111
      endif

      fft_sign%switch = 2

      fft_sign%N = N
      IER = 0
  !!    IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
  !!      IER = 2
  !!      CALL XERFFT ('COSQ1I', 3)
  !!      GO TO 300
  !!    ENDIF
!
!!      LNSV = N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4


      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(N))


      !!PIH = 2._rk*ATAN(1._rk)
      PIH = 0.5_dk*PI_long
      DT = PIH/REAL(N,kind=rk)
      FK = 0._rk
      DO 101 K=1,N
         FK = FK+1._rk
        !! WSAVE(K) = COS(FK*DT)
         fft_sign%WSAVE(K)=COS(FK*DT)
  101 CONTINUE

      !!CALL RFFT1I_orig (N, WSAVE(N+1), LNSV, IER1)
      CALL all_FFT1I(N, fft_sign)
      fft_sign%sign_type=COSQ1_sign_type

      LENWRK=N
      if(allocated(fft_sign%work) ) deallocate(fft_sign%work)
      allocate(fft_sign%work(LENWRK))
      IER=0


  300 CONTINUE
  1111 continue
  if( debug > 0 ) write(*,*) 'COSQ1I end '
      END subroutine COSQ1I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQB1 (N,INC,X,WSAVE,WORK, fft_sign,IER)
      SUBROUTINE COSQB1 (N,INC,X,WORK, fft_sign,IER)
      integer       :: N,INC,IER
      integer       :: NS2,NP2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IER1
      real(kind=rk)     :: XIM1 
      !!real(kind=rk) :: WSAVE(*)
      real(kind=rk)     :: WORK(*)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign


      IER = 0
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 I=3,N,2
         XIM1 = X(1,I-1)+X(1,I)
         X(1,I) = .5_rk*(X(1,I-1)-X(1,I))
         X(1,I-1) = .5_rk*XIM1
  101 CONTINUE
      X(1,1) = .5_rk*X(1,1)
      MODN = MOD(N,2)
      IF (MODN .NE. 0) GO TO 302
      X(1,N) = .5_rk*X(1,N)
  302 LENX = INC*(N-1)  + 1
      LNSV = N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
      LNWK = N
!
      !!CALL RFFT1B(N,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK, fft_sign,IER1)
      !!CALL RFFT1B(N,INC,X,LENX,WORK,LNWK, fft_sign,IER1)
      CALL RFFT1B(N,INC,X,LENX,fft_sign,IER1)

      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQB1',-5)
        GO TO 400
      ENDIF
!
      DO 102 K=2,NS2
         KC = NP2-K
         !!WORK(K) = WSAVE(K-1)*X(1,KC)+WSAVE(KC-1)*X(1,K)
         !!WORK(KC) = WSAVE(K-1)*X(1,K)-WSAVE(KC-1)*X(1,KC)
         WORK(K)  = fft_sign%WSAVE(K-1)*X(1,KC)+fft_sign%WSAVE(KC-1)*X(1,K)
         WORK(KC) = fft_sign%WSAVE(K-1)*X(1,K )-fft_sign%WSAVE(KC-1)*X(1,KC)
  102 CONTINUE
      IF (MODN .NE. 0) GO TO 305
      X(1,NS2+1) = fft_sign%WSAVE(NS2)*(X(1,NS2+1)+X(1,NS2+1))
  305 DO 103 K=2,NS2
         KC = NP2-K
         X(1,K) = WORK(K)+WORK(KC)
         X(1,KC) = WORK(K)-WORK(KC)
  103 CONTINUE
      X(1,1) = X(1,1)+X(1,1)
  400 continue
      END subroutine COSQB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQF1 (N,INC,X,WSAVE,WORK,fft_sign,IER)
      SUBROUTINE COSQF1 (N,INC,X,WORK,fft_sign,IER)
      integer       :: N,INC,IER
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: XIM1 
      !!real(kind=rk) :: WSAVE(*)
      real(kind=rk)     :: WORK(*)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign

      if( debug > 0 ) write(*,*) 'COSQF1 start '
      IER = 0
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 K=2,NS2
         KC = NP2-K
         WORK(K)  = X(1,K)+X(1,KC)
         WORK(KC) = X(1,K)-X(1,KC)
  101 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) GO TO 301
      WORK(NS2+1) = X(1,NS2+1)+X(1,NS2+1)
  301 DO 102 K=2,NS2
         KC = NP2-K
         X(1,K)  = fft_sign%WSAVE(K-1)*WORK(KC)+fft_sign%WSAVE(KC-1)*WORK(K)
         X(1,KC) = fft_sign%WSAVE(K-1)*WORK(K) -fft_sign%WSAVE(KC-1)*WORK(KC)
  102 CONTINUE
      IF (MODN .NE. 0) GO TO 303
      X(1,NS2+1) = fft_sign%WSAVE(NS2)*WORK(NS2+1)
  303 LENX = INC*(N-1)  + 1
      LNSV = N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
      LNWK = N
!
      !!CALL RFFT1F(N,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK,fft_sign,IER1)
      !!CALL RFFT1F(N,INC,X,LENX,WORK,LNWK,fft_sign,IER1)
      CALL RFFT1F(N,INC,X,LENX,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQF1',-5)
        GO TO 400
      ENDIF
!
      DO 103 I=3,N,2
         XIM1 = .5_rk*(X(1,I-1)+X(1,I))
         X(1,I) = .5_rk*(X(1,I-1)-X(1,I))
         X(1,I-1) = XIM1
  103 CONTINUE
  400 continue
  if( debug > 0 ) write(*,*) 'COSQF1 end'
      END subroutine COSQF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQMB (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      SUBROUTINE COSQMB (LOT, JUMP, N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: X1
      real(kind=rk)     :: SSQRT2
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSQMB', 6)
        GO TO 300
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('COSQMB', 8)
!!        GO TO 300
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('COSQMB', 10)
        GO TO 300
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('COSQMB', -1)
        GO TO 300
      ENDIF
!
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 101,102,103
 101  DO 201 M=1,LJ,JUMP
      X(M,1) = X(M,1)
 201  CONTINUE
      RETURN
 102  SSQRT2 = 1._rk/SQRT(2._rk)
      DO 202 M=1,LJ,JUMP
      X1 = X(M,1)+X(M,2)
      X(M,2) = SSQRT2*(X(M,1)-X(M,2))
      X(M,1) = X1
 202  CONTINUE
      RETURN
  103 continue 
      !!CALL MCSQB1 (LOT,JUMP,N,INC,X,WSAVE,WORK,fft_sign,IER1)
      CALL MCSQB1 (LOT,JUMP,N,INC,X,WORK,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQMB',-5)
      ENDIF
!
  300 CONTINUE
      END subroutine COSQMB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQMF (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      SUBROUTINE COSQMF (LOT, JUMP, N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: SSQRT2
      real(kind=rk)     :: TSQX
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSQMF', 6)
        GO TO 300
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('COSQMF', 8)
!!        GO TO 300
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('COSQMF', 10)
        GO TO 300
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('COSQMF', -1)
        GO TO 300
      ENDIF
!
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 102,101,103
  101 SSQRT2 = 1._rk/SQRT(2._rk)
      DO 201 M=1,LJ,JUMP
      TSQX = SSQRT2*X(M,2)
      X(M,2) = .5_rk*X(M,1)-TSQX
      X(M,1) = .5_rk*X(M,1)+TSQX
  201 CONTINUE
  102 RETURN
  103 continue
      !!CALL MCSQF1 (LOT,JUMP,N,INC,X,WSAVE,WORK, fft_sign, IER1)
      CALL MCSQF1 (LOT,JUMP,N,INC,X,WORK, fft_sign, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQMF',-5)
      ENDIF
!
  300 CONTINUE
      END subroutine COSQMF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSQMI (N, WSAVE, LENSAV, fft_sign, IER)
      SUBROUTINE COSQMI (N, fft_sign, IER)
      integer     :: N
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: PIH,DT,FK
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
!!      IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('COSQMI', 3)
!!        GO TO 300
!!      ENDIF
!
      PIH = 2._rk*ATAN(1._rk)
      PIH = 0.5_dk*PI_long
      DT = PIH/REAL(N,kind=rk)

      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(N))

      FK = 0._rk
      DO 101 K=1,N
         FK = FK+1._rk
         !!WSAVE(K) = COS(FK*DT)
         fft_sign%WSAVE(K) = COS(FK*DT)
  101 CONTINUE

      !!LNSV = N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4
      !!CALL RFFTMI (N, WSAVE(N+1), LNSV,fft_sign, IER1)
      CALL RFFTMI (N, fft_sign, IER1)

      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQMI',-5)
      ENDIF
  300 CONTINUE
      END subroutine COSQMI
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COST1B ( N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE COST1B ( N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE COST1B ( N, INC, X, LENX, fft_sign, IER)
      integer     :: N, INC, LENX
      !!integer    ::  LENWRK
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IER1
     !! REAL(kind=rk)     :: WSAVE(LENSAV)
     !! real(kind=rk)     :: WORK(LENWRK)
      !!real(kind=rk)     :: X(INC,*)
      real(kind=rk)     :: X(*)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COST1B', 6)
        GO TO 100
  !!    ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
  !!      IER = 2
  !!      CALL XERFFT ('COST1B', 8)
  !!      GO TO 100
  !!    ELSEIF (LENWRK .LT. N-1) THEN
  !!      IER = 3
  !!      CALL XERFFT ('COST1B', 10)
  !!      GO TO 100
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!CALL COSTB1 (N,INC,X,WSAVE,WORK,fft_sign,IER1)
      !!CALL COSTB1 (N,INC,X,WORK,fft_sign,IER1)
      CALL COSTB1 (N,INC,X,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COST1B',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine COST1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COST1F ( N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE COST1F ( N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE COST1F ( N, INC, X, LENX, fft_sign, IER)
      integer     :: N, INC, LENX
      !!integer    ::  LENWRK
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IER1
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COST1F', 6)
        GO TO 100
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('COST1F', 8)
!!        GO TO 100
  !!    ELSEIF (LENWRK .LT. N-1) THEN
  !!      IER = 3
  !!      CALL XERFFT ('COST1F', 10)
  !!      GO TO 100
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!CALL COSTF1(N,INC,X,WSAVE,WORK, fft_sign,IER1)
      !!CALL COSTF1(N,INC,X,WORK, fft_sign,IER1)
      CALL COSTF1(N,INC,X,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COST1F',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine COST1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COST1I (N, WSAVE, LENSAV,fft_sign, IER)
      SUBROUTINE COST1I (N, fft_sign, IER)
      integer     :: N
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: PI,DT,FK
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      integer     :: LENWRK
      type(fft_real_sign_type) :: fft_sign
      integer                  :: size_WA
!
  !!    if( fft_sign%switch == 2 .and. fft_sign%N == N .and. allocated(fft_sign%FAC) .and. allocated(fft_sign%WA)) then
  !!         goto 1111
  !!    endif

      fft_sign%switch = 2

      fft_sign%N = N

      IER = 0
!
  !!    IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
  !!      IER = 2
  !!      CALL XERFFT ('COST1I', 3)
  !!      GO TO 300
  !!    ENDIF
!
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
     !! LNSV = NM1 + INT(LOG(REAL(NM1))/LOG(2._rk)) +4


      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(NM1))



      !!PI = 4._rk*ATAN(1._rk)
      PI = PI_long
      DT = PI/REAL(NM1,kind=rk)
      FK = 0._rk
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1._rk
       !!  WSAVE(K) = 2._rk*SIN(FK*DT)
       !!  WSAVE(KC) = 2._rk*COS(FK*DT)
         fft_sign%WSAVE(K) = 2._rk*SIN(FK*DT)
         fft_sign%WSAVE(KC) = 2._rk*COS(FK*DT)
  101 CONTINUE

      !!CALL RFFT1I_orig (NM1, WSAVE(N+1), LNSV, IER1)
      CALL all_FFT1I(NM1, fft_sign)
      fft_sign%sign_type=COST1_sign_type

      LENWRK=N-1
      if(allocated(fft_sign%work) ) deallocate(fft_sign%work)
      allocate(fft_sign%work(LENWRK))

      IER=0

  300 CONTINUE
  1111 continue
      END subroutine COST1I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSTB1(N,INC,X,WSAVE,WORK,fft_sign,IER)
      !!SUBROUTINE COSTB1(N,INC,X,WORK,fft_sign,IER)
      SUBROUTINE COSTB1(N,INC,X,fft_sign,IER)
      INTEGER           :: N,INC,IER
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      !!REAL(kind=rk)     :: WSAVE(*)
      REAL(kind=rk)     :: X(INC,*)
      !!real(kind=rk)     :: WORK(*)  ! Kuester
      real(kind=rk)     :: T1,T2
      real(kind=rk)     :: XI
      real(kind=rk)     :: X1H
      real(kind=rk)     :: X1P3,X2
      real(kind=rk)     :: FNM1S2,FNM1S4
      real(kind=dk)    :: DSUM
      type(fft_real_sign_type) :: fft_sign
      IER = 0
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      IF (N-2) 106,101,102
  101 X1H = X(1,1)+X(1,2)
      X(1,2) = X(1,1)-X(1,2)
      X(1,1) = X1H
      RETURN
  102 IF (N .GT. 3) GO TO 103
      X1P3 = X(1,1)+X(1,3)
      X2 = X(1,2)
      X(1,2) = X(1,1)-X(1,3)
      X(1,1) = X1P3+X2
      X(1,3) = X1P3-X2
      RETURN
  103 X(1,1) = X(1,1)+X(1,1)
      X(1,N) = X(1,N)+X(1,N)
      DSUM = X(1,1)-X(1,N)
      X(1,1) = X(1,1)+X(1,N)
      DO 104 K=2,NS2
         KC = NP1-K
         T1 = X(1,K)+X(1,KC)
         T2 = X(1,K)-X(1,KC)
         DSUM = DSUM+fft_sign%WSAVE(KC)*T2
         T2 = fft_sign%WSAVE(K)*T2
         X(1,K) = T1-T2
         X(1,KC) = T1+T2
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
      X(1,NS2+1) = X(1,NS2+1)+X(1,NS2+1)
  124 LENX = INC*(NM1-1)  + 1
      LNSV = NM1 + INT(LOG(REAL(NM1))/LOG(2._rk)) + 4
      LNWK = NM1
!
      !!CALL RFFT1F(NM1,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK,fft_sign,IER1)
      !!CALL RFFT1F(NM1,INC,X,LENX,WORK,LNWK,fft_sign,IER1)
      CALL RFFT1F(NM1,INC,X,LENX,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSTB1',-5)
        RETURN
      ENDIF
!
      FNM1S2 = REAL(NM1,kind=rk)/2._rk
      DSUM = .5_dk*DSUM
      X(1,1) = FNM1S2*X(1,1)
      IF(MOD(NM1,2) .NE. 0) GO TO 30
      X(1,NM1) = X(1,NM1)+X(1,NM1)
   30 FNM1S4 = REAL(NM1,kind=rk)/4._rk
      DO 105 I=3,N,2
         XI = FNM1S4*X(1,I)
         X(1,I) = FNM1S4*X(1,I-1)
         X(1,I-1) = DSUM
         DSUM = DSUM+XI
  105 CONTINUE
      IF (MODN .NE. 0) RETURN
         X(1,N) = DSUM
  106 continue
      END subroutine COSTB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSTF1(N,INC,X,WSAVE,WORK,fft_sign, IER)
      !!SUBROUTINE COSTF1(N,INC,X,WORK,fft_sign, IER)
      SUBROUTINE COSTF1(N,INC,X,fft_sign, IER)
      INTEGER           :: N,INC,IER
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      REAL(kind=dk)     :: SNM1
     !! REAL(kind=rk)     :: WSAVE(*)
      REAL(kind=rk)     :: X(INC,*)
      real(kind=rk)     :: T1,T2
      !!real(kind=rk)     :: WORK(*)    ! Kuester
      real(kind=rk)     :: TX2,X1P3
      real(kind=rk)     :: XI
      real(kind=rk)     :: X1H
      real(kind=dk)    :: DSUM
      type(fft_real_sign_type) :: fft_sign


      IER = 0
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      IF (N-2) 200,101,102
  101 X1H = X(1,1)+X(1,2)
      X(1,2) = .5_rk*(X(1,1)-X(1,2))
      X(1,1) = .5_rk*X1H
      GO TO 200
  102 IF (N .GT. 3) GO TO 103
      X1P3 = X(1,1)+X(1,3)
      TX2 = X(1,2)+X(1,2)
      X(1,2) = .5_rk*(X(1,1)-X(1,3))
      X(1,1) = .25_rk*(X1P3+TX2)
      X(1,3) = .25_rk*(X1P3-TX2)
      GO TO 200
  103 DSUM = X(1,1)-X(1,N)
      X(1,1) = X(1,1)+X(1,N)
      DO 104 K=2,NS2
         KC = NP1-K
         T1 = X(1,K)+X(1,KC)
         T2 = X(1,K)-X(1,KC)
         DSUM = DSUM+fft_sign%WSAVE(KC)*T2
         T2 = fft_sign%WSAVE(K)*T2
         X(1,K) = T1-T2
         X(1,KC) = T1+T2
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
      X(1,NS2+1) = X(1,NS2+1)+X(1,NS2+1)
  124 LENX = INC*(NM1-1)  + 1
      LNSV = NM1 + INT(LOG(REAL(NM1))/LOG(2._rk)) + 4
      LNWK = NM1
!
      !!CALL RFFT1F(NM1,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK,fft_sign,IER1)
      !!CALL RFFT1F(NM1,INC,X,LENX,WORK,LNWK,fft_sign,IER1)
      CALL RFFT1F(NM1,INC,X,LENX,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSTF1',-5)
        GO TO 200
      ENDIF
!
      SNM1 = 1._rk/REAL(NM1,kind=dk)
      DSUM = SNM1*DSUM
      IF(MOD(NM1,2) .NE. 0) GO TO 30
      X(1,NM1) = X(1,NM1)+X(1,NM1)
   30 DO 105 I=3,N,2
         XI = .5_rk*X(1,I)
         X(1,I) = .5_rk*X(1,I-1)
         X(1,I-1) = DSUM
         DSUM = DSUM+XI
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 117
      X(1,N) = DSUM
  117 X(1,1) = .5_rk*X(1,1)
      X(1,N) = .5_rk*X(1,N)
  200 continue
      END subroutine COSTF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      ! Entrance
      !!SUBROUTINE COSTMB (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, DSUM, LENWRK,fft_sign, IER)
      SUBROUTINE COSTMB (LOT, JUMP, N, INC, X, LENX, WORK, DSUM, LENWRK,fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IW1,IW2,IER1
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      real(kind=dk)    :: DSUM(*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSTMB', 6)
        GO TO 100
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('COSTMB', 8)
!!        GO TO 100
      ELSEIF (LENWRK .LT. LOT*(N+1)) THEN
        IER = 3
        CALL XERFFT ('COSTMB', 10)
        GO TO 100
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('COSTMB', -1)
        GO TO 100
      ENDIF
!
      IW1 = LOT+LOT+1
      !Kuester CALL MCSTB1(LOT,JUMP,N,INC,X,WSAVE,WORK,WORK(IW1),IER1)
      !!CALL MCSTB1(LOT,JUMP,N,INC,X,WSAVE,DSUM,WORK(IW1),fft_sign,IER1)
      CALL MCSTB1(LOT,JUMP,N,INC,X,DSUM,WORK(IW1),fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSTMB',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine COSTMB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      ! Entrance
      !!SUBROUTINE COSTMF (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, DSUM, LENWRK, fft_sign, IER)
      SUBROUTINE COSTMF (LOT, JUMP, N, INC, X, LENX, WORK, DSUM, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IW1,IW2,IER1
 !!     real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      real(kind=dk)    :: DSUM(*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSTMF', 6)
        GO TO 100
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('COSTMF', 8)
!!        GO TO 100
      ELSEIF (LENWRK .LT. LOT*(N+1)) THEN
        IER = 3
        CALL XERFFT ('COSTMF', 10)
        GO TO 100
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('COSTMF', -1)
        GO TO 100
      ENDIF
!
      IW1 = LOT+LOT+1
      !Kuester CALL MCSTF1(LOT,JUMP,N,INC,X,WSAVE,WORK,WORK(IW1),IER1)
      !!CALL MCSTF1(LOT,JUMP,N,INC,X,WSAVE,DSUM,WORK(IW1), fft_sign,IER1)
      CALL MCSTF1(LOT,JUMP,N,INC,X,DSUM,WORK(IW1), fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSTMF',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine COSTMF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE COSTMI (N, WSAVE, LENSAV, IER)
      SUBROUTINE COSTMI (N, IER)
      integer     :: N
      integer     ::  IER
      integer     :: LENSAV
      integer           :: NM1,NP1,NS2,K,KC,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: PI,DT,FK
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
 !!     IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
 !!       IER = 2
 !!       CALL XERFFT ('COSTMI', 3)
 !!       GO TO 300
 !!     ENDIF
!
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2

      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(NM1))

      PI = 4._rk*ATAN(1._rk)
      PI = PI_long
      DT = PI/REAL(NM1,kind=rk)
      FK = 0._rk
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1._rk
         !!WSAVE(K) = 2._rk*SIN(FK*DT)
         !!WSAVE(KC) = 2._rk*COS(FK*DT)
         fft_sign%WSAVE(K) = 2._rk*SIN(FK*DT)
         fft_sign%WSAVE(KC) = 2._rk*COS(FK*DT)
  101 CONTINUE

      !!LNSV = NM1 + INT(LOG(REAL(NM1))/LOG(2._rk)) +4
      !!CALL RFFTMI (NM1, WSAVE(N+1), LNSV,fft_sign, IER1)
      CALL RFFTMI (NM1, fft_sign, IER1)

      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSTMI',-5)
      ENDIF
  300 CONTINUE
      END subroutine COSTMI
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE FACTOR_integer (N,NF,FAC)
      integer            :: N,NF
      integer            :: NL,NTRY,NQ,NR
      integer            :: J
      integer              :: FAC(*)
      INTEGER NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
!
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      !!write(*,'(*(a,i0))') 'J ',J,' NQ ',NQ,' NR ',NR,' NF ',NF
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF) = NTRY
      NL = NQ
      IF (NL .NE. 1) GO TO 104

      END subroutine FACTOR_integer
!************************************************************************************************************************
      SUBROUTINE FACTOR_real (N,NF,FAC)
      integer            :: N,NF
      integer            :: NL,NTRY,NQ,NR
      integer            :: J
      REAL(kind=rk)     :: FAC(*)
      INTEGER NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
!
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF) = NTRY
      NL = NQ
      IF (NL .NE. 1) GO TO 104

      END subroutine FACTOR_real
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MCFTI1_orig (N,WA,FNF,FAC)
      integer           :: N
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: fnf  
      REAL(kind=rk)     :: WA(*),FAC(*)
      integer           :: ll
!


      CALL FACTOR (N,NF,FAC)
         write(*,*)
         write(*,'(a,i0)') 'MCFTI1_orig: prime factors of ',n
      do ll=1,NF 
         write(*,'(2i5)') ll,int(fac(ll))
      enddo
         write(*,*)



      FNF = NF
      IW = 1
      L1 = 1
      DO 110 K1=1,NF
         IP = FAC(K1)
         L2 = L1*IP
         IDO = N/L2
         CALL TABLES (IDO,IP,WA(IW))
         IW = IW+(IP-1)*(IDO+IDO)
         L1 = L2
  110 CONTINUE
      END subroutine MCFTI1_orig
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MCFTI1 (N,WA,FNF,FAC,fft_sign)
      SUBROUTINE MCFTI1 (N,fft_sign)
      integer           :: N
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: fnf  
      !!REAL(kind=rk)     :: WA(*),FAC(*)
      integer           :: ll
      integer             :: ibase2
      type(fft_real_sign_type) :: fft_sign
!
      if( fft_sign%switch == 1 .and. fft_sign%N == N .and. allocated(fft_sign%FAC) .and. allocated(fft_sign%WA)) then
           goto 1111
      endif

      fft_sign%switch = 1

      fft_sign%N = N

      ibase2 = NINT(LOG(real(N,kind=dk))/LOG(2.0_dk))+1
      if(allocated(fft_sign%FAC)) deallocate(fft_sign%FAC)
      allocate(fft_sign%FAC(ibase2))
      write(*,'(a,i0)') 'MCFTI1: ibase2 ',ibase2

      CALL FACTOR (N,fft_sign%NF,fft_sign%FAC)

     !! CALL FACTOR (N,NF,FAC)
         write(*,*)
         write(*,'(a,i0)') 'MCFTI1: prime factors of ',n
      do ll=1,fft_sign%nf
         write(*,'(i5)',advance='no') int(fft_sign%fac(ll))
      enddo
         write(*,*)

      IW = 1
      L1 = 1
      do K1=1,fft_sign%NF
         IP = fft_sign%FAC(K1)
         L2 = L1*IP
         IDO = N/L2
         IW = IW+(IP-1)*(IDO+IDO)
         L1 = L2
      enddo

      !!print*,'MCFTI1: iw=',iw
      if(allocated(fft_sign%WA)) deallocate(fft_sign%WA)
      allocate(fft_sign%WA(IW))

      IW = 1
      L1 = 1
      do K1=1,fft_sign%NF
         IP = fft_sign%FAC(K1)
         L2 = L1*IP
         IDO = N/L2
         !!write(*,*) 'MCFTI1: ',k1,ido,ip,l1,l2
         CALL TABLES(IDO,IP,fft_sign%WA(IW))
         IW = IW+(IP-1)*(IDO+IDO)
         L1 = L2
      enddo


!!      FNF = NF
!!      IW = 1
!!      L1 = 1
!!      DO 110 K1=1,NF
!!         IP = FAC(K1)
!!         L2 = L1*IP
!!         IDO = N/L2
!!         CALL TABLES (IDO,IP,WA(IW))
!!         IW = IW+(IP-1)*(IDO+IDO)
!!         L1 = L2
!!  110 CONTINUE
      1111 continue

      END subroutine MCFTI1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MCSQB1 (LOT,JUMP,N,INC,X,WSAVE,WORK, fft_sign,IER)
      SUBROUTINE MCSQB1 (LOT,JUMP,N,INC,X,WORK, fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: IER1
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: XIM1 
      !!real(kind=rk)     :: WSAVE(*)
      real(kind=rk)     :: X(INC,*),WORK(LOT,*)
      type(fft_real_sign_type) :: fft_sign
      IER = 0
      LJ = (LOT-1)*JUMP+1
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 I=3,N,2
         DO 201 M=1,LJ,JUMP
         XIM1 = X(M,I-1)+X(M,I)
         X(M,I) = .5_rk*(X(M,I-1)-X(M,I))
         X(M,I-1) = .5_rk*XIM1
 201     CONTINUE
  101 CONTINUE
      DO 301 M=1,LJ,JUMP
      X(M,1) = .5_rk*X(M,1)
 301  CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) GO TO 302
      DO 303 M=1,LJ,JUMP
      X(M,N) = .5_rk*X(M,N)
 303  CONTINUE
 302  CONTINUE
      LENX = (LOT-1)*JUMP + INC*(N-1)  + 1
      LNSV = N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
      LNWK = LOT*N
!
      !!CALL RFFTMB(LOT,JUMP,N,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK, fft_sign,IER1)
      CALL RFFTMB(LOT,JUMP,N,INC,X,LENX,WORK,LNWK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('MCSQB1',-5)
        GO TO 400
      ENDIF
!
      DO 102 K=2,NS2
         KC = NP2-K
         M1 = 0
         DO 202 M=1,LJ,JUMP
         M1 = M1 + 1
         WORK(M1,K) = fft_sign%WSAVE(K-1)*X(M,KC)+fft_sign%WSAVE(KC-1)*X(M,K)
         WORK(M1,KC) = fft_sign%WSAVE(K-1)*X(M,K)-fft_sign%WSAVE(KC-1)*X(M,KC)
 202     CONTINUE
  102 CONTINUE
      IF (MODN .NE. 0) GO TO 305
      DO 304 M=1,LJ,JUMP
         X(M,NS2+1) = fft_sign%WSAVE(NS2)*(X(M,NS2+1)+X(M,NS2+1))
 304     CONTINUE
 305  DO 103 K=2,NS2
         KC = NP2-K
         M1 = 0
         DO 203 M=1,LJ,JUMP
            M1 = M1 + 1
            X(M,K) = WORK(M1,K)+WORK(M1,KC)
            X(M,KC) = WORK(M1,K)-WORK(M1,KC)
 203     CONTINUE
  103 CONTINUE
      DO 104 M=1,LJ,JUMP
      X(M,1) = X(M,1)+X(M,1)
 104  CONTINUE
  400 CONTINUE
      END subroutine MCSQB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MCSQF1 (LOT,JUMP,N,INC,X,WSAVE,WORK, fft_sign,IER)
      SUBROUTINE MCSQF1 (LOT,JUMP,N,INC,X,WORK, fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IER1
      real(kind=rk)     :: XIM1 
     !! real(kind=rk)     :: WSAVE(*)
      real(kind=rk)     :: X(INC,*),WORK(LOT,*)
      type(fft_real_sign_type) :: fft_sign

      IER = 0
      LJ = (LOT-1)*JUMP+1
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 K=2,NS2
         KC = NP2-K
         M1 = 0
         DO 201 M=1,LJ,JUMP
         M1 = M1 + 1
         WORK(M1,K)  = X(M,K)+X(M,KC)
         WORK(M1,KC) = X(M,K)-X(M,KC)
 201     CONTINUE
  101 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) GO TO 301
         M1 = 0
         DO 202 M=1,LJ,JUMP
         M1 = M1 + 1
         WORK(M1,NS2+1) = X(M,NS2+1)+X(M,NS2+1)
 202     CONTINUE
 301     DO 102 K=2,NS2
         KC = NP2-K
         M1 = 0
         DO 302 M=1,LJ,JUMP
         M1 = M1 + 1
         X(M,K)  = fft_sign%WSAVE(K-1)*WORK(M1,KC)+fft_sign%WSAVE(KC-1)*WORK(M1,K)
         X(M,KC) = fft_sign%WSAVE(K-1)*WORK(M1,K) -fft_sign%WSAVE(KC-1)*WORK(M1,KC)
 302     CONTINUE
  102 CONTINUE
      IF (MODN .NE. 0) GO TO 303
      M1 = 0
      DO 304 M=1,LJ,JUMP
         M1 = M1 + 1
         X(M,NS2+1) = fft_sign%WSAVE(NS2)*WORK(M1,NS2+1)
 304  CONTINUE
 303  CONTINUE
      LENX = (LOT-1)*JUMP + INC*(N-1)  + 1
      LNSV = N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) + 4
      LNWK = LOT*N
!
      !!CALL RFFTMF(LOT,JUMP,N,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK, fft_sign,IER1)
      CALL RFFTMF(LOT,JUMP,N,INC,X,LENX,WORK,LNWK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('MCSQF1',-5)
        GO TO 400
      ENDIF
!
      DO 103 I=3,N,2
         DO 203 M=1,LJ,JUMP
            XIM1 = .5_rk*(X(M,I-1)+X(M,I))
            X(M,I) = .5_rk*(X(M,I-1)-X(M,I))
            X(M,I-1) = XIM1
 203     CONTINUE
  103 CONTINUE
  400 CONTINUE
      END subroutine MCSQF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MCSTB1(LOT,JUMP,N,INC,X,WSAVE,DSUM,WORK, fft_sign,IER)
      SUBROUTINE MCSTB1(LOT,JUMP,N,INC,X,DSUM,WORK, fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: IER1
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      !!REAL(kind=rk)     :: WSAVE(*)
      REAL(kind=rk)     :: X(INC,*)
      real(kind=rk)     :: WORK(*)   ! Kuester
      real(kind=rk)     :: T1,T2
      real(kind=rk)     :: XI
      real(kind=rk)     :: X1H
      real(kind=rk)     :: X1P3,X2
      real(kind=rk)     :: FNM1S2,FNM1S4
      real(kind=dk)    :: DSUM(*) ! no value on entrance
      type(fft_real_sign_type) :: fft_sign

      IER = 0
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 106,101,102
  101 DO 111 M=1,LJ,JUMP
         X1H = X(M,1)+X(M,2)
         X(M,2) = X(M,1)-X(M,2)
         X(M,1) = X1H
  111 CONTINUE
      RETURN
  102 IF (N .GT. 3) GO TO 103
      DO 112 M=1,LJ,JUMP
         X1P3 = X(M,1)+X(M,3)
         X2 = X(M,2)
         X(M,2) = X(M,1)-X(M,3)
         X(M,1) = X1P3+X2
         X(M,3) = X1P3-X2
  112 CONTINUE
      RETURN
 103  DO 118 M=1,LJ,JUMP
      X(M,1) = X(M,1)+X(M,1)
      X(M,N) = X(M,N)+X(M,N)
 118  CONTINUE
      M1 = 0
      DO 113 M=1,LJ,JUMP
         M1 = M1+1
         DSUM(M1) = X(M,1)-X(M,N)
         X(M,1) = X(M,1)+X(M,N)
  113 CONTINUE
      DO 104 K=2,NS2
         M1 = 0
         DO 114 M=1,LJ,JUMP
         M1 = M1+1
         KC = NP1-K
         T1 = X(M,K)+X(M,KC)
         T2 = X(M,K)-X(M,KC)
         DSUM(M1) = DSUM(M1)+fft_sign%WSAVE(KC)*T2
         T2 = fft_sign%WSAVE(K)*T2
         X(M,K) = T1-T2
         X(M,KC) = T1+T2
  114    CONTINUE
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
         DO 123 M=1,LJ,JUMP
         X(M,NS2+1) = X(M,NS2+1)+X(M,NS2+1)
  123    CONTINUE
 124  CONTINUE
      LENX = (LOT-1)*JUMP + INC*(NM1-1)  + 1
      LNSV = NM1 + INT(LOG(REAL(NM1))/LOG(2._rk)) + 4
      LNWK = LOT*NM1
!
      !!CALL RFFTMF(LOT,JUMP,NM1,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK, fft_sign,IER1)
      CALL RFFTMF(LOT,JUMP,NM1,INC,X,LENX,WORK,LNWK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('MCSTB1',-5)
        GO TO 106
      ENDIF
!
      FNM1S2 = REAL(NM1,kind=rk)/2._rk
      M1 = 0
      DO 10 M=1,LJ,JUMP
      M1 = M1+1
      DSUM(M1) = .5_dk*DSUM(M1)
      X(M,1) = FNM1S2*X(M,1)
   10 CONTINUE
      IF(MOD(NM1,2) .NE. 0) GO TO 30
      DO 20 M=1,LJ,JUMP
      X(M,NM1) = X(M,NM1)+X(M,NM1)
   20 CONTINUE
 30   FNM1S4 = REAL(NM1,kind=rk)/4._rk
      DO 105 I=3,N,2
         M1 = 0
         DO 115 M=1,LJ,JUMP
            M1 = M1+1
            XI = FNM1S4*X(M,I)
            X(M,I) = FNM1S4*X(M,I-1)
            X(M,I-1) = DSUM(M1)
            DSUM(M1) = DSUM(M1)+XI
  115 CONTINUE
  105 CONTINUE
      IF (MODN .NE. 0) RETURN
      M1 = 0
      DO 116 M=1,LJ,JUMP
         M1 = M1+1
         X(M,N) = DSUM(M1)
  116 CONTINUE
  106 CONTINUE
      END subroutine MCSTB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MCSTF1(LOT,JUMP,N,INC,X,WSAVE,DSUM,WORK, fft_sign,IER)
      SUBROUTINE MCSTF1(LOT,JUMP,N,INC,X,DSUM,WORK, fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: IER1
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      !!REAL(kind=rk)     :: WSAVE(*)
      REAL(kind=rk)     :: X(INC,*)
      real(kind=rk)     :: WORK(*) ! Kuester
      real(kind=rk)     :: T1,T2
      REAL(kind=dk)     :: SNM1
      real(kind=rk)     :: TX2,X1P3
      real(kind=rk)     :: XI
      real(kind=rk)     :: X1H
      real(kind=dk)    :: DSUM(*) ! no value on entrance
      type(fft_real_sign_type) :: fft_sign
      IER = 0
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 200,101,102
  101 DO 111 M=1,LJ,JUMP
         X1H = X(M,1)+X(M,2)
         X(M,2) = .5_rk*(X(M,1)-X(M,2))
         X(M,1) = .5_rk*X1H
  111 CONTINUE
      GO TO 200
  102 IF (N .GT. 3) GO TO 103
      DO 112 M=1,LJ,JUMP
         X1P3 = X(M,1)+X(M,3)
         TX2 = X(M,2)+X(M,2)
         X(M,2) = .5_rk*(X(M,1)-X(M,3))
         X(M,1) = .25_rk*(X1P3+TX2)
         X(M,3) = .25_rk*(X1P3-TX2)
  112 CONTINUE
      GO TO 200
  103 M1 = 0
      DO 113 M=1,LJ,JUMP
         M1 = M1+1
         DSUM(M1) = X(M,1)-X(M,N)
         X(M,1) = X(M,1)+X(M,N)
  113 CONTINUE
      DO 104 K=2,NS2
         M1 = 0
         DO 114 M=1,LJ,JUMP
         M1 = M1+1
         KC = NP1-K
         T1 = X(M,K)+X(M,KC)
         T2 = X(M,K)-X(M,KC)
         DSUM(M1) = DSUM(M1)+fft_sign%WSAVE(KC)*T2
         T2 = fft_sign%WSAVE(K)*T2
         X(M,K) = T1-T2
         X(M,KC) = T1+T2
  114    CONTINUE
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
         DO 123 M=1,LJ,JUMP
         X(M,NS2+1) = X(M,NS2+1)+X(M,NS2+1)
  123    CONTINUE
 124  CONTINUE
      LENX = (LOT-1)*JUMP + INC*(NM1-1)  + 1
      LNSV = NM1 + INT(LOG(REAL(NM1))/LOG(2._rk)) + 4
      LNWK = LOT*NM1
!
      !!CALL RFFTMF(LOT,JUMP,NM1,INC,X,LENX,WSAVE(N+1),LNSV,WORK,LNWK, fft_sign,IER1)
      CALL RFFTMF(LOT,JUMP,NM1,INC,X,LENX,WORK,LNWK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('MCSTF1',-5)
        GO TO 200
      ENDIF
!
      SNM1 = 1._rk/REAL(NM1,kind=dk)
      DO 10 M=1,LOT
      DSUM(M) = SNM1*DSUM(M)
   10 CONTINUE
      IF(MOD(NM1,2) .NE. 0) GO TO 30
      DO 20 M=1,LJ,JUMP
      X(M,NM1) = X(M,NM1)+X(M,NM1)
   20 CONTINUE
 30   DO 105 I=3,N,2
         M1 = 0
         DO 115 M=1,LJ,JUMP
            M1 = M1+1
            XI = .5_rk*X(M,I)
            X(M,I) = .5_rk*X(M,I-1)
            X(M,I-1) = DSUM(M1)
            DSUM(M1) = DSUM(M1)+XI
  115 CONTINUE
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 117
      M1 = 0
      DO 116 M=1,LJ,JUMP
         M1 = M1+1
         X(M,N) = DSUM(M1)
  116 CONTINUE
 117  DO 118 M=1,LJ,JUMP
      X(M,1) = .5_rk*X(M,1)
      X(M,N) = .5_rk*X(M,N)
 118  CONTINUE
!
 200  CONTINUE
      END subroutine MCSTF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADB2 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      REAL(kind=rk)     :: CC(IN1,IDO,2,L1), CH(IN2,IDO,L1,2), WA1(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      DO 101 K=1,L1
          M2 = M2S
          DO 1001 M1=1,M1D,IM1
          M2 = M2+IM2
         CH(M2,1,K,1) = CC(M1,1,1,K)+CC(M1,IDO,2,K)
         CH(M2,1,K,2) = CC(M1,1,1,K)-CC(M1,IDO,2,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               M2 = M2S
               DO 1002 M1=1,M1D,IM1
               M2 = M2+IM2
        CH(M2,I-1,K,1) = CC(M1,I-1,1,K)+CC(M1,IC-1,2,K)
        CH(M2,I,K,1) = CC(M1,I,1,K)-CC(M1,IC,2,K)
        CH(M2,I-1,K,2) = WA1(I-2)*(CC(M1,I-1,1,K)-CC(M1,IC-1,2,K))-WA1(I-1)*(CC(M1,I,1,K)+CC(M1,IC,2,K))
        CH(M2,I,K,2) = WA1(I-2)*(CC(M1,I,1,K)+CC(M1,IC,2,K))+WA1(I-1)*(CC(M1,I-1,1,K)-CC(M1,IC-1,2,K))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
          M2 = M2S
          DO 1003 M1=1,M1D,IM1
          M2 = M2+IM2
         CH(M2,IDO,K,1) = CC(M1,IDO,1,K)+CC(M1,IDO,1,K)
         CH(M2,IDO,K,2) = -(CC(M1,1,2,K)+CC(M1,1,2,K))
 1003     CONTINUE
  106 CONTINUE
  107 continue
      END subroutine MRADB2
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADB3 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1,WA2)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TAUR,TAUI
      REAL(kind=rk)     :: CC(IN1,IDO,3,L1),CH(IN2,IDO,L1,3),WA1(IDO),WA2(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      ARG = 2._dk*4._dk*ATAN(1.0_dk)/3._dk
      ARG = 2._dk*PI_long/3._dk
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
          M2 = M2S
          DO 1001 M1=1,M1D,IM1
          M2 = M2+IM2
         CH(M2,1,K,1) = CC(M1,1,1,K)+2._rk*CC(M1,IDO,2,K)
         CH(M2,1,K,2) = CC(M1,1,1,K)+(2._rk*TAUR)*CC(M1,IDO,2,K)-(2._rk*TAUI)*CC(M1,1,3,K)
         CH(M2,1,K,3) = CC(M1,1,1,K)+(2._rk*TAUR)*CC(M1,IDO,2,K)+2._rk*TAUI*CC(M1,1,3,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
               M2 = M2S
               DO 1002 M1=1,M1D,IM1
               M2 = M2+IM2
        CH(M2,I-1,K,1) = CC(M1,I-1,1,K)+(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))
        CH(M2,I,K,1) = CC(M1,I,1,K)+(CC(M1,I,3,K)-CC(M1,IC,2,K))
        CH(M2,I-1,K,2) = WA1(I-2)*&
     & ((CC(M1,I-1,1,K)+TAUR*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K)))-&
     & (TAUI*(CC(M1,I,3,K)+CC(M1,IC,2,K))))&
     &                   -WA1(I-1)*&
     & ((CC(M1,I,1,K)+TAUR*(CC(M1,I,3,K)-CC(M1,IC,2,K)))+&
     & (TAUI*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))))
            CH(M2,I,K,2) = WA1(I-2)*&
     & ((CC(M1,I,1,K)+TAUR*(CC(M1,I,3,K)-CC(M1,IC,2,K)))+&
     & (TAUI*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))))&
     &                  +WA1(I-1)*&
     & ((CC(M1,I-1,1,K)+TAUR*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K)))-&
     & (TAUI*(CC(M1,I,3,K)+CC(M1,IC,2,K))))
              CH(M2,I-1,K,3) = WA2(I-2)*&
     & ((CC(M1,I-1,1,K)+TAUR*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K)))+&
     & (TAUI*(CC(M1,I,3,K)+CC(M1,IC,2,K))))&
     &   -WA2(I-1)*&
     & ((CC(M1,I,1,K)+TAUR*(CC(M1,I,3,K)-CC(M1,IC,2,K)))-&
     & (TAUI*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))))
            CH(M2,I,K,3) = WA2(I-2)*&
     & ((CC(M1,I,1,K)+TAUR*(CC(M1,I,3,K)-CC(M1,IC,2,K)))-&
     & (TAUI*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))))&
     &                 +WA2(I-1)*&
     & ((CC(M1,I-1,1,K)+TAUR*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K)))+&
     & (TAUI*(CC(M1,I,3,K)+CC(M1,IC,2,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      END subroutine MRADB3
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADB4 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1,WA2,WA3)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SQRT2
      REAL(kind=rk)     :: CC(IN1,IDO,4,L1),CH(IN2,IDO,L1,4),WA1(IDO),WA2(IDO),WA3(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      SQRT2=SQRT(2._rk)
      DO 101 K=1,L1
          M2 = M2S
          DO 1001 M1=1,M1D,IM1
          M2 = M2+IM2
         CH(M2,1,K,3) = (CC(M1,1,1,K)+CC(M1,IDO,4,K))-(CC(M1,IDO,2,K)+CC(M1,IDO,2,K))
         CH(M2,1,K,1) = (CC(M1,1,1,K)+CC(M1,IDO,4,K))+(CC(M1,IDO,2,K)+CC(M1,IDO,2,K))
         CH(M2,1,K,4) = (CC(M1,1,1,K)-CC(M1,IDO,4,K))+(CC(M1,1,3,K)+CC(M1,1,3,K))
         CH(M2,1,K,2) = (CC(M1,1,1,K)-CC(M1,IDO,4,K))-(CC(M1,1,3,K)+CC(M1,1,3,K))
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               M2 = M2S
               DO 1002 M1=1,M1D,IM1
               M2 = M2+IM2
        CH(M2,I-1,K,1) = (CC(M1,I-1,1,K)+CC(M1,IC-1,4,K))+(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))
        CH(M2,I,K,1) = (CC(M1,I,1,K)-CC(M1,IC,4,K))+(CC(M1,I,3,K)-CC(M1,IC,2,K))
        CH(M2,I-1,K,2)=WA1(I-2)*((CC(M1,I-1,1,K)-CC(M1,IC-1,4,K))&
     &  -(CC(M1,I,3,K)+CC(M1,IC,2,K)))-WA1(I-1)&
     &  *((CC(M1,I,1,K)+CC(M1,IC,4,K))+(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K)))
        CH(M2,I,K,2)=WA1(I-2)*((CC(M1,I,1,K)+CC(M1,IC,4,K))&
     &  +(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K)))+WA1(I-1)&
     &  *((CC(M1,I-1,1,K)-CC(M1,IC-1,4,K))-(CC(M1,I,3,K)+CC(M1,IC,2,K)))
        CH(M2,I-1,K,3)=WA2(I-2)*((CC(M1,I-1,1,K)+CC(M1,IC-1,4,K))&
     &  -(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K)))-WA2(I-1)&
     &  *((CC(M1,I,1,K)-CC(M1,IC,4,K))-(CC(M1,I,3,K)-CC(M1,IC,2,K)))
        CH(M2,I,K,3)=WA2(I-2)*((CC(M1,I,1,K)-CC(M1,IC,4,K))&
     &  -(CC(M1,I,3,K)-CC(M1,IC,2,K)))+WA2(I-1)&
     &  *((CC(M1,I-1,1,K)+CC(M1,IC-1,4,K))-(CC(M1,I-1,3,K)&
     &  +CC(M1,IC-1,2,K)))
        CH(M2,I-1,K,4)=WA3(I-2)*((CC(M1,I-1,1,K)-CC(M1,IC-1,4,K))&
     &  +(CC(M1,I,3,K)+CC(M1,IC,2,K)))-WA3(I-1)&
     & *((CC(M1,I,1,K)+CC(M1,IC,4,K))-(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K)))
        CH(M2,I,K,4)=WA3(I-2)*((CC(M1,I,1,K)+CC(M1,IC,4,K))&
     &  -(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K)))+WA3(I-1)&
     &  *((CC(M1,I-1,1,K)-CC(M1,IC-1,4,K))+(CC(M1,I,3,K)+CC(M1,IC,2,K)))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
               M2 = M2S
               DO 1003 M1=1,M1D,IM1
               M2 = M2+IM2
         CH(M2,IDO,K,1) = (CC(M1,IDO,1,K)+CC(M1,IDO,3,K))+(CC(M1,IDO,1,K)+CC(M1,IDO,3,K))
         CH(M2,IDO,K,2) = SQRT2*((CC(M1,IDO,1,K)-CC(M1,IDO,3,K))-(CC(M1,1,2,K)+CC(M1,1,4,K)))
         CH(M2,IDO,K,3) = (CC(M1,1,4,K)-CC(M1,1,2,K))+(CC(M1,1,4,K)-CC(M1,1,2,K))
         CH(M2,IDO,K,4) = -SQRT2*((CC(M1,IDO,1,K)-CC(M1,IDO,3,K))+(CC(M1,1,2,K)+CC(M1,1,4,K)))
 1003          CONTINUE
  106 CONTINUE
  107 continue
      END subroutine MRADB4
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADB5 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1,WA2,WA3,WA4)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      REAL(kind=rk)     :: CC(IN1,IDO,5,L1),CH(IN2,IDO,L1,5),WA1(IDO),WA2(IDO),WA3(IDO),WA4(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/5._dk
      ARG = 2._dk*PI_long/5._dk
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2._rk*ARG)
      TI12=SIN(2._rk*ARG)
      DO 101 K=1,L1
      M2 = M2S
      DO 1001 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(M2,1,K,1) = CC(M1,1,1,K)+2._rk*CC(M1,IDO,2,K)+2._rk*CC(M1,IDO,4,K)
         CH(M2,1,K,2) = (CC(M1,1,1,K)+TR11*2._rk*CC(M1,IDO,2,K)&
     &   +TR12*2._rk*CC(M1,IDO,4,K))-(TI11*2._rk*CC(M1,1,3,K)&
     &   +TI12*2._rk*CC(M1,1,5,K))
         CH(M2,1,K,3) = (CC(M1,1,1,K)+TR12*2._rk*CC(M1,IDO,2,K)&
     &   +TR11*2._rk*CC(M1,IDO,4,K))-(TI12*2._rk*CC(M1,1,3,K)&
     &   -TI11*2._rk*CC(M1,1,5,K))
         CH(M2,1,K,4) = (CC(M1,1,1,K)+TR12*2._rk*CC(M1,IDO,2,K)&
     &   +TR11*2._rk*CC(M1,IDO,4,K))+(TI12*2._rk*CC(M1,1,3,K)&
     &   -TI11*2._rk*CC(M1,1,5,K))
         CH(M2,1,K,5) = (CC(M1,1,1,K)+TR11*2._rk*CC(M1,IDO,2,K)&
     &   +TR12*2._rk*CC(M1,IDO,4,K))+(TI11*2._rk*CC(M1,1,3,K)&
     &   +TI12*2._rk*CC(M1,1,5,K))
 1001          CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            M2 = M2S
      DO 1002 M1=1,M1D,IM1
        M2 = M2+IM2
        CH(M2,I-1,K,1) = CC(M1,I-1,1,K)+(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))+(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K))
        CH(M2,I,K,1) = CC(M1,I,1,K)+(CC(M1,I,3,K)-CC(M1,IC,2,K))+(CC(M1,I,5,K)-CC(M1,IC,4,K))
        CH(M2,I-1,K,2) = WA1(I-2)*((CC(M1,I-1,1,K)+TR11*&
     &  (CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))+TR12&
     &  *(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))-(TI11*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))+TI12*(CC(M1,I,5,K)+CC(M1,IC,4,K))))&
     &  -WA1(I-1)*((CC(M1,I,1,K)+TR11*(CC(M1,I,3,K)-CC(M1,IC,2,K))&
     &  +TR12*(CC(M1,I,5,K)-CC(M1,IC,4,K)))+(TI11*(CC(M1,I-1,3,K)&
     &  -CC(M1,IC-1,2,K))+TI12*(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))
        CH(M2,I,K,2) = WA1(I-2)*((CC(M1,I,1,K)+TR11*(CC(M1,I,3,K)&
     &  -CC(M1,IC,2,K))+TR12*(CC(M1,I,5,K)-CC(M1,IC,4,K)))&
     &  +(TI11*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))+TI12&
     &  *(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))+WA1(I-1)&
     &  *((CC(M1,I-1,1,K)+TR11*(CC(M1,I-1,3,K)&
     &  +CC(M1,IC-1,2,K))+TR12*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))&
     &  -(TI11*(CC(M1,I,3,K)+CC(M1,IC,2,K))+TI12&
     &  *(CC(M1,I,5,K)+CC(M1,IC,4,K))))
        CH(M2,I-1,K,3) = WA2(I-2)&
     &  *((CC(M1,I-1,1,K)+TR12*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))&
     &  +TR11*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))-(TI12*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))-TI11*(CC(M1,I,5,K)+CC(M1,IC,4,K))))&
     & -WA2(I-1)&
     & *((CC(M1,I,1,K)+TR12*(CC(M1,I,3,K)-&
     &  CC(M1,IC,2,K))+TR11*(CC(M1,I,5,K)-CC(M1,IC,4,K)))&
     &  +(TI12*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))-TI11&
     &  *(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))
        CH(M2,I,K,3) = WA2(I-2)&
     & *((CC(M1,I,1,K)+TR12*(CC(M1,I,3,K)-&
     &  CC(M1,IC,2,K))+TR11*(CC(M1,I,5,K)-CC(M1,IC,4,K)))&
     &  +(TI12*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))-TI11&
     &  *(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))&
     &  +WA2(I-1)&
     &  *((CC(M1,I-1,1,K)+TR12*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))&
     &  +TR11*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))-(TI12*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))-TI11*(CC(M1,I,5,K)+CC(M1,IC,4,K))))
        CH(M2,I-1,K,4) = WA3(I-2)&
     &  *((CC(M1,I-1,1,K)+TR12*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))&
     &  +TR11*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))+(TI12*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))-TI11*(CC(M1,I,5,K)+CC(M1,IC,4,K))))&
     &  -WA3(I-1)&
     & *((CC(M1,I,1,K)+TR12*(CC(M1,I,3,K)-&
     &  CC(M1,IC,2,K))+TR11*(CC(M1,I,5,K)-CC(M1,IC,4,K)))&
     &  -(TI12*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))-TI11&
     &  *(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))
        CH(M2,I,K,4) = WA3(I-2)&
     & *((CC(M1,I,1,K)+TR12*(CC(M1,I,3,K)-&
     &  CC(M1,IC,2,K))+TR11*(CC(M1,I,5,K)-CC(M1,IC,4,K)))&
     &  -(TI12*(CC(M1,I-1,3,K)-CC(M1,IC-1,2,K))-TI11&
     &  *(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))&
     &  +WA3(I-1)&
     &  *((CC(M1,I-1,1,K)+TR12*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))&
     &  +TR11*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))+(TI12*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))-TI11*(CC(M1,I,5,K)+CC(M1,IC,4,K))))
        CH(M2,I-1,K,5) = WA4(I-2)&
     &  *((CC(M1,I-1,1,K)+TR11*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))&
     &  +TR12*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))+(TI11*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))+TI12*(CC(M1,I,5,K)+CC(M1,IC,4,K))))&
     &  -WA4(I-1)&
     &  *((CC(M1,I,1,K)+TR11*(CC(M1,I,3,K)-CC(M1,IC,2,K))&
     &  +TR12*(CC(M1,I,5,K)-CC(M1,IC,4,K)))-(TI11*(CC(M1,I-1,3,K)&
     &  -CC(M1,IC-1,2,K))+TI12*(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))
        CH(M2,I,K,5) = WA4(I-2)&
     &  *((CC(M1,I,1,K)+TR11*(CC(M1,I,3,K)-CC(M1,IC,2,K))&
     &  +TR12*(CC(M1,I,5,K)-CC(M1,IC,4,K)))-(TI11*(CC(M1,I-1,3,K)&
     &  -CC(M1,IC-1,2,K))+TI12*(CC(M1,I-1,5,K)-CC(M1,IC-1,4,K))))&
     &  +WA4(I-1)&
     &  *((CC(M1,I-1,1,K)+TR11*(CC(M1,I-1,3,K)+CC(M1,IC-1,2,K))&
     &  +TR12*(CC(M1,I-1,5,K)+CC(M1,IC-1,4,K)))+(TI11*(CC(M1,I,3,K)&
     &  +CC(M1,IC,2,K))+TI12*(CC(M1,I,5,K)+CC(M1,IC,4,K))))
 1002      CONTINUE
  102    CONTINUE
  103 CONTINUE
      END subroutine MRADB5
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADBG (M,IDO,IP,L1,IDL1,CC,C1,C2,IM1,IN1,CH,CH2,IM2,IN2,WA)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: NBD,IK
      real(kind=rk)     :: AR1,AI1,AR2,AI2
      real(kind=rk)     :: DC2,DCP,DS2,DSP
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TPI           
      REAL(kind=rk)  :: CH(IN2,IDO,L1,IP),CC(IN1,IDO,IP,L1),C1(IN1,IDO,L1,IP),C2(IN1,IDL1,IP),CH2(IN2,IDL1,IP),WA(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      !!TPI=2._dk*4._dk*ATAN(1.0_dk)
      TPI=2._dk*PI_long
      ARG = TPI/REAL(IP,kind=rk)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            M2 = M2S
            DO 1001 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,I,K,1) = CC(M1,I,1,K)
 1001       CONTINUE
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            M2 = M2S
            DO 1004 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,I,K,1) = CC(M1,I,1,K)
 1004       CONTINUE
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            M2 = M2S
            DO 1007 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,1,K,J) = CC(M1,IDO,J2-2,K)+CC(M1,IDO,J2-2,K)
            CH(M2,1,K,JC) = CC(M1,1,J2-1,K)+CC(M1,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               M2 = M2S
               DO 1009 M1=1,M1D,IM1
               M2 = M2+IM2
               CH(M2,I-1,K,J) = CC(M1,I-1,2*J-1,K)+CC(M1,IC-1,2*J-2,K)
               CH(M2,I-1,K,JC) = CC(M1,I-1,2*J-1,K)-CC(M1,IC-1,2*J-2,K)
               CH(M2,I,K,J) = CC(M1,I,2*J-1,K)-CC(M1,IC,2*J-2,K)
               CH(M2,I,K,JC) = CC(M1,I,2*J-1,K)+CC(M1,IC,2*J-2,K)
 1009          CONTINUE
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               M2 = M2S
               DO 1013 M1=1,M1D,IM1
               M2 = M2+IM2
               CH(M2,I-1,K,J) = CC(M1,I-1,2*J-1,K)+CC(M1,IC-1,2*J-2,K)
               CH(M2,I-1,K,JC) = CC(M1,I-1,2*J-1,K)-CC(M1,IC-1,2*J-2,K)
               CH(M2,I,K,J) = CC(M1,I,2*J-1,K)-CC(M1,IC,2*J-2,K)
               CH(M2,I,K,JC) = CC(M1,I,2*J-1,K)+CC(M1,IC,2*J-2,K)
 1013          CONTINUE
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1._rk
      AI1 = 0._rk
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            M2 = M2S
            DO 1017 M1=1,M1D,IM1
            M2 = M2+IM2
            C2(M1,IK,L) = CH2(M2,IK,1)+AR1*CH2(M2,IK,2)
            C2(M1,IK,LC) = AI1*CH2(M2,IK,IP)
 1017       CONTINUE
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               M2 = M2S
               DO 1018 M1=1,M1D,IM1
               M2 = M2+IM2
               C2(M1,IK,L) = C2(M1,IK,L)+AR2*CH2(M2,IK,J)
               C2(M1,IK,LC) = C2(M1,IK,LC)+AI2*CH2(M2,IK,JC)
 1018          CONTINUE
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            M2 = M2S
            DO 1021 M1=1,M1D,IM1
            M2 = M2+IM2
            CH2(M2,IK,1) = CH2(M2,IK,1)+CH2(M2,IK,J)
 1021       CONTINUE
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            M2 = M2S
            DO 1023 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,1,K,J) = C1(M1,1,K,J)-C1(M1,1,K,JC)
            CH(M2,1,K,JC) = C1(M1,1,K,J)+C1(M1,1,K,JC)
 1023       CONTINUE
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               M2 = M2S
               DO 1025 M1=1,M1D,IM1
               M2 = M2+IM2
               CH(M2,I-1,K,J) = C1(M1,I-1,K,J)-C1(M1,I,K,JC)
               CH(M2,I-1,K,JC) = C1(M1,I-1,K,J)+C1(M1,I,K,JC)
               CH(M2,I,K,J) = C1(M1,I,K,J)+C1(M1,I-1,K,JC)
               CH(M2,I,K,JC) = C1(M1,I,K,J)-C1(M1,I-1,K,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               M2 = M2S
               DO 1029 M1=1,M1D,IM1
               M2 = M2+IM2
               CH(M2,I-1,K,J) = C1(M1,I-1,K,J)-C1(M1,I,K,JC)
               CH(M2,I-1,K,JC) = C1(M1,I-1,K,J)+C1(M1,I,K,JC)
               CH(M2,I,K,J) = C1(M1,I,K,J)+C1(M1,I-1,K,JC)
               CH(M2,I,K,JC) = C1(M1,I,K,J)-C1(M1,I-1,K,JC)
 1029          CONTINUE
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         M2 = M2S
         DO 1033 M1=1,M1D,IM1
         M2 = M2+IM2
         C2(M1,IK,1) = CH2(M2,IK,1)
 1033    CONTINUE
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            M2 = M2S
            DO 1034 M1=1,M1D,IM1
            M2 = M2+IM2
            C1(M1,1,K,J) = CH(M2,1,K,J)
 1034       CONTINUE
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               M2 = M2S
               DO 1036 M1=1,M1D,IM1
               M2 = M2+IM2
               C1(M1,I-1,K,J) = WA(IDIJ-1)*CH(M2,I-1,K,J)-WA(IDIJ)*CH(M2,I,K,J)
               C1(M1,I,K,J) = WA(IDIJ-1)*CH(M2,I,K,J)+WA(IDIJ)*CH(M2,I-1,K,J)
 1036          CONTINUE
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               M2 = M2S
               DO 1040 M1=1,M1D,IM1
               M2 = M2+IM2
               C1(M1,I-1,K,J) = WA(IDIJ-1)*CH(M2,I-1,K,J)-WA(IDIJ)*CH(M2,I,K,J)
               C1(M1,I,K,J) = WA(IDIJ-1)*CH(M2,I,K,J)+WA(IDIJ)*CH(M2,I-1,K,J)
 1040          CONTINUE
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 continue
      END subroutine MRADBG
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADF2 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      REAL(kind=rk)     :: CH(IN2,IDO,2,L1),CC(IN1,IDO,L1,2), WA1(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      DO 101 K=1,L1
         M2 = M2S
         DO 1001 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(M2,1,1,K) = CC(M1,1,K,1)+CC(M1,1,K,2)
         CH(M2,IDO,2,K) = CC(M1,1,K,1)-CC(M1,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            M2 = M2S
            DO 1003 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,I,1,K) = CC(M1,I,K,1)+(WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2))
            CH(M2,IC,2,K) = (WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2))-CC(M1,I,K,1)
            CH(M2,I-1,1,K) = CC(M1,I-1,K,1)+(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2))
            CH(M2,IC-1,2,K) = CC(M1,I-1,K,1)-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         M2 = M2S
         DO 1006 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(M2,1,2,K) = -CC(M1,IDO,K,2)
         CH(M2,IDO,1,K) = CC(M1,IDO,K,1)
 1006    CONTINUE
  106 CONTINUE
  107 continue
      END subroutine MRADF2
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADF3 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1,WA2)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TAUR,TAUI
      REAL(kind=rk)     :: CH(IN2,IDO,3,L1),CC(IN1,IDO,L1,3),WA1(IDO),WA2(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/3._dk
      ARG = 2._dk*PI_long/3._dk
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         M2 = M2S
         DO 1001 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(M2,1,1,K) = CC(M1,1,K,1)+(CC(M1,1,K,2)+CC(M1,1,K,3))
         CH(M2,1,3,K) = TAUI*(CC(M1,1,K,3)-CC(M1,1,K,2))
         CH(M2,IDO,2,K) = CC(M1,1,K,1)+TAUR*(CC(M1,1,K,2)+CC(M1,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            M2 = M2S
            DO 1002 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,I-1,1,K) = CC(M1,I-1,K,1)+((WA1(I-2)*CC(M1,I-1,K,2)+&
     &       WA1(I-1)*CC(M1,I,K,2))+(WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3)))
            CH(M2,I,1,K) = CC(M1,I,K,1)+((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3)))
            CH(M2,I-1,3,K) = (CC(M1,I-1,K,1)+TAUR*((WA1(I-2)*&
     &       CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2))+(WA2(I-2)*&
     &       CC(M1,I-1,K,3)+WA2(I-1)*CC(M1,I,K,3))))+(TAUI*((WA1(I-2)*&
     &       CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2))-(WA2(I-2)*&
     &       CC(M1,I,K,3)-WA2(I-1)*CC(M1,I-1,K,3))))
            CH(M2,IC-1,2,K) = (CC(M1,I-1,K,1)+TAUR*((WA1(I-2)*&
     &       CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2))+(WA2(I-2)*&
     &       CC(M1,I-1,K,3)+WA2(I-1)*CC(M1,I,K,3))))-(TAUI*((WA1(I-2)*&
     &       CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2))-(WA2(I-2)*&
     &       CC(M1,I,K,3)-WA2(I-1)*CC(M1,I-1,K,3))))
            CH(M2,I,3,K) = (CC(M1,I,K,1)+TAUR*((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2))))
            CH(M2,IC,2,K) = (TAUI*((WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2))))-(CC(M1,I,K,1)+TAUR*((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      END subroutine MRADF3
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADF4 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1,WA2,WA3)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: HSQT2
      REAL(kind=rk)     :: CC(IN1,IDO,L1,4),CH(IN2,IDO,4,L1),WA1(IDO),WA2(IDO),WA3(IDO)
!
      HSQT2=SQRT(2._rk)/2._rk
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      DO 101 K=1,L1
         M2 = M2S
         DO 1001 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(M2,1,1,K) = (CC(M1,1,K,2)+CC(M1,1,K,4))+(CC(M1,1,K,1)+CC(M1,1,K,3))
         CH(M2,IDO,4,K) = (CC(M1,1,K,1)+CC(M1,1,K,3))-(CC(M1,1,K,2)+CC(M1,1,K,4))
         CH(M2,IDO,2,K) = CC(M1,1,K,1)-CC(M1,1,K,3)
         CH(M2,1,3,K) = CC(M1,1,K,4)-CC(M1,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            M2 = M2S
            DO 1003 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,I-1,1,K) = ((WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2))+(WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4)))+(CC(M1,I-1,K,1)+(WA2(I-2)*CC(M1,I-1,K,3)+&
     &       WA2(I-1)*CC(M1,I,K,3)))
            CH(M2,IC-1,4,K) = (CC(M1,I-1,K,1)+(WA2(I-2)*CC(M1,I-1,K,3)+&
     &       WA2(I-1)*CC(M1,I,K,3)))-((WA1(I-2)*CC(M1,I-1,K,2)+&
     &       WA1(I-1)*CC(M1,I,K,2))+(WA3(I-2)*CC(M1,I-1,K,4)+&
     &       WA3(I-1)*CC(M1,I,K,4)))
            CH(M2,I,1,K) = ((WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*&
     &       CC(M1,I-1,K,2))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4)))+(CC(M1,I,K,1)+(WA2(I-2)*CC(M1,I,K,3)-&
     &       WA2(I-1)*CC(M1,I-1,K,3)))
            CH(M2,IC,4,K) = ((WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*&
     &       CC(M1,I-1,K,2))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4)))-(CC(M1,I,K,1)+(WA2(I-2)*CC(M1,I,K,3)-&
     &       WA2(I-1)*CC(M1,I-1,K,3)))
            CH(M2,I-1,3,K) = ((WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*&
     &       CC(M1,I-1,K,2))-(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4)))+(CC(M1,I-1,K,1)-(WA2(I-2)*CC(M1,I-1,K,3)+&
     &       WA2(I-1)*CC(M1,I,K,3)))
            CH(M2,IC-1,2,K) = (CC(M1,I-1,K,1)-(WA2(I-2)*CC(M1,I-1,K,3)+&
     &       WA2(I-1)*CC(M1,I,K,3)))-((WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*&
     &       CC(M1,I-1,K,2))-(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4)))
            CH(M2,I,3,K) = ((WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2)))+(CC(M1,I,K,1)-(WA2(I-2)*CC(M1,I,K,3)-&
     &       WA2(I-1)*CC(M1,I-1,K,3)))
            CH(M2,IC,2,K) = ((WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2)))-(CC(M1,I,K,1)-(WA2(I-2)*CC(M1,I,K,3)-&
     &       WA2(I-1)*CC(M1,I-1,K,3)))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         M2 = M2S
         DO 1006 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,IDO,1,K) = (HSQT2*(CC(M1,IDO,K,2)-CC(M1,IDO,K,4)))+CC(M1,IDO,K,1)
            CH(M2,IDO,3,K) = CC(M1,IDO,K,1)-(HSQT2*(CC(M1,IDO,K,2)-CC(M1,IDO,K,4)))
            CH(M2,1,2,K) = (-HSQT2*(CC(M1,IDO,K,2)+CC(M1,IDO,K,4)))-CC(M1,IDO,K,3)
            CH(M2,1,4,K) = (-HSQT2*(CC(M1,IDO,K,2)+CC(M1,IDO,K,4)))+CC(M1,IDO,K,3)
 1006    CONTINUE
  106 CONTINUE
  107 continue
      END subroutine MRADF4
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADF5 (M,IDO,L1,CC,IM1,IN1,CH,IM2,IN2,WA1,WA2,WA3,WA4)
      integer           :: IM1,IM2,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      REAL(kind=rk)     :: CC(IN1,IDO,L1,5),CH(IN2,IDO,5,L1),WA1(IDO),WA2(IDO),WA3(IDO),WA4(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/5._dk
      ARG = 2._dk*PI_long/5._dk
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2._rk*ARG)
      TI12=SIN(2._rk*ARG)
      DO 101 K=1,L1
         M2 = M2S
         DO 1001 M1=1,M1D,IM1
         M2 = M2+IM2
         CH(M2,1,1,K) = CC(M1,1,K,1)+(CC(M1,1,K,5)+CC(M1,1,K,2))+(CC(M1,1,K,4)+CC(M1,1,K,3))
         CH(M2,IDO,2,K) = CC(M1,1,K,1)+TR11*(CC(M1,1,K,5)+CC(M1,1,K,2))+TR12*(CC(M1,1,K,4)+CC(M1,1,K,3))
         CH(M2,1,3,K) = TI11*(CC(M1,1,K,5)-CC(M1,1,K,2))+TI12*(CC(M1,1,K,4)-CC(M1,1,K,3))
         CH(M2,IDO,4,K) = CC(M1,1,K,1)+TR12*(CC(M1,1,K,5)+CC(M1,1,K,2))+TR11*(CC(M1,1,K,4)+CC(M1,1,K,3))
         CH(M2,1,5,K) = TI12*(CC(M1,1,K,5)-CC(M1,1,K,2))-TI11*(CC(M1,1,K,4)-CC(M1,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            M2 = M2S
            DO 1002 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,I-1,1,K) = CC(M1,I-1,K,1)+((WA1(I-2)*CC(M1,I-1,K,2)+&
     &       WA1(I-1)*CC(M1,I,K,2))+(WA4(I-2)*CC(M1,I-1,K,5)+WA4(I-1)*&
     &       CC(M1,I,K,5)))+((WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))+(WA3(I-2)*CC(M1,I-1,K,4)+&
     &       WA3(I-1)*CC(M1,I,K,4)))
            CH(M2,I,1,K) = CC(M1,I,K,1)+((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*&
     &       CC(M1,I-1,K,5)))+((WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4)))
            CH(M2,I-1,3,K) = CC(M1,I-1,K,1)+TR11*&
     &      ( WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2)&
     &       +WA4(I-2)*CC(M1,I-1,K,5)+WA4(I-1)*CC(M1,I,K,5))+TR12*&
     &      ( WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*CC(M1,I,K,3)&
     &       +WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*CC(M1,I,K,4))+TI11*&
     &      ( WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2)&
     &       -(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*CC(M1,I-1,K,5)))+TI12*&
     &      ( WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*CC(M1,I-1,K,3)&
     &       -(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*CC(M1,I-1,K,4)))
            CH(M2,IC-1,2,K) = CC(M1,I-1,K,1)+TR11*&
     &      ( WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2)&
     &       +WA4(I-2)*CC(M1,I-1,K,5)+WA4(I-1)*CC(M1,I,K,5))+TR12*&
     &     ( WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*CC(M1,I,K,3)&
     &      +WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*CC(M1,I,K,4))-(TI11*&
     &      ( WA1(I-2)*CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2)&
     &       -(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*CC(M1,I-1,K,5)))+TI12*&
     &      ( WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*CC(M1,I-1,K,3)&
     &       -(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*CC(M1,I-1,K,4))))
            CH(M2,I,3,K) = (CC(M1,I,K,1)+TR11*((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*&
     &       CC(M1,I-1,K,5)))+TR12*((WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M1,I-1,K,5)+&
     &       WA4(I-1)*CC(M1,I,K,5))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2)))+TI12*((WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4))-(WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))))
            CH(M2,IC,2,K) = (TI11*((WA4(I-2)*CC(M1,I-1,K,5)+WA4(I-1)*&
     &       CC(M1,I,K,5))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2)))+TI12*((WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4))-(WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))))-(CC(M1,I,K,1)+TR11*((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*&
     &       CC(M1,I-1,K,5)))+TR12*((WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4))))
            CH(M2,I-1,5,K) = (CC(M1,I-1,K,1)+TR12*((WA1(I-2)*&
     &       CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2))+(WA4(I-2)*&
     &       CC(M1,I-1,K,5)+WA4(I-1)*CC(M1,I,K,5)))+TR11*((WA2(I-2)*&
     &       CC(M1,I-1,K,3)+WA2(I-1)*CC(M1,I,K,3))+(WA3(I-2)*&
     &       CC(M1,I-1,K,4)+WA3(I-1)*CC(M1,I,K,4))))+(TI12*((WA1(I-2)*&
     &       CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2))-(WA4(I-2)*&
     &       CC(M1,I,K,5)-WA4(I-1)*CC(M1,I-1,K,5)))-TI11*((WA2(I-2)*&
     &       CC(M1,I,K,3)-WA2(I-1)*CC(M1,I-1,K,3))-(WA3(I-2)*&
     &       CC(M1,I,K,4)-WA3(I-1)*CC(M1,I-1,K,4))))
            CH(M2,IC-1,4,K) = (CC(M1,I-1,K,1)+TR12*((WA1(I-2)*&
     &       CC(M1,I-1,K,2)+WA1(I-1)*CC(M1,I,K,2))+(WA4(I-2)*&
     &       CC(M1,I-1,K,5)+WA4(I-1)*CC(M1,I,K,5)))+TR11*((WA2(I-2)*&
     &       CC(M1,I-1,K,3)+WA2(I-1)*CC(M1,I,K,3))+(WA3(I-2)*&
     &       CC(M1,I-1,K,4)+WA3(I-1)*CC(M1,I,K,4))))-(TI12*((WA1(I-2)*&
     &       CC(M1,I,K,2)-WA1(I-1)*CC(M1,I-1,K,2))-(WA4(I-2)*&
     &       CC(M1,I,K,5)-WA4(I-1)*CC(M1,I-1,K,5)))-TI11*((WA2(I-2)*&
     &       CC(M1,I,K,3)-WA2(I-1)*CC(M1,I-1,K,3))-(WA3(I-2)*&
     &       CC(M1,I,K,4)-WA3(I-1)*CC(M1,I-1,K,4))))
            CH(M2,I,5,K) = (CC(M1,I,K,1)+TR12*((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*&
     &       CC(M1,I-1,K,5)))+TR11*((WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M1,I-1,K,5)+&
     &       WA4(I-1)*CC(M1,I,K,5))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2)))-TI11*((WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4))-(WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))))
            CH(M2,IC,4,K) = (TI12*((WA4(I-2)*CC(M1,I-1,K,5)+WA4(I-1)*&
     &       CC(M1,I,K,5))-(WA1(I-2)*CC(M1,I-1,K,2)+WA1(I-1)*&
     &       CC(M1,I,K,2)))-TI11*((WA3(I-2)*CC(M1,I-1,K,4)+WA3(I-1)*&
     &       CC(M1,I,K,4))-(WA2(I-2)*CC(M1,I-1,K,3)+WA2(I-1)*&
     &       CC(M1,I,K,3))))-(CC(M1,I,K,1)+TR12*((WA1(I-2)*CC(M1,I,K,2)-&
     &       WA1(I-1)*CC(M1,I-1,K,2))+(WA4(I-2)*CC(M1,I,K,5)-WA4(I-1)*&
     &       CC(M1,I-1,K,5)))+TR11*((WA2(I-2)*CC(M1,I,K,3)-WA2(I-1)*&
     &       CC(M1,I-1,K,3))+(WA3(I-2)*CC(M1,I,K,4)-WA3(I-1)*&
     &       CC(M1,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      END subroutine MRADF5
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRADFG (M,IDO,IP,L1,IDL1,CC,C1,C2,IM1,IN1,CH,CH2,IM2,IN2,WA)
      integer           :: N,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: NBD,IK
      real(kind=rk)     :: AR1,AI1,AR2,AI2
      real(kind=rk)     :: DC2,DCP,DS2,DSP
      real(kind=rk)     :: TPI           
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      REAL(kind=rk)     :: CH(IN2,IDO,L1,IP),CC(IN1,IDO,IP,L1),C1(IN1,IDO,L1,IP),C2(IN1,IDL1,IP),CH2(IN2,IDL1,IP),WA(IDO)
!
      M1D = (M-1)*IM1+1
      M2S = 1-IM2
      !!TPI=2._dk*4._dk*ATAN(1.0_dk)
      TPI=2._dk*PI_long
      ARG = TPI/REAL(IP,kind=rk)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         M2 = M2S
         DO 1001 M1=1,M1D,IM1
         M2 = M2+IM2
         CH2(M2,IK,1) = C2(M1,IK,1)
 1001    CONTINUE
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            M2 = M2S
            DO 1002 M1=1,M1D,IM1
            M2 = M2+IM2
            CH(M2,1,K,J) = C1(M1,1,K,J)
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               M2 = M2S
               DO 1004 M1=1,M1D,IM1
               M2 = M2+IM2
               CH(M2,I-1,K,J) = WA(IDIJ-1)*C1(M1,I-1,K,J)+WA(IDIJ)*C1(M1,I,K,J)
               CH(M2,I,K,J) = WA(IDIJ-1)*C1(M1,I,K,J)-WA(IDIJ)*C1(M1,I-1,K,J)
 1004          CONTINUE
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               M2 = M2S
               DO 1008 M1=1,M1D,IM1
               M2 = M2+IM2
               CH(M2,I-1,K,J) = WA(IDIJ-1)*C1(M1,I-1,K,J)+WA(IDIJ)*C1(M1,I,K,J)
               CH(M2,I,K,J) = WA(IDIJ-1)*C1(M1,I,K,J)-WA(IDIJ)*C1(M1,I-1,K,J)
 1008          CONTINUE
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               M2 = M2S
               DO 1012 M1=1,M1D,IM1
               M2 = M2+IM2
               C1(M1,I-1,K,J) = CH(M2,I-1,K,J)+CH(M2,I-1,K,JC)
               C1(M1,I-1,K,JC) = CH(M2,I,K,J)-CH(M2,I,K,JC)
               C1(M1,I,K,J) = CH(M2,I,K,J)+CH(M2,I,K,JC)
               C1(M1,I,K,JC) = CH(M2,I-1,K,JC)-CH(M2,I-1,K,J)
 1012          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               M2 = M2S
               DO 1016 M1=1,M1D,IM1
               M2 = M2+IM2
               C1(M1,I-1,K,J) = CH(M2,I-1,K,J)+CH(M2,I-1,K,JC)
               C1(M1,I-1,K,JC) = CH(M2,I,K,J)-CH(M2,I,K,JC)
               C1(M1,I,K,J) = CH(M2,I,K,J)+CH(M2,I,K,JC)
               C1(M1,I,K,JC) = CH(M2,I-1,K,JC)-CH(M2,I-1,K,J)
 1016          CONTINUE
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         M2 = M2S
         DO 1020 M1=1,M1D,IM1
         M2 = M2+IM2
         C2(M1,IK,1) = CH2(M2,IK,1)
 1020    CONTINUE
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            M2 = M2S
            DO 1022 M1=1,M1D,IM1
            M2 = M2+IM2
            C1(M1,1,K,J) = CH(M2,1,K,J)+CH(M2,1,K,JC)
            C1(M1,1,K,JC) = CH(M2,1,K,JC)-CH(M2,1,K,J)
 1022       CONTINUE
  122    CONTINUE
  123 CONTINUE
!
      AR1 = 1._rk
      AI1 = 0._rk
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            M2 = M2S
            DO 1024 M1=1,M1D,IM1
            M2 = M2+IM2
            CH2(M2,IK,L) = C2(M1,IK,1)+AR1*C2(M1,IK,2)
            CH2(M2,IK,LC) = AI1*C2(M1,IK,IP)
 1024       CONTINUE
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               M2 = M2S
               DO 1025 M1=1,M1D,IM1
               M2 = M2+IM2
               CH2(M2,IK,L) = CH2(M2,IK,L)+AR2*C2(M1,IK,J)
               CH2(M2,IK,LC) = CH2(M2,IK,LC)+AI2*C2(M1,IK,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            M2 = M2S
            DO 1028 M1=1,M1D,IM1
            M2 = M2+IM2
            CH2(M2,IK,1) = CH2(M2,IK,1)+C2(M1,IK,J)
 1028       CONTINUE
  128    CONTINUE
  129 CONTINUE
!
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            M2 = M2S
            DO 1030 M1=1,M1D,IM1
            M2 = M2+IM2
            CC(M1,I,1,K) = CH(M2,I,K,1)
 1030       CONTINUE
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            M2 = M2S
            DO 1033 M1=1,M1D,IM1
            M2 = M2+IM2
            CC(M1,I,1,K) = CH(M2,I,K,1)
 1033       CONTINUE
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            M2 = M2S
            DO 1036 M1=1,M1D,IM1
            M2 = M2+IM2
            CC(M1,IDO,J2-2,K) = CH(M2,1,K,J)
            CC(M1,1,J2-1,K) = CH(M2,1,K,JC)
 1036       CONTINUE
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               M2 = M2S
               DO 1038 M1=1,M1D,IM1
               M2 = M2+IM2
               CC(M1,I-1,J2-1,K) = CH(M2,I-1,K,J)+CH(M2,I-1,K,JC)
               CC(M1,IC-1,J2-2,K) = CH(M2,I-1,K,J)-CH(M2,I-1,K,JC)
               CC(M1,I,J2-1,K) = CH(M2,I,K,J)+CH(M2,I,K,JC)
               CC(M1,IC,J2-2,K) = CH(M2,I,K,JC)-CH(M2,I,K,J)
 1038          CONTINUE
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               M2 = M2S
               DO 1042 M1=1,M1D,IM1
               M2 = M2+IM2
               CC(M1,I-1,J2-1,K) = CH(M2,I-1,K,J)+CH(M2,I-1,K,JC)
               CC(M1,IC-1,J2-2,K) = CH(M2,I-1,K,J)-CH(M2,I-1,K,JC)
               CC(M1,I,J2-1,K) = CH(M2,I,K,J)+CH(M2,I,K,JC)
               CC(M1,IC,J2-2,K) = CH(M2,I,K,JC)-CH(M2,I,K,J)
 1042          CONTINUE
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      END subroutine MRADFG
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRFTB1 (M,IM,N,IN,C,CH,WA,FAC)
      integer           :: N,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IX2,IX3,IX4
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: HALF,HALFM 
      REAL(kind=rk)     :: CH(M,*), C(IN,*), WA(N)
      integer           :: FAC(*)
!
      NF = FAC(2)
      NA = 0
      DO 10 K1=1,NF
      IP = FAC(K1+2)
      NA = 1-NA
      IF(IP .LE. 5) GO TO 10
      IF(K1 .EQ. NF) GO TO 10
      NA = 1-NA
   10 CONTINUE
      HALF = .5_rk
      HALFM = -.5_rk
      MODN = MOD(N,2)
      NL = N-2
      IF(MODN .NE. 0) NL = N-1
      IF (NA .EQ. 0) GO TO 120
      M2 = 1-IM
      DO 117 I=1,M
      M2 = M2+IM
      CH(I,1) = C(M2,1)
      CH(I,N) = C(M2,N)
  117 CONTINUE
      DO 118 J=2,NL,2
      M2 = 1-IM
      DO 118 I=1,M
         M2 = M2+IM
        CH(I,J) = HALF*C(M2,J)
        CH(I,J+1) = HALFM*C(M2,J+1)
  118 CONTINUE
      GO TO 124
  120 CONTINUE
      DO 122 J=2,NL,2
      M2 = 1-IM
      DO 122 I=1,M
         M2 = M2+IM
        C(M2,J) = HALF*C(M2,J)
        C(M2,J+1) = HALFM*C(M2,J+1)
  122 CONTINUE
  124 L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL MRADB4 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL MRADB4 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
        CALL MRADB2 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW))
         GO TO 105
  104    CALL MRADB2 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
        CALL MRADB3 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW),WA(IX2))
         GO TO 108
  107    CALL MRADB3 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL MRADB5 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL MRADB5 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
        CALL MRADBG (M,IDO,IP,L1,IDL1,C,C,C,IM,IN,CH,CH,1,M,WA(IW))
         GO TO 114
  113    CALL MRADBG (M,IDO,IP,L1,IDL1,CH,CH,CH,1,M,C,C,IM,IN,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      END subroutine MRFTB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRFTF1 (M,IM,N,IN,C,CH,WA,FAC)
      integer           :: N,IM1,IM2,IN,IN1,IN2
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IX2,IX3,IX4
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: TSN,SN,TSNM
      integer           :: FAC(15)
      REAL(kind=rk)     :: CH(M,*),C(IN,*),WA(N)
!
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
        CALL MRADF4 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL MRADF4 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
        CALL MRADF2 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW))
         GO TO 110
  103    CALL MRADF2 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
        CALL MRADF3 (M,IDO,L1,C,IM,IN,CH,1,M,WA(IW),WA(IX2))
         GO TO 110
  105    CALL MRADF3 (M,IDO,L1,CH,1,M,C,IM,IN,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL MRADF5(M,IDO,L1,C,IM,IN,CH,1,M,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL MRADF5(M,IDO,L1,CH,1,M,C,IM,IN,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
        CALL MRADFG (M,IDO,IP,L1,IDL1,C,C,C,IM,IN,CH,CH,1,M,WA(IW))
         NA = 1
         GO TO 110
  109    CALL MRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,1,M,C,C,IM,IN,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SN = 1._rk/real(N,kind=rk)
      TSN = 2._rk/real(N,kind=rk)
      TSNM = -TSN
      MODN = MOD(N,2)
      NL = N-2
      IF(MODN .NE. 0) NL = N-1
      IF (NA .NE. 0) GO TO 120
      M2 = 1-IM
      DO 117 I=1,M
         M2 = M2+IM
         C(M2,1) = SN*CH(I,1)
  117 CONTINUE
      DO 118 J=2,NL,2
      M2 = 1-IM
      DO 118 I=1,M
         M2 = M2+IM
        C(M2,J) = TSN*CH(I,J)
        C(M2,J+1) = TSNM*CH(I,J+1)
  118 CONTINUE
      IF(MODN .NE. 0) RETURN
      M2 = 1-IM
      DO 119 I=1,M
         M2 = M2+IM
         C(M2,N) = SN*CH(I,N)
  119 CONTINUE
      RETURN
  120 M2 = 1-IM
      DO 121 I=1,M
         M2 = M2+IM
         C(M2,1) = SN*C(M2,1)
  121 CONTINUE
      DO 122 J=2,NL,2
      M2 = 1-IM
      DO 122 I=1,M
         M2 = M2+IM
        C(M2,J) = TSN*C(M2,J)
        C(M2,J+1) = TSNM*C(M2,J+1)
  122 CONTINUE
      IF(MODN .NE. 0) RETURN
      M2 = 1-IM
      DO 123 I=1,M
         M2 = M2+IM
         C(M2,N) = SN*C(M2,N)
  123 CONTINUE
      END subroutine MRFTF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE MRFTI1 (N,WA,FAC)
      ! no difference to routine RFFTI1
      integer           :: N,IN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: NTRY
      integer           :: IPM
      REAL(kind=rk)     :: WA(N),FAC(15)
      real(kind=dk)     :: FI
      INTEGER    NTRYH(4)
      real(kind=dk)    :: TPI,ARGH,ARGLD,ARG
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
!
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      !!TPI = 8._dk*DATAN(1._dk)
      TPI=2._dk*PI_long
      ARGH = TPI/REAL(N,kind=dk)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = REAL(LD,kind=dk)*ARGH
            FI = 0._dk
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1._dk
               ARG = FI*ARGLD
              !!WA(I-1) = DCOS(ARG)
              !!WA(I) = DSIN(ARG)
              WA(I-1) = COS(ARG)
              WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      END subroutine MRFTI1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MSNTB1(LOT,JUMP,N,INC,X,WSAVE,DSUM,XH,WORK, fft_sign,IER)
      SUBROUTINE MSNTB1(LOT,JUMP,N,INC,X,DSUM,XH,WORK, fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: IER1
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,XHOLD,FNP1S4,DT,PI
      !!real(kind=rk)     :: WSAVE(*)
      real(kind=rk)     :: X(INC,*),XH(LOT,*)
      real(kind=rk)     :: WORK(*)   ! Kuester
      real(kind=rk)     :: T1,T2
      real(kind=dk)    :: DSUM(*) ! no value on entrance
      type(fft_real_sign_type) :: fft_sign
      IER = 0
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 200,102,103
 102  SRT3S2 = SQRT(3._rk)/2._rk
      DO 112 M=1,LJ,JUMP
         XHOLD = SRT3S2*(X(M,1)+X(M,2))
         X(M,2) = SRT3S2*(X(M,1)-X(M,2))
         X(M,1) = XHOLD
  112 CONTINUE
      GO TO 200
  103 NP1 = N+1
      NS2 = N/2
      DO 104 K=1,NS2
         KC = NP1-K
         M1 = 0
         DO 114 M=1,LJ,JUMP
         M1 = M1+1
         T1 = X(M,K)-X(M,KC)
         T2 = fft_sign%WSAVE(K)*(X(M,K)+X(M,KC))
         XH(M1,K+1) = T1+T2
         XH(M1,KC+1) = T2-T1
  114    CONTINUE
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
      M1 = 0
      DO 123 M=1,LJ,JUMP
         M1 = M1+1
         XH(M1,NS2+2) = 4._rk*X(M,NS2+1)
  123 CONTINUE
  124 DO 127 M=1,LOT
         XH(M,1) = 0._rk
  127 CONTINUE
      LNXH = LOT-1 + LOT*(NP1-1) + 1
      LNSV = NP1 + INT(LOG(REAL(NP1))/LOG(2._rk)) + 4
      LNWK = LOT*NP1
!
      !!CALL RFFTMF(LOT,1,NP1,LOT,XH,LNXH,WSAVE(NS2+1),LNSV,WORK,LNWK, fft_sign,IER1)
      CALL RFFTMF(LOT,1,NP1,LOT,XH,LNXH,WORK,LNWK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('MSNTB1',-5)
        GO TO 200
      ENDIF
!
      IF(MOD(NP1,2) .NE. 0) GO TO 30
      DO 20 M=1,LOT
      XH(M,NP1) = XH(M,NP1)+XH(M,NP1)
   20 CONTINUE
 30   FNP1S4 = REAL(NP1,kind=rk)/4._rk
      M1 = 0
      DO 125 M=1,LJ,JUMP
         M1 = M1+1
         X(M,1) = FNP1S4*XH(M1,1)
         DSUM(M1) = X(M,1)
  125 CONTINUE
      DO 105 I=3,N,2
         M1 = 0
         DO 115 M=1,LJ,JUMP
            M1 = M1+1
            X(M,I-1) = FNP1S4*XH(M1,I)
            DSUM(M1) = DSUM(M1)+FNP1S4*XH(M1,I-1)
            X(M,I) = DSUM(M1)
  115    CONTINUE
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 200
      M1 = 0
      DO 116 M=1,LJ,JUMP
         M1 = M1+1
         X(M,N) = FNP1S4*XH(M1,N+1)
  116 CONTINUE
!
  200 CONTINUE
      END subroutine MSNTB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE MSNTF1(LOT,JUMP,N,INC,X,WSAVE,DSUM,XH,WORK, fft_sign,IER)
      SUBROUTINE MSNTF1(LOT,JUMP,N,INC,X,DSUM,XH,WORK, fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: IER1
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,XHOLD,FNP1S4,DT,PI
      !!real(kind=rk)     :: WSAVE(*)
      real(kind=rk)     :: X(INC,*),XH(LOT,*)
      real(kind=rk)     :: WORK(*)  ! Kuester
      real(kind=rk)     :: T1,T2
      real(kind=dk)    :: DSUM(*)
      type(fft_real_sign_type) :: fft_sign
      IER = 0
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 101,102,103
 102  SSQRT3 = 1._rk/SQRT(3._rk)
      DO 112 M=1,LJ,JUMP
         XHOLD = SSQRT3*(X(M,1)+X(M,2))
         X(M,2) = SSQRT3*(X(M,1)-X(M,2))
         X(M,1) = XHOLD
  112 CONTINUE
  101  GO TO 200
  103 NP1 = N+1
      NS2 = N/2
      DO 104 K=1,NS2
         KC = NP1-K
         M1 = 0
         DO 114 M=1,LJ,JUMP
         M1 = M1 + 1
         T1 = X(M,K)-X(M,KC)
         T2 = fft_sign%WSAVE(K)*(X(M,K)+X(M,KC))
         XH(M1,K+1) = T1+T2
         XH(M1,KC+1) = T2-T1
  114    CONTINUE
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
      M1 = 0
      DO 123 M=1,LJ,JUMP
         M1 = M1 + 1
         XH(M1,NS2+2) = 4._rk*X(M,NS2+1)
  123 CONTINUE
  124 DO 127 M=1,LOT
         XH(M,1) = 0._rk
  127 CONTINUE
      LNXH = LOT-1 + LOT*(NP1-1) + 1
      LNSV = NP1 + INT(LOG(REAL(NP1))/LOG(2._rk)) + 4
      LNWK = LOT*NP1
!
      !!CALL RFFTMF(LOT,1,NP1,LOT,XH,LNXH,WSAVE(NS2+1),LNSV,WORK,LNWK, fft_sign,IER1)
      CALL RFFTMF(LOT,1,NP1,LOT,XH,LNXH,WORK,LNWK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('MSNTF1',-5)
        GO TO 200
      ENDIF
!
      IF(MOD(NP1,2) .NE. 0) GO TO 30
      DO 20 M=1,LOT
      XH(M,NP1) = XH(M,NP1)+XH(M,NP1)
   20 CONTINUE
   30 SFNP1 = 1._rk/REAL(NP1,kind=rk)
      M1 = 0
      DO 125 M=1,LJ,JUMP
         M1 = M1+1
         X(M,1) = .5_rk*XH(M1,1)
         DSUM(M1) = X(M,1)
  125 CONTINUE
      DO 105 I=3,N,2
         M1 = 0
         DO 115 M=1,LJ,JUMP
            M1 = M1+1
            X(M,I-1) = .5_rk*XH(M1,I)
            DSUM(M1) = DSUM(M1)+.5_rk*XH(M1,I-1)
            X(M,I) = DSUM(M1)
  115    CONTINUE
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 200
      M1 = 0
      DO 116 M=1,LJ,JUMP
         M1 = M1+1
         X(M,N) = .5_rk*XH(M1,N+1)
  116 CONTINUE
  200 CONTINUE
      END subroutine MSNTF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F2KB (IDO,L1,CC,IN1,CH,IN2,WA1)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      REAL(kind=rk)     :: CC(IN1,IDO,2,L1), CH(IN2,IDO,L1,2), WA1(IDO)
!
      DO 101 K=1,L1
         CH(1,1,K,1) = CC(1,1,1,K)+CC(1,IDO,2,K)
         CH(1,1,K,2) = CC(1,1,1,K)-CC(1,IDO,2,K)
 101  CONTINUE
      IF (IDO-2) 107,105,102
 102  IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I

            CH(1,I-1,K,1) = CC(1,I-1,1,K)+CC(1,IC-1,2,K)
            CH(1,I,K,1) = CC(1,I,1,K)-CC(1,IC,2,K)

            CH(1,I-1,K,2) = WA1(I-2)*(CC(1,I-1,1,K)-CC(1,IC-1,2,K))-WA1(I-1)*(CC(1,I,1,K)+CC(1,IC,2,K))
            CH(1,I,K,2) = WA1(I-2)*(CC(1,I,1,K)+CC(1,IC,2,K))+WA1(I-1)*(CC(1,I-1,1,K)-CC(1,IC-1,2,K))

 103     CONTINUE
 104  CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
 105  DO 106 K=1,L1
         CH(1,IDO,K,1) = CC(1,IDO,1,K)+CC(1,IDO,1,K)
         CH(1,IDO,K,2) = -(CC(1,1,2,K)+CC(1,1,2,K))
 106  CONTINUE
 107  continue
      END subroutine R1F2KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F2KF (IDO,L1,CC,IN1,CH,IN2,WA1)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      REAL(kind=rk)     :: CH(IN2,IDO,2,L1),CC(IN1,IDO,L1,2), WA1(IDO)
!
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+CC(1,1,K,2)
         CH(1,IDO,2,K) = CC(1,1,K,1)-CC(1,1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(1,I,1,K) = CC(1,I,K,1)+(WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))
            CH(1,IC,2,K) = (WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-CC(1,I,K,1)
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))
            CH(1,IC-1,2,K) = CC(1,I-1,K,1)-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,1,2,K) = -CC(1,IDO,K,2)
         CH(1,IDO,1,K) = CC(1,IDO,K,1)
  106 CONTINUE
  107 continue
      END subroutine R1F2KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F3KB (IDO,L1,CC,IN1,CH,IN2,WA1,WA2)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TAUR,TAUI
      REAL(kind=rk)     :: CC(IN1,IDO,3,L1),CH(IN2,IDO,L1,3),WA1(IDO),WA2(IDO)
!
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/3._dk
      ARG = 2._dk*PI_long/3._dk
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         CH(1,1,K,1) = CC(1,1,1,K)+2._rk*CC(1,IDO,2,K)
         CH(1,1,K,2) = CC(1,1,1,K)+(2._rk*TAUR)*CC(1,IDO,2,K)-(2._rk*TAUI)*CC(1,1,3,K)
         CH(1,1,K,3) = CC(1,1,1,K)+(2._rk*TAUR)*CC(1,IDO,2,K)+2._rk*TAUI*CC(1,1,3,K)
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
        CH(1,I-1,K,1) = CC(1,I-1,1,K)+(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
        CH(1,I,K,1) = CC(1,I,1,K)+(CC(1,I,3,K)-CC(1,IC,2,K))
        CH(1,I-1,K,2) = WA1(I-2)*&
     & ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))-&
     & (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))&
     &                   -WA1(I-1)*&
     & ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))+&
     & (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))
            CH(1,I,K,2) = WA1(I-2)*&
     & ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))+&
     & (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))&
     &                  +WA1(I-1)*&
     & ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))-&
     & (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))
              CH(1,I-1,K,3) = WA2(I-2)*&
     & ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))+&
     & (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))&
     &   -WA2(I-1)*&
     & ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))-&
     & (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))
            CH(1,I,K,3) = WA2(I-2)*&
     & ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))-&
     & (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))&
     &                 +WA2(I-1)*&
     & ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))+&
     & (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))
  102    CONTINUE
  103 CONTINUE
      END subroutine R1F3KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F3KF (IDO,L1,CC,IN1,CH,IN2,WA1,WA2)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: TAUR,TAUI
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      REAL(kind=rk)     :: CH(IN2,IDO,3,L1),CC(IN1,IDO,L1,3),WA1(IDO),WA2(IDO)
!
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/3._dk
      ARG = 2._dk*PI_long/3._dk
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+(CC(1,1,K,2)+CC(1,1,K,3))
         CH(1,1,3,K) = TAUI*(CC(1,1,K,3)-CC(1,1,K,2))
         CH(1,IDO,2,K) = CC(1,1,K,1)+TAUR*(CC(1,1,K,2)+CC(1,1,K,3))
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+((WA1(I-2)*CC(1,I-1,K,2)+&
     &       WA1(I-1)*CC(1,I,K,2))+(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3)))
            CH(1,I,1,K) = CC(1,I,K,1)+((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3)))
            CH(1,I-1,3,K) = (CC(1,I-1,K,1)+TAUR*((WA1(I-2)*&
     &       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA2(I-2)*&
     &       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))))+(TAUI*((WA1(I-2)*&
     &       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA2(I-2)*&
     &       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))))
            CH(1,IC-1,2,K) = (CC(1,I-1,K,1)+TAUR*((WA1(I-2)*&
     &       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA2(I-2)*&
     &       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))))-(TAUI*((WA1(I-2)*&
     &       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA2(I-2)*&
     &       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))))
            CH(1,I,3,K) = (CC(1,I,K,1)+TAUR*((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2))))
            CH(1,IC,2,K) = (TAUI*((WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2))))-(CC(1,I,K,1)+TAUR*((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))))
  102    CONTINUE
  103 CONTINUE
      END subroutine R1F3KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F4KB (IDO,L1,CC,IN1,CH,IN2,WA1,WA2,WA3)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SQRT2
      REAL(kind=rk)     :: CC(IN1,IDO,4,L1),CH(IN2,IDO,L1,4),WA1(IDO), WA2(IDO), WA3(IDO)
!
      SQRT2=SQRT(2._rk)
      DO 101 K=1,L1
         CH(1,1,K,3) = (CC(1,1,1,K)+CC(1,IDO,4,K))-(CC(1,IDO,2,K)+CC(1,IDO,2,K))
         CH(1,1,K,1) = (CC(1,1,1,K)+CC(1,IDO,4,K))+(CC(1,IDO,2,K)+CC(1,IDO,2,K))
         CH(1,1,K,4) = (CC(1,1,1,K)-CC(1,IDO,4,K))+(CC(1,1,3,K)+CC(1,1,3,K))
         CH(1,1,K,2) = (CC(1,1,1,K)-CC(1,IDO,4,K))-(CC(1,1,3,K)+CC(1,1,3,K))
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
        CH(1,I-1,K,1) = (CC(1,I-1,1,K)+CC(1,IC-1,4,K))+(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
        CH(1,I,K,1) = (CC(1,I,1,K)-CC(1,IC,4,K))+(CC(1,I,3,K)-CC(1,IC,2,K))
        CH(1,I-1,K,2)=WA1(I-2)*((CC(1,I-1,1,K)-CC(1,IC-1,4,K))&
     &  -(CC(1,I,3,K)+CC(1,IC,2,K)))-WA1(I-1)&
     &  *((CC(1,I,1,K)+CC(1,IC,4,K))+(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))
        CH(1,I,K,2)=WA1(I-2)*((CC(1,I,1,K)+CC(1,IC,4,K))&
     &  +(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))+WA1(I-1)&
     &  *((CC(1,I-1,1,K)-CC(1,IC-1,4,K))-(CC(1,I,3,K)+CC(1,IC,2,K)))
        CH(1,I-1,K,3)=WA2(I-2)*((CC(1,I-1,1,K)+CC(1,IC-1,4,K))&
     &  -(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))-WA2(I-1)&
     &  *((CC(1,I,1,K)-CC(1,IC,4,K))-(CC(1,I,3,K)-CC(1,IC,2,K)))
        CH(1,I,K,3)=WA2(I-2)*((CC(1,I,1,K)-CC(1,IC,4,K))&
     &  -(CC(1,I,3,K)-CC(1,IC,2,K)))+WA2(I-1)&
     &  *((CC(1,I-1,1,K)+CC(1,IC-1,4,K))-(CC(1,I-1,3,K)&
     &  +CC(1,IC-1,2,K)))
        CH(1,I-1,K,4)=WA3(I-2)*((CC(1,I-1,1,K)-CC(1,IC-1,4,K))&
     &  +(CC(1,I,3,K)+CC(1,IC,2,K)))-WA3(I-1)&
     & *((CC(1,I,1,K)+CC(1,IC,4,K))-(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))
        CH(1,I,K,4)=WA3(I-2)*((CC(1,I,1,K)+CC(1,IC,4,K))&
     &  -(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))+WA3(I-1)&
     &  *((CC(1,I-1,1,K)-CC(1,IC-1,4,K))+(CC(1,I,3,K)+CC(1,IC,2,K)))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         CH(1,IDO,K,1) = (CC(1,IDO,1,K)+CC(1,IDO,3,K))+(CC(1,IDO,1,K)+CC(1,IDO,3,K))
         CH(1,IDO,K,2) = SQRT2*((CC(1,IDO,1,K)-CC(1,IDO,3,K))-(CC(1,1,2,K)+CC(1,1,4,K)))
         CH(1,IDO,K,3) = (CC(1,1,4,K)-CC(1,1,2,K))+(CC(1,1,4,K)-CC(1,1,2,K))
         CH(1,IDO,K,4) = -SQRT2*((CC(1,IDO,1,K)-CC(1,IDO,3,K))+(CC(1,1,2,K)+CC(1,1,4,K)))
  106 CONTINUE
  107 continue
      END subroutine R1F4KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F4KF (IDO,L1,CC,IN1,CH,IN2,WA1,WA2,WA3)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: HSQT2
      REAL(kind=rk)     :: CC(IN1,IDO,L1,4),CH(IN2,IDO,4,L1),WA1(IDO),WA2(IDO),WA3(IDO)
!
      HSQT2=SQRT(2._rk)/2._rk
      DO 101 K=1,L1
         CH(1,1,1,K) = (CC(1,1,K,2)+CC(1,1,K,4))+(CC(1,1,K,1)+CC(1,1,K,3))
         CH(1,IDO,4,K) = (CC(1,1,K,1)+CC(1,1,K,3))-(CC(1,1,K,2)+CC(1,1,K,4))
         CH(1,IDO,2,K) = CC(1,1,K,1)-CC(1,1,K,3)
         CH(1,1,3,K) = CC(1,1,K,4)-CC(1,1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(1,I-1,1,K) = ((WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2))+(WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4)))+(CC(1,I-1,K,1)+(WA2(I-2)*CC(1,I-1,K,3)+&
     &       WA2(I-1)*CC(1,I,K,3)))
            CH(1,IC-1,4,K) = (CC(1,I-1,K,1)+(WA2(I-2)*CC(1,I-1,K,3)+&
     &       WA2(I-1)*CC(1,I,K,3)))-((WA1(I-2)*CC(1,I-1,K,2)+&
     &       WA1(I-1)*CC(1,I,K,2))+(WA3(I-2)*CC(1,I-1,K,4)+&
     &       WA3(I-1)*CC(1,I,K,4)))
            CH(1,I,1,K) = ((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*&
     &       CC(1,I-1,K,2))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4)))+(CC(1,I,K,1)+(WA2(I-2)*CC(1,I,K,3)-&
     &       WA2(I-1)*CC(1,I-1,K,3)))
            CH(1,IC,4,K) = ((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*&
     &       CC(1,I-1,K,2))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4)))-(CC(1,I,K,1)+(WA2(I-2)*CC(1,I,K,3)-&
     &       WA2(I-1)*CC(1,I-1,K,3)))
            CH(1,I-1,3,K) = ((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*&
     &       CC(1,I-1,K,2))-(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4)))+(CC(1,I-1,K,1)-(WA2(I-2)*CC(1,I-1,K,3)+&
     &       WA2(I-1)*CC(1,I,K,3)))
            CH(1,IC-1,2,K) = (CC(1,I-1,K,1)-(WA2(I-2)*CC(1,I-1,K,3)+&
     &       WA2(I-1)*CC(1,I,K,3)))-((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*&
     &       CC(1,I-1,K,2))-(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4)))
            CH(1,I,3,K) = ((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2)))+(CC(1,I,K,1)-(WA2(I-2)*CC(1,I,K,3)-&
     &       WA2(I-1)*CC(1,I-1,K,3)))
            CH(1,IC,2,K) = ((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2)))-(CC(1,I,K,1)-(WA2(I-2)*CC(1,I,K,3)-&
     &       WA2(I-1)*CC(1,I-1,K,3)))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
            CH(1,IDO,1,K) = (HSQT2*(CC(1,IDO,K,2)-CC(1,IDO,K,4)))+CC(1,IDO,K,1)
            CH(1,IDO,3,K) = CC(1,IDO,K,1)-(HSQT2*(CC(1,IDO,K,2)-CC(1,IDO,K,4)))
            CH(1,1,2,K) = (-HSQT2*(CC(1,IDO,K,2)+CC(1,IDO,K,4)))-CC(1,IDO,K,3)
            CH(1,1,4,K) = (-HSQT2*(CC(1,IDO,K,2)+CC(1,IDO,K,4)))+CC(1,IDO,K,3)
  106 CONTINUE
  107 continue
      END subroutine R1F4KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F5KB (IDO,L1,CC,IN1,CH,IN2,WA1,WA2,WA3,WA4)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      REAL(kind=rk)     :: CC(IN1,IDO,5,L1),CH(IN2,IDO,L1,5),WA1(IDO),WA2(IDO),WA3(IDO),WA4(IDO)
!
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/5._dk
      ARG = 2._dk*PI_long/5._dk
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2._rk*ARG)
      TI12=SIN(2._rk*ARG)
      DO 101 K=1,L1
         CH(1,1,K,1) = CC(1,1,1,K)+2._rk*CC(1,IDO,2,K)+2._rk*CC(1,IDO,4,K)
         CH(1,1,K,2) = (CC(1,1,1,K)+TR11*2._rk*CC(1,IDO,2,K)&
     &   +TR12*2._rk*CC(1,IDO,4,K))-(TI11*2._rk*CC(1,1,3,K)&
     &   +TI12*2._rk*CC(1,1,5,K))
         CH(1,1,K,3) = (CC(1,1,1,K)+TR12*2._rk*CC(1,IDO,2,K)&
     &   +TR11*2._rk*CC(1,IDO,4,K))-(TI12*2._rk*CC(1,1,3,K)&
     &   -TI11*2._rk*CC(1,1,5,K))
         CH(1,1,K,4) = (CC(1,1,1,K)+TR12*2._rk*CC(1,IDO,2,K)&
     &   +TR11*2._rk*CC(1,IDO,4,K))+(TI12*2._rk*CC(1,1,3,K)&
     &   -TI11*2._rk*CC(1,1,5,K))
         CH(1,1,K,5) = (CC(1,1,1,K)+TR11*2._rk*CC(1,IDO,2,K)&
     &   +TR12*2._rk*CC(1,IDO,4,K))+(TI11*2._rk*CC(1,1,3,K)&
     &   +TI12*2._rk*CC(1,1,5,K))
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
        CH(1,I-1,K,1) = CC(1,I-1,1,K)+(CC(1,I-1,3,K)+CC(1,IC-1,2,K))+(CC(1,I-1,5,K)+CC(1,IC-1,4,K))
        CH(1,I,K,1) = CC(1,I,1,K)+(CC(1,I,3,K)-CC(1,IC,2,K))+(CC(1,I,5,K)-CC(1,IC,4,K))
        CH(1,I-1,K,2) = WA1(I-2)*((CC(1,I-1,1,K)+TR11*&
     &  (CC(1,I-1,3,K)+CC(1,IC-1,2,K))+TR12&
     &  *(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))-(TI11*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))+TI12*(CC(1,I,5,K)+CC(1,IC,4,K))))&
     &  -WA1(I-1)*((CC(1,I,1,K)+TR11*(CC(1,I,3,K)-CC(1,IC,2,K))&
     &  +TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))+(TI11*(CC(1,I-1,3,K)&
     &  -CC(1,IC-1,2,K))+TI12*(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,2) = WA1(I-2)*((CC(1,I,1,K)+TR11*(CC(1,I,3,K)&
     &  -CC(1,IC,2,K))+TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))&
     &  +(TI11*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))+TI12&
     &  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))+WA1(I-1)&
     &  *((CC(1,I-1,1,K)+TR11*(CC(1,I-1,3,K)&
     &  +CC(1,IC-1,2,K))+TR12*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))&
     &  -(TI11*(CC(1,I,3,K)+CC(1,IC,2,K))+TI12&
     &  *(CC(1,I,5,K)+CC(1,IC,4,K))))
        CH(1,I-1,K,3) = WA2(I-2)&
     &  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))&
     &  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))-(TI12*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))&
     & -WA2(I-1)&
     & *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-&
     &  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))&
     &  +(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11&
     &  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,3) = WA2(I-2)&
     & *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-&
     &  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))&
     &  +(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11&
     &  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))&
     &  +WA2(I-1)&
     &  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))&
     &  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))-(TI12*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))
        CH(1,I-1,K,4) = WA3(I-2)&
     &  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))&
     &  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI12*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))&
     &  -WA3(I-1)&
     & *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-&
     &  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))&
     &  -(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11&
     &  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,4) = WA3(I-2)&
     & *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-&
     &  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))&
     &  -(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11&
     &  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))&
     &  +WA3(I-1)&
     &  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))&
     &  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI12*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))
        CH(1,I-1,K,5) = WA4(I-2)&
     &  *((CC(1,I-1,1,K)+TR11*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))&
     &  +TR12*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI11*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))+TI12*(CC(1,I,5,K)+CC(1,IC,4,K))))&
     &  -WA4(I-1)&
     &  *((CC(1,I,1,K)+TR11*(CC(1,I,3,K)-CC(1,IC,2,K))&
     &  +TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))-(TI11*(CC(1,I-1,3,K)&
     &  -CC(1,IC-1,2,K))+TI12*(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,5) = WA4(I-2)&
     &  *((CC(1,I,1,K)+TR11*(CC(1,I,3,K)-CC(1,IC,2,K))&
     &  +TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))-(TI11*(CC(1,I-1,3,K)&
     &  -CC(1,IC-1,2,K))+TI12*(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))&
     &  +WA4(I-1)&
     &  *((CC(1,I-1,1,K)+TR11*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))&
     &  +TR12*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI11*(CC(1,I,3,K)&
     &  +CC(1,IC,2,K))+TI12*(CC(1,I,5,K)+CC(1,IC,4,K))))
  102    CONTINUE
  103 CONTINUE
      END subroutine R1F5KB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1F5KF (IDO,L1,CC,IN1,CH,IN2,WA1,WA2,WA3,WA4)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      real(kind=rk)     :: TR11,TI11,TR12,TI12
      REAL(kind=rk)     :: CC(IN1,IDO,L1,5),CH(IN2,IDO,5,L1),WA1(IDO),WA2(IDO),WA3(IDO),WA4(IDO)
!
      !!ARG=2._dk*4._dk*ATAN(1.0_dk)/5._dk
      ARG = 2._dk*PI_long/5._dk
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2._rk*ARG)
      TI12=SIN(2._rk*ARG)
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+(CC(1,1,K,5)+CC(1,1,K,2))+(CC(1,1,K,4)+CC(1,1,K,3))
         CH(1,IDO,2,K) = CC(1,1,K,1)+TR11*(CC(1,1,K,5)+CC(1,1,K,2))+TR12*(CC(1,1,K,4)+CC(1,1,K,3))
         CH(1,1,3,K) = TI11*(CC(1,1,K,5)-CC(1,1,K,2))+TI12*(CC(1,1,K,4)-CC(1,1,K,3))
         CH(1,IDO,4,K) = CC(1,1,K,1)+TR12*(CC(1,1,K,5)+CC(1,1,K,2))+TR11*(CC(1,1,K,4)+CC(1,1,K,3))
         CH(1,1,5,K) = TI12*(CC(1,1,K,5)-CC(1,1,K,2))-TI11*(CC(1,1,K,4)-CC(1,1,K,3))
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+((WA1(I-2)*CC(1,I-1,K,2)+&
     &       WA1(I-1)*CC(1,I,K,2))+(WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*&
     &       CC(1,I,K,5)))+((WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))+(WA3(I-2)*CC(1,I-1,K,4)+&
     &       WA3(I-1)*CC(1,I,K,4)))
            CH(1,I,1,K) = CC(1,I,K,1)+((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*&
     &       CC(1,I-1,K,5)))+((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4)))
            CH(1,I-1,3,K) = CC(1,I-1,K,1)+TR11*&
     &      ( WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2)&
     &       +WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5))+TR12*&
     &      ( WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3)&
     &       +WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))+TI11*&
     &      ( WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2)&
     &       -(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))+TI12*&
     &      ( WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3)&
     &       -(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4)))
            CH(1,IC-1,2,K) = CC(1,I-1,K,1)+TR11*&
     &      ( WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2)&
     &       +WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5))+TR12*&
     &     ( WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3)&
     &      +WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))-(TI11*&
     &      ( WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2)&
     &       -(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))+TI12*&
     &      ( WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3)&
     &       -(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4))))
            CH(1,I,3,K) = (CC(1,I,K,1)+TR11*((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*&
     &       CC(1,I-1,K,5)))+TR12*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4))))+(TI11*((WA4(I-2)*CC(1,I-1,K,5)+&
     &       WA4(I-1)*CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2)))+TI12*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))))
            CH(1,IC,2,K) = (TI11*((WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*&
     &       CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2)))+TI12*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))))-(CC(1,I,K,1)+TR11*((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*&
     &       CC(1,I-1,K,5)))+TR12*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4))))
            CH(1,I-1,5,K) = (CC(1,I-1,K,1)+TR12*((WA1(I-2)*&
     &       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA4(I-2)*&
     &       CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5)))+TR11*((WA2(I-2)*&
     &       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))+(WA3(I-2)*&
     &       CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))))+(TI12*((WA1(I-2)*&
     &       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA4(I-2)*&
     &       CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))-TI11*((WA2(I-2)*&
     &       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))-(WA3(I-2)*&
     &       CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4))))
            CH(1,IC-1,4,K) = (CC(1,I-1,K,1)+TR12*((WA1(I-2)*&
     &       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA4(I-2)*&
     &       CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5)))+TR11*((WA2(I-2)*&
     &       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))+(WA3(I-2)*&
     &       CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))))-(TI12*((WA1(I-2)*&
     &       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA4(I-2)*&
     &       CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))-TI11*((WA2(I-2)*&
     &       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))-(WA3(I-2)*&
     &       CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4))))
            CH(1,I,5,K) = (CC(1,I,K,1)+TR12*((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*&
     &       CC(1,I-1,K,5)))+TR11*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4))))+(TI12*((WA4(I-2)*CC(1,I-1,K,5)+&
     &       WA4(I-1)*CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2)))-TI11*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))))
            CH(1,IC,4,K) = (TI12*((WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*&
     &       CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*&
     &       CC(1,I,K,2)))-TI11*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*&
     &       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*&
     &       CC(1,I,K,3))))-(CC(1,I,K,1)+TR12*((WA1(I-2)*CC(1,I,K,2)-&
     &       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*&
     &       CC(1,I-1,K,5)))+TR11*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*&
     &       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*&
     &       CC(1,I-1,K,4))))
  102    CONTINUE
  103 CONTINUE
      END subroutine R1F5KF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1FGKB (IDO,IP,L1,IDL1,CC,C1,C2,IN1,CH,CH2,IN2,WA)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: NBD,IK
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: AR1,AI1,AR2,AI2
      real(kind=rk)     :: DC2,DCP,DS2,DSP
      real(kind=rk)     :: TPI           
      real(kind=rk)     :: AR,ARG,AR1H,AR2H
      REAL(kind=rk)     :: CH(IN2,IDO,L1,IP),CC(IN1,IDO,IP,L1),C1(IN1,IDO,L1,IP),C2(IN1,IDL1,IP),CH2(IN2,IDL1,IP)
      !!REAL(kind=rk)     :: WA(IDO) !Error!
      REAL(kind=rk)     :: WA(*)
!
      !!TPI=2._dk*4._dk*ATAN(1.0_dk)
      TPI=2._dk*PI_long
      ARG = TPI/REAL(IP,kind=rk)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(1,I,K,1) = CC(1,I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(1,I,K,1) = CC(1,I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,1,K,J) = CC(1,IDO,J2-2,K)+CC(1,IDO,J2-2,K)
            CH(1,1,K,JC) = CC(1,1,J2-1,K)+CC(1,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(1,I-1,K,J) = CC(1,I-1,2*J-1,K)+CC(1,IC-1,2*J-2,K)
               CH(1,I-1,K,JC) = CC(1,I-1,2*J-1,K)-CC(1,IC-1,2*J-2,K)
               CH(1,I,K,J) = CC(1,I,2*J-1,K)-CC(1,IC,2*J-2,K)
               CH(1,I,K,JC) = CC(1,I,2*J-1,K)+CC(1,IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(1,I-1,K,J) = CC(1,I-1,2*J-1,K)+CC(1,IC-1,2*J-2,K)
               CH(1,I-1,K,JC) = CC(1,I-1,2*J-1,K)-CC(1,IC-1,2*J-2,K)
               CH(1,I,K,J) = CC(1,I,2*J-1,K)-CC(1,IC,2*J-2,K)
               CH(1,I,K,JC) = CC(1,I,2*J-1,K)+CC(1,IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1._rk
      AI1 = 0._rk
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(1,IK,L) = CH2(1,IK,1)+AR1*CH2(1,IK,2)
            C2(1,IK,LC) = AI1*CH2(1,IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(1,IK,L) = C2(1,IK,L)+AR2*CH2(1,IK,J)
               C2(1,IK,LC) = C2(1,IK,LC)+AI2*CH2(1,IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(1,IK,1) = CH2(1,IK,1)+CH2(1,IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,1,K,J) = C1(1,1,K,J)-C1(1,1,K,JC)
            CH(1,1,K,JC) = C1(1,1,K,J)+C1(1,1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(1,I-1,K,J) = C1(1,I-1,K,J)-C1(1,I,K,JC)
               CH(1,I-1,K,JC) = C1(1,I-1,K,J)+C1(1,I,K,JC)
               CH(1,I,K,J) = C1(1,I,K,J)+C1(1,I-1,K,JC)
               CH(1,I,K,JC) = C1(1,I,K,J)-C1(1,I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(1,I-1,K,J) = C1(1,I-1,K,J)-C1(1,I,K,JC)
               CH(1,I-1,K,JC) = C1(1,I-1,K,J)+C1(1,I,K,JC)
               CH(1,I,K,J) = C1(1,I,K,J)+C1(1,I-1,K,JC)
               CH(1,I,K,JC) = C1(1,I,K,J)-C1(1,I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(1,IK,1) = CH2(1,IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,1,K,J) = CH(1,1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(1,I-1,K,J) = WA(IDIJ-1)*CH(1,I-1,K,J)-WA(IDIJ)*CH(1,I,K,J)
               C1(1,I,K,J) = WA(IDIJ-1)*CH(1,I,K,J)+WA(IDIJ)*CH(1,I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(1,I-1,K,J) = WA(IDIJ-1)*CH(1,I-1,K,J)-WA(IDIJ)*CH(1,I,K,J)
               C1(1,I,K,J) = WA(IDIJ-1)*CH(1,I,K,J)+WA(IDIJ)*CH(1,I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 continue
      END subroutine R1FGKB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE R1FGKF (IDO,IP,L1,IDL1,CC,C1,C2,IN1,CH,CH2,IN2,WA)
      integer           :: N,M,IDL1,IM1,IM2,IN,IN1,IN2
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: NBD,IK
      real(kind=rk)     :: AR1,AI1,AR2,AI2
      real(kind=rk)     :: DC2,DCP,DS2,DSP
      real(kind=rk)     :: ARG,AR1H,AR2H
      real(kind=rk)     :: TPI           
      REAL(kind=rk)     :: CH(IN2,IDO,L1,IP),CC(IN1,IDO,IP,L1),C1(IN1,IDO,L1,IP),C2(IN1,IDL1,IP),CH2(IN2,IDL1,IP)
      REAL(kind=rk)     :: WA(:)  ! falsch?
     !! REAL(kind=rk)     :: WA(IDO)  ! wrong dimension!
!
      !!TPI=2._dk*4._dk*ATAN(1.0_dk)
      TPI=2._dk*PI_long
      ARG = TPI/REAL(IP,kind=rk)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(1,IK,1) = C2(1,IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,1,K,J) = C1(1,1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
 !!     IS = -IDO
 !!     DO J=2,IP
 !!        IS = IS+IDO
 !!        IDIJ = IS
 !!        DO I=3,IDO,2
 !!           IDIJ = IDIJ+2
 !!           write(*,*) 'J,IS,I,IDIJ,WA(IDIJ-1),WA(IDIJ)=',J,IS,I,IDIJ,WA(IDIJ-1),WA(IDIJ)
 !!        enddo   
 !!     enddo   

      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(1,I-1,K,J) = WA(IDIJ-1)*C1(1,I-1,K,J)+WA(IDIJ)*C1(1,I,K,J)
               CH(1,I,K,J) = WA(IDIJ-1)*C1(1,I,K,J)-WA(IDIJ)*C1(1,I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(1,I-1,K,J) = WA(IDIJ-1)*C1(1,I-1,K,J)+WA(IDIJ)*C1(1,I,K,J)
               CH(1,I,K,J) = WA(IDIJ-1)*C1(1,I,K,J)-WA(IDIJ)*C1(1,I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               C1(1,I-1,K,J) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               C1(1,I-1,K,JC) = CH(1,I,K,J)-CH(1,I,K,JC)
               C1(1,I,K,J) = CH(1,I,K,J)+CH(1,I,K,JC)
               C1(1,I,K,JC) = CH(1,I-1,K,JC)-CH(1,I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(1,I-1,K,J) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               C1(1,I-1,K,JC) = CH(1,I,K,J)-CH(1,I,K,JC)
               C1(1,I,K,J) = CH(1,I,K,J)+CH(1,I,K,JC)
               C1(1,I,K,JC) = CH(1,I-1,K,JC)-CH(1,I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(1,IK,1) = CH2(1,IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,1,K,J) = CH(1,1,K,J)+CH(1,1,K,JC)
            C1(1,1,K,JC) = CH(1,1,K,JC)-CH(1,1,K,J)
  122    CONTINUE
  123 CONTINUE
!
      AR1 = 1._rk
      AI1 = 0._rk
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(1,IK,L) = C2(1,IK,1)+AR1*C2(1,IK,2)
            CH2(1,IK,LC) = AI1*C2(1,IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(1,IK,L) = CH2(1,IK,L)+AR2*C2(1,IK,J)
               CH2(1,IK,LC) = CH2(1,IK,LC)+AI2*C2(1,IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(1,IK,1) = CH2(1,IK,1)+C2(1,IK,J)
  128    CONTINUE
  129 CONTINUE
!
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(1,I,1,K) = CH(1,I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(1,I,1,K) = CH(1,I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(1,IDO,J2-2,K) = CH(1,1,K,J)
            CC(1,1,J2-1,K) = CH(1,1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(1,I-1,J2-1,K) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               CC(1,IC-1,J2-2,K) = CH(1,I-1,K,J)-CH(1,I-1,K,JC)
               CC(1,I,J2-1,K) = CH(1,I,K,J)+CH(1,I,K,JC)
               CC(1,IC,J2-2,K) = CH(1,I,K,JC)-CH(1,I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(1,I-1,J2-1,K) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               CC(1,IC-1,J2-2,K) = CH(1,I-1,K,J)-CH(1,I-1,K,JC)
               CC(1,I,J2-1,K) = CH(1,I,K,J)+CH(1,I,K,JC)
               CC(1,IC,J2-2,K) = CH(1,I,K,JC)-CH(1,I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      END subroutine R1FGKF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      subroutine r2w(ldr,ldw,l,m,r,w)
      integer       :: ldr,ldw,l,m,i,j
      real(kind=rk) :: r(ldr,*),w(ldw,*)
      do j=1,m
      do i=1,l
      w(i,j) = r( i,j)
      end do
      end do
      return
      end subroutine r2w
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFT1B ( N, INC, R, LENR, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE RFFT1B ( N, INC, R, LENR, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE RFFT1B ( N, INC, R, LENR, fft_sign, IER)
      integer     :: N, INC, LENR
      !!integer    ::  LENWRK
      integer    ::  IER
      !!integer   :: LENSAV
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: R(LENR)
      type(fft_real_sign_type) :: fft_sign
!
      if( debug > 0 ) write(*,*) 'RFFT1B start'
      IER = 0
!
      IF (LENR .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('RFFT1B ', 6)
 !!     ELSEIF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
 !!       IER = 2
 !!       CALL XERFFT ('RFFT1B ', 8)
   !!   ELSEIF (LENWRK .LT. N) THEN
   !!     IER = 3
   !!     CALL XERFFT ('RFFT1B ', 10)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!CALL RFFTB1 (N,INC,R,WORK,WSAVE,WSAVE(N+1))
      CALL RFFTB1 (N,INC,R,fft_sign%WORK,fft_sign%WA,fft_sign%FAC)
      if( debug > 0 ) write(*,*) 'RFFT1B end'
      END subroutine RFFT1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFT1F (N,INC,R,LENR,WSAVE,LENSAV,WORK,LENWRK,fft_sign,IER)
      !!SUBROUTINE RFFT1F (N,INC,R,LENR,WORK,LENWRK,fft_sign,IER)
      SUBROUTINE RFFT1F (N,INC,R,LENR,fft_sign,IER)
      integer     :: N, INC, LENR
      !!integer    ::  LENWRK
      integer    ::  IER
      !!integer   :: LENSAV
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: R(LENR)
      type(fft_real_sign_type) :: fft_sign
!
      if( debug > 0 ) write(*,*) 'RFFT1F_no_if( debug > 0 ) work_alloc start'
      IER = 0
!
      IF (LENR .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('RFFT1F ', 6)
    !!  ELSEIF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
    !!    IER = 2
    !!    CALL XERFFT ('RFFT1F ', 8)
    !!    goto 100
     !! ELSEIF (LENWRK .LT. N) THEN
     !!   IER = 3
     !!   CALL XERFFT ('RFFT1F ', 10)
     !!   goto 100
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!CALL RFFTF1 (N,INC,R,WORK,WSAVE,WSAVE(N+1))
      CALL RFFTF1 (N,INC,R,fft_sign%WORK,fft_sign%WA,fft_sign%FAC)
      100 continue
      if( debug > 0 ) write(*,*) 'RFFT1F_no_if( debug > 0 ) work_alloc end'
      END subroutine RFFT1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFT1I ( N, WSAVE, LENSAV, fft_sign, IER )
      SUBROUTINE RFFT1I ( N, fft_sign, IER )
      integer     :: N
      integer    ::  IER
      !!integer   :: LENSAV
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      integer                  :: ibase2
      integer                  :: NF
      integer                  :: k
      integer                  :: size_WA
      integer     :: LENWRK
      type(fft_real_sign_type) :: fft_sign
!
      if( debug > 0 ) write(*,*) 'RFFT1I start'

      IER = 0
!
     !! IF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
     !!   IER = 2
     !!   CALL XERFFT ('RFFT1I ', 3)
     !! ENDIF
!
  !!    if( fft_sign%switch == 2 .and. fft_sign%N == N .and. allocated(fft_sign%FAC) .and. allocated(fft_sign%WA)) then
  !!         goto 1111
  !!    endif

      IF (N .EQ. 1) RETURN
!
      call all_FFT1I ( N, fft_sign)

      fft_sign%sign_type=RFFT1_sign_type

      LENWRK=N
      if(allocated(fft_sign%work) ) deallocate(fft_sign%work)
      allocate(fft_sign%work(LENWRK))

      1111 continue
      if( debug > 0 ) write(*,*) 'RFFT1I end'
      END subroutine RFFT1I
!*******************************************************************************************************************
      SUBROUTINE all_FFT1I ( N, fft_sign)
      integer                  :: N
      integer                  ::  IER
      integer                  :: ibase2
      integer                  :: NF
      integer                  :: k
      integer                  :: size_WA
      type(fft_real_sign_type) :: fft_sign
!
      if( debug > 0 ) write(*,*) 'all_FFT1I start'

      IER = 0
!
      fft_sign%switch = 2

      fft_sign%N = N
!
      call factor_RFFTI1 (N,NF,fft_sign%FAC)

      call size_WA_RFFTI1 (N,NF,fft_sign%FAC,size_WA)

       !!write(*,'(a,i0)') 'all_FFT1I: size_WA=',size_WA

      if(allocated(fft_sign%WA) ) deallocate(fft_sign%WA)
      allocate(fft_sign%WA(size_WA))

      call WA_RFFTI1 (N,NF,fft_sign%WA,fft_sign%FAC)

      if( debug > 0 ) write(*,*) 'all_FFT1I end'
      END subroutine all_FFT1I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE RFFT1I_orig ( N, WSAVE, LENSAV, IER )
      integer     :: N
      integer    ::  IER
      integer   :: LENSAV
      REAL(kind=rk)     :: WSAVE(LENSAV)
!
      IER = 0
!
      IF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
        IER = 2
        CALL XERFFT ('RFFT1I ', 3)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      CALL RFFTI1 (N,WSAVE(1),WSAVE(N+1))
      END subroutine RFFT1I_orig
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFT2B (LDIM, L, M, R, WSAVE, LENSAV, WORK,LENWRK, fft_sign, IER)
      SUBROUTINE RFFT2B (LDIM, L, M, R, WORK,LENWRK, fft_sign, IER)
      integer     :: LDIM, L, M
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: I,J,K               
      integer           :: LDH,LDW,LDX         
      integer           :: LWSAV,MWSAV,MMSAV         
      integer           :: MODL,MODM                 
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: R(LDIM,M)
      type(fft_real_sign_type),allocatable,dimension(:),intent(in) :: fft_sign
!
!
! INITIALIZE IER
!
      IER = 0
!
! VERIFY LENSAV
!
      LWSAV =   L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4
      MWSAV =   2*M+INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
      MMSAV =   M+INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
      MODL = MOD(L,2)
      MODM = MOD(M,2)
!
!!      IF (LENSAV .LT. LWSAV+MWSAV+MMSAV) THEN
!!        IER = 2
!!        CALL XERFFT ('RFFT2B', 6)
!!        GO TO 100
!!      ENDIF
!
! VERIFY LENWRK
!
      IF (LENWRK .LT. (L+1)*M) THEN
        IER = 3
        CALL XERFFT ('RFFT2B', 8)
        GO TO 100
      ENDIF
!
! VERIFY LDIM IS AS BIG AS L
!
      IF (LDIM .LT. L) THEN
        IER = 5
        CALL XERFFT ('RFFT2B', -6)
        GO TO 100
      ENDIF
!
! TRANSFORM SECOND DIMENSION OF ARRAY
!
      DO J=2,2*((M+1)/2)-1
      R(1,J) = R(1,J)+R(1,J)
      END DO
      DO J=3,M,2
      R(1,J) = -R(1,J)
      END DO
      !!CALL RFFTMB(1,1,M,LDIM,R,M*LDIM,WSAVE(LWSAV+MWSAV+1),MMSAV,WORK,LENWRK,fft_sign(3),IER1)
      CALL RFFTMB(1,1,M,LDIM,R,M*LDIM,WORK,LENWRK,fft_sign(3),IER1)
      LDH = INT((L+1)/2)
      IF(LDH .GT. 1) THEN
      LDW = LDH+LDH
!
!     R AND WORK ARE SWITCHED BECAUSE THE THE FIRST DIMENSION
!     OF THE INPUT TO COMPLEX CFFTMF MUST BE EVEN.
!
      CALL R2W(LDIM,LDW,L,M,R,WORK)
      !!CALL CFFTMB(LDH-1,1,M,LDH,WORK(2),LDH*M,WSAVE(LWSAV+1),MWSAV,R,L*M, IER1)
      CALL CFFTMB(LDH-1,1,M,LDH,WORK(2),LDH*M,R,L*M, fft_sign(2), IER1)
      IF(IER1 .NE. 0) THEN
         IER=20
         CALL XERFFT('RFFT2B',-5)
         GO TO 100
      END IF
      CALL W2R(LDIM,LDW,L,M,R,WORK)
      END IF
!
      IF(MODL .EQ. 0) THEN
      DO J=2,2*((M+1)/2)-1
      R(L,J) = R(L,J)+R(L,J)
      END DO
      DO J=3,M,2
      R(L,J) = -R(L,J)
      END DO
      !!CALL RFFTMB(1,1,M,LDIM,R(L,1),M*LDIM,WSAVE(LWSAV+MWSAV+1),MMSAV,WORK,LENWRK,fft_sign(3),IER1)
      CALL RFFTMB(1,1,M,LDIM,R(L,1),M*LDIM,WORK,LENWRK,fft_sign(3),IER1)
      END IF
!
!     PRINT*, 'BACKWARD TRANSFORM IN THE J DIRECTION'
!     DO I=1,L
!       PRINT*, (R(I,J),J=1,M)
!     END DO
!
! TRANSFORM FIRST DIMENSION OF ARRAY
!
      LDX = 2*INT((L+1)/2)-1
      DO I=2,LDX
      DO J=1,M
      R(I,J) = R(I,J)+R(I,J)
      END DO
      END DO
      DO J=1,M
      DO I=3,LDX,2
      R(I,J) = -R(I,J)
      END DO
      END DO
      !!CALL RFFTMB(M,LDIM,L,1,R,M*LDIM,WSAVE(1),L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4,WORK,LENWRK,fft_sign(1),IER1)
      CALL RFFTMB(M,LDIM,L,1,R,M*LDIM,WORK,LENWRK,fft_sign(1),IER1)
!

!
!     PRINT*, 'BACKWARD TRANSFORM IN THE I DIRECTION'
!     DO I=1,L
!       PRINT*, (R(I,J),J=1,M)
!     END DO
!
      IF(IER1 .NE. 0 ) THEN
         IER=20
         CALL XERFFT('RFFT2B',-5)
         GO TO 100
      ENDIF
!
      IF(IER1 .NE. 0 ) THEN
         IER=20
         CALL XERFFT('RFFT2B',-5)
         GO TO 100
      ENDIF
!
  100 CONTINUE
!
      END subroutine RFFT2B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFT2F (LDIM, L, M, R, WSAVE, LENSAV, WORK,LENWRK, fft_sign, IER)
      SUBROUTINE RFFT2F (LDIM, L, M, R, WORK,LENWRK, fft_sign, IER)
      integer     :: LDIM, L, M
      integer    ::  LENWRK, IER, IDX, MODL, MODM,IDH, IDW
      !!integer   :: LENSAV
      integer           :: IER1
      integer           :: I,J,K               
      integer           :: LDH,LDW,LDX         
      integer           :: LWSAV,MWSAV,MMSAV         
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: R(LDIM,M)
      type(fft_real_sign_type),allocatable,dimension(:),intent(in) :: fft_sign
!
!
! INITIALIZE IER
!
      IER = 0
!
! VERIFY LENSAV
!
      LWSAV =   L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4
      MWSAV =   2*M+INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
      MMSAV =   M+INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
!
!!      IF (LENSAV .LT. LWSAV+MWSAV+MMSAV) THEN
!!        IER = 2
!!        CALL XERFFT ('RFFT2F', 6)
!!        GO TO 100
!!      ENDIF
!
! VERIFY LENWRK
!
      IF (LENWRK .LT. (L+1)*M) THEN
        IER = 3
        CALL XERFFT ('RFFT2F', 8)
        GO TO 100
      ENDIF
!
! VERIFY LDIM IS AS BIG AS L
!
      IF (LDIM .LT. L) THEN
        IER = 5
        CALL XERFFT ('RFFT2F', -6)
        GO TO 100
      ENDIF
!
! TRANSFORM FIRST DIMENSION OF ARRAY
!
      !!CALL RFFTMF(M,LDIM,L,1,R,M*LDIM,WSAVE(1),L+INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4,WORK,LENWRK,fft_sign(1),IER1)
      CALL RFFTMF(M,LDIM,L,1,R,M*LDIM,WORK,LENWRK,fft_sign(1),IER1)
!
      IF(IER1 .NE. 0 ) THEN
         IER=20
         CALL XERFFT('RFFT2F',-5)
         GO TO 100
      ENDIF
!
      LDX = 2*INT((L+1)/2)-1
      DO I=2,LDX
      DO J=1,M
      R(I,J) = .5_rk*R(I,J)
      END DO
      END DO
      DO J=1,M
      DO I=3,LDX,2
      R(I,J) = -R(I,J)
      END DO
      END DO
!
!     PRINT*, 'FORWARD TRANSFORM IN THE I DIRECTION'
!     DO I=1,L
!       PRINT*, (R(I,J),J=1,M)
!     END DO
!
! RESHUFFLE TO ADD IN NYQUIST IMAGINARY COMPONENTS
!
      MODL = MOD(L,2)
      MODM = MOD(M,2)
!
! TRANSFORM SECOND DIMENSION OF ARRAY
!
      !!CALL RFFTMF(1,1,M,LDIM,R,M*LDIM,WSAVE(LWSAV+MWSAV+1),MMSAV,WORK,LENWRK,fft_sign(3),IER1)
      CALL RFFTMF(1,1,M,LDIM,R,M*LDIM,WORK,LENWRK,fft_sign(3),IER1)
      DO J=2,2*((M+1)/2)-1
      R(1,J) = .5_rk*R(1,J)
      END DO
      DO J=3,M,2
      R(1,J) = -R(1,J)
      END DO
      LDH = INT((L+1)/2)
      IF(LDH .GT. 1) THEN
      LDW = LDH+LDH
!
!     R AND WORK ARE SWITCHED BECAUSE THE THE FIRST DIMENSION
!     OF THE INPUT TO COMPLEX CFFTMF MUST BE EVEN.
!
      CALL R2W(LDIM,LDW,L,M,R,WORK)
      !!CALL CFFTMF(LDH-1,1,M,LDH,WORK(2),LDH*M,WSAVE(LWSAV+1),MWSAV,R,L*M, IER1)
      CALL CFFTMF(LDH-1,1,M,LDH,WORK(2),LDH*M,R,L*M, fft_sign(2), IER1)
      IF(IER1 .NE. 0 ) THEN
         IER=20
         CALL XERFFT('RFFT2F',-5)
         GO TO 100
      ENDIF
      CALL W2R(LDIM,LDW,L,M,R,WORK)
      END IF
!
      IF(MODL .EQ. 0) THEN
      !!CALL RFFTMF(1,1,M,LDIM,R(L,1),M*LDIM,WSAVE(LWSAV+MWSAV+1),MMSAV,WORK,LENWRK,fft_sign(3),IER1)
      CALL RFFTMF(1,1,M,LDIM,R(L,1),M*LDIM,WORK,LENWRK,fft_sign(3),IER1)
      DO J=2,2*((M+1)/2)-1
      R(L,J) = .5_rk*R(L,J)
      END DO
      DO J=3,M,2
      R(L,J) = -R(L,J)
      END DO
      END IF
!
!     PRINT*, 'FORWARD TRANSFORM IN THE J DIRECTION'
!     DO I=1,L
!       PRINT*, (R(I,J),J=1,M)
!     END DO
!
      IF(IER1 .NE. 0 ) THEN
         IER=20
         CALL XERFFT('RFFT2F',-5)
         GO TO 100
      ENDIF
!
!
  100 CONTINUE
      END subroutine RFFT2F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFT2I (L, M, WSAVE, LENSAV, fft_sign, IER)
      SUBROUTINE RFFT2I (L, M, fft_sign, IER)
      integer     :: L, M
      integer     ::  IER
      integer     :: LENSAV
      INTEGER           :: LWSAV,MWSAV,MMSAV
      integer           :: IER1
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type),allocatable,dimension(:),intent(out) :: fft_sign
      
!
! INITIALIZE IER
!
      IER = 0
!
! VERIFY LENSAV
!
      LWSAV =   L  +INT(LOG(REAL(L,kind=rk))/LOG(2._rk))+4
      MWSAV =   2*M+INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
      MMSAV =   M  +INT(LOG(REAL(M,kind=rk))/LOG(2._rk))+4
!!      IF (LENSAV .LT. LWSAV+MWSAV+MMSAV) THEN
!!        IER = 2
!!        CALL XERFFT ('RFFT2I', 4)
!!        GO TO 100
!!      ENDIF

      allocate(fft_sign(3))
!
      !!CALL RFFTMI (L, WSAVE(1), LWSAV,fft_sign(1), IER1)
      CALL RFFTMI (L, fft_sign(1), IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('RFFT2I',-5)
        GO TO 100
      ENDIF

      !!CALL CFFTMI_orig (M, WSAVE(LWSAV+1),MWSAV,IER1)
      CALL CFFTMI(M, fft_sign(2))
     !! IF (IER1 .NE. 0) THEN
     !!   IER = 20
     !!   CALL XERFFT ('RFFT2I',-5)
     !! ENDIF
!
      !!CALL RFFTMI (M,WSAVE(LWSAV+MWSAV+1),MMSAV,fft_sign(3), IER1)
      CALL RFFTMI (M,fft_sign(3), IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('RFFT2I',-5)
        GO TO 100
      END IF
!
  100 CONTINUE
      END subroutine RFFT2I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE RFFTB1 (N,IN,C,CH,WA,FAC)
      integer           :: N,IN
      integer           :: NF,NA,IP,NL,KI,K1,L1,IX2,IX3,IX4
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IDO,IDL1,IW,L2
      integer           :: I,J,K               
      real(kind=rk)     :: HALF,HALFM 
      REAL(kind=rk)     :: CH(*), C(IN,*), WA(N)
      integer           :: FAC(*)
!
      NF = FAC(2)
      NA = 0
      DO 10 K1=1,NF
      IP = FAC(K1+2)
      NA = 1-NA
      IF(IP .LE. 5) GO TO 10
      IF(K1 .EQ. NF) GO TO 10
      NA = 1-NA
   10 CONTINUE
      HALF = .5_rk
      HALFM = -.5_rk
      MODN = MOD(N,2)
      NL = N-2
      IF(MODN .NE. 0) NL = N-1
      IF (NA .EQ. 0) GO TO 120
      CH(1) = C(1,1)
      CH(N) = C(1,N)
      DO 118 J=2,NL,2
        CH(J) = HALF*C(1,J)
        CH(J+1) = HALFM*C(1,J+1)
  118 CONTINUE
      GO TO 124
  120 DO 122 J=2,NL,2
        C(1,J) = HALF*C(1,J)
        C(1,J+1) = HALFM*C(1,J+1)
  122 CONTINUE
  124 L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL R1F4KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL R1F4KB (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
        CALL R1F2KB (IDO,L1,C,IN,CH,1,WA(IW))
         GO TO 105
  104    CALL R1F2KB (IDO,L1,CH,1,C,IN,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
        CALL R1F3KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2))
         GO TO 108
  107    CALL R1F3KB (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL R1F5KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL R1F5KB (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
        CALL R1FGKB (IDO,IP,L1,IDL1,C,C,C,IN,CH,CH,1,WA(IW))
         GO TO 114
  113    CALL R1FGKB (IDO,IP,L1,IDL1,CH,CH,CH,1,C,C,IN,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      END subroutine RFFTB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE RFFTF1 (N,IN,C,CH,WA,FAC)
      integer           :: N,IN
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      integer           :: IX2,IX3,IX4
      real(kind=rk)     :: HALF,HALFM 
      real(kind=rk)     :: TSN,SN,TSNM
      REAL(kind=rk)     :: CH(*),C(IN,*)
      !!REAL(kind=rk)     :: WA(N)
      REAL(kind=rk)     :: WA(:)
      integer           :: FAC(:)
!
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         !!write(*,*) 'K1,KH+3,IP,IDO,IW=',K1,KH+3,IP,IDO,IW
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL R1F4KF (IDO,L1,C,IN,CH,1,WA(IW:),WA(IX2:),WA(IX3:))
         GO TO 110
  101    continue
         CALL R1F4KF (IDO,L1,CH,1,C,IN,WA(IW:),WA(IX2:),WA(IX3:))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL R1F2KF (IDO,L1,C,IN,CH,1,WA(IW:))
         GO TO 110
  103    continue
         CALL R1F2KF (IDO,L1,CH,1,C,IN,WA(IW:))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL R1F3KF (IDO,L1,C,IN,CH,1,WA(IW:),WA(IX2:))
         GO TO 110
  105    continue
         CALL R1F3KF (IDO,L1,CH,1,C,IN,WA(IW:),WA(IX2:))
         GO TO 110
  106    continue
         IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL R1F5KF (IDO,L1,C,IN,CH,1,WA(IW:),WA(IX2:),WA(IX3:),WA(IX4:))
         GO TO 110
  107    continue
         CALL R1F5KF (IDO,L1,CH,1,C,IN,WA(IW:),WA(IX2:),WA(IX3:),WA(IX4:))
         GO TO 110
  108    continue
         IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL R1FGKF (IDO,IP,L1,IDL1,C,C,C,IN,CH,CH,1,WA(IW:))
         NA = 1
         GO TO 110
  109    continue
         CALL R1FGKF (IDO,IP,L1,IDL1,CH,CH,CH,1,C,C,IN,WA(IW:))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SN = 1._rk/real(N,kind=rk)
      TSN = 2._rk/real(N,kind=rk)
      TSNM = -TSN
      MODN = MOD(N,2)
      NL = N-2
      IF(MODN .NE. 0) NL = N-1
      IF (NA .NE. 0) GO TO 120
      C(1,1) = SN*CH(1)
      DO 118 J=2,NL,2
        C(1,J) = TSN*CH(J)
        C(1,J+1) = TSNM*CH(J+1)
  118 CONTINUE
      IF(MODN .NE. 0) RETURN

      C(1,N) = SN*CH(N)
      RETURN

  120 C(1,1) = SN*C(1,1)
      DO 122 J=2,NL,2
        C(1,J) = TSN*C(1,J)
        C(1,J+1) = TSNM*C(1,J+1)
  122 CONTINUE
      IF(MODN .NE. 0) RETURN

      C(1,N) = SN*C(1,N)
      END subroutine RFFTF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE RFFTI1 (N,WA,FAC)
      integer           :: N,IN
      REAL(kind=rk)     :: WA(N),FAC(*)
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IPM
      real(kind=dk)     :: FI
      real(kind=dk)     :: TPI,ARGH,ARGLD,ARG
      integer           :: NTRY
      INTEGER           :: NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
!
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF


      !!TPI = 8._dk*DATAN(1._dk)
      TPI=2._dk*PI_long
      ARGH = TPI/REAL(N,kind=dk)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = REAL(LD,kind=dk)*ARGH
            FI = 0._dk
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1._dk
               ARG = FI*ARGLD
              !!WA(I-1) = DCOS(ARG)
              !!WA(I) = DSIN(ARG)
              WA(I-1) = COS(ARG)
              WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      END subroutine RFFTI1
!******************************************************************************************************************
      SUBROUTINE factor_RFFTI1 (N,NF,FAC)
      integer           :: N,IN
      integer,allocatable,dimension(:),intent(out) :: FAC
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IPM
      real(kind=dk)     :: FI
      real(kind=dk)     :: TPI,ARGH,ARGLD,ARG
      integer           :: NTRY
      INTEGER           :: NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      integer           :: ibase2
      integer           :: ll
!
      ibase2 = NINT(LOG(real(N,kind=dk))/LOG(2.0_dk))+1

      allocate(FAC(ibase2+2))
      write(*,'(*(a,i0))') 'factor_RFFTI1: N ',N,' ibase2 ',ibase2
      if(N==1) then

      NF = 0

      else

      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
 !! write(*,'(*(a,i0))') 'J ',J,' NQ ',NQ,' NR ',NR,' NF ',NF
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      endif
      FAC(1) = N   ! factored number
      FAC(2) = NF  ! count of factors

      FAC=FAC(1:NF+2)

         write(*,*)
         write(*,'(*(a,i0))') 'factor_RFFTI1: ',nf,' prime factors of ',n
      do ll=2+1,NF+2 
         write(*,'(i5)',advance='no') fac(ll)
      enddo
         write(*,*)

end subroutine factor_RFFTI1
!******************************************************************************************************************
      SUBROUTINE size_WA_RFFTI1 (N,NF,FAC,size_WA)
      integer           :: N,IN
      integer           :: FAC(*)
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IPM
      real(kind=dk)     :: FI
      real(kind=dk)     :: TPI,ARGH,ARGLD,ARG
      integer           :: NTRY
      INTEGER           :: NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      integer           :: size_WA
!

      size_WA=0
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      do K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         do J=1,IPM
            LD = LD+L1
            I = IS
            DO II=3,IDO,2
               I = I+2
               size_WA=max(size_WA,I)
            enddo   
            IS = IS+IDO
         enddo
         L1 = L2
      enddo

      END subroutine size_WA_RFFTI1
!******************************************************************************************************************
      SUBROUTINE WA_RFFTI1 (N,NF,WA,FAC)
      integer           :: N,IN
      REAL(kind=rk)     :: WA(N)
      integer           :: FAC(*)
      integer           :: IW,IDO,IPP2,IPPH,IS,IDIJ,I,IB,IC,II,IP,IM,I1,I2,J,JC,J1,J2,K,KC,KH,KI,K1,K2,L,LC,LD,L1,L2      
      integer           :: M,M1,M2,M2S,M1D,NS2,NF,NA,NL,NR,NQ,NP1,NP2,NM1,NFM1,IDL1,LJ
      integer           :: IPM
      real(kind=dk)     :: FI
      real(kind=dk)     :: TPI,ARGH,ARGLD,ARG
      integer           :: NTRY
      INTEGER           :: NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
!


      !!TPI = 8._dk*DATAN(1._dk)
      TPI=2._dk*PI_long
      ARGH = TPI/REAL(N,kind=dk)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
           !! write(*,'(*(a,i0))') 'K1 ',K1,' FAC(K1+2) ',IP 
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            !!write(*,'(*(a,i0))') '          J ',J,' IS ',IS,' LD ',LD
            ARGLD = REAL(LD,kind=dk)*ARGH
            FI = 0._dk
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1._dk
               ARG = FI*ARGLD
              !!WA(I-1) = DCOS(ARG)
              !!WA(I) = DSIN(ARG)
              WA(I-1) = COS(ARG)
              WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
           !! write(*,'(*(a,i0))') '                         K1 ',K1,' IP ',IP,' IDO ',IDO,' IS ',IS
  109    CONTINUE
         L1 = L2
  110 CONTINUE

      END subroutine WA_RFFTI1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFTMB (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV,WORK, LENWRK,fft_sign, IER)
      SUBROUTINE RFFTMB (LOT, JUMP, N, INC, R, LENR, WORK, LENWRK,fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENR
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: R(LENR)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL  XERCON
!
      IER = 0
!
      IF (LENR .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('RFFTMB ', 6)
!!      ELSEIF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('RFFTMB ', 8)
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('RFFTMB ', 10)
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('RFFTMB ', -1)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!CALL MRFTB1 (LOT,JUMP,N,INC,R,WORK,WSAVE,WSAVE(N+1))
      CALL MRFTB1 (LOT,JUMP,N,INC,R,WORK,fft_sign%WA,fft_sign%FAC)
      END subroutine RFFTMB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFTMF (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      SUBROUTINE RFFTMF (LOT, JUMP, N, INC, R, LENR, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENR
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: R(LENR)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL  XERCON
!
      IER = 0
!
      IF (LENR .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('RFFTMF ', 6)
!!      ELSEIF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('RFFTMF ', 8)
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('RFFTMF ', 10)
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('RFFTMF ', -1)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      !!CALL MRFTF1 (LOT,JUMP,N,INC,R,WORK,WSAVE,WSAVE(N+1))
      CALL MRFTF1 (LOT,JUMP,N,INC,R,WORK,fft_sign%WA,fft_sign%FAC)
      END subroutine RFFTMF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE RFFTMI (N, WSAVE, LENSAV,fft_sign, IER)
      SUBROUTINE RFFTMI (N, fft_sign, IER)
      integer     :: N
      integer     ::  IER
      integer     :: LENSAV
      !!REAL(kind=rk)            :: WSAVE(LENSAV)
      integer                  :: NF
      integer                  :: size_WA
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
 !!     IF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
 !!       IER = 2
 !!       CALL XERFFT ('RFFTMI ', 3)
 !!     ENDIF
!
  !!    if( fft_sign%switch == 2 .and. fft_sign%N == N .and. allocated(fft_sign%FAC) .and. allocated(fft_sign%WA)) then
  !!         goto 1111
  !!    endif

      fft_sign%switch = 2

      fft_sign%N = N

      IF (N .EQ. 1) RETURN
!
      call factor_RFFTI1 (N,NF,fft_sign%FAC)

      call size_WA_RFFTI1 (N,NF,fft_sign%FAC,size_WA)

      if(allocated(fft_sign%WA) ) deallocate(fft_sign%WA)
      allocate(fft_sign%WA(size_WA))

      call WA_RFFTI1 (N,NF,fft_sign%WA,fft_sign%FAC)
!
      !!CALL MRFTI1 (N,WSAVE(1),WSAVE(N+1))
      END subroutine RFFTMI
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE RFFTMI_orig (N, WSAVE, LENSAV, IER)
      integer     :: N
      integer    ::  IER
      integer   :: LENSAV
      REAL(kind=rk)     :: WSAVE(LENSAV)
!
      IER = 0
!
      IF (LENSAV .LT. N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
        IER = 2
        CALL XERFFT ('RFFTMI ', 3)
      ENDIF
!
      IF (N .EQ. 1) RETURN
!
      CALL MRFTI1 (N,WSAVE(1),WSAVE(N+1))
      END subroutine RFFTMI_orig
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!************************************************************************************************************************
      !!SUBROUTINE SINQ1B ( N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK,fft_sign, IER)
      !!SUBROUTINE SINQ1B ( N, INC, X, LENX, WORK, LENWRK,fft_sign, IER)
      SUBROUTINE SINQ1B ( N, INC, X, LENX,fft_sign, IER)
      integer     :: N, INC, LENX
      !!integer     ::  LENWRK
      integer     ::  IER
      integer     :: LENSAV
      integer           :: LJ,K,KC,M,NS2,IER1
      real(kind=rk)     :: XHOLD
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINQ1B', 6)
   !!   ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
   !!     IER = 2
   !!     CALL XERFFT ('SINQ1B', 8)
    !!  ELSEIF (LENWRK .LT. N) THEN
    !!    IER = 3
    !!    CALL XERFFT ('SINQ1B', 10)
      ENDIF
!
      IF (N .GT. 1) GO TO 101
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         X(1,K) = -X(1,K)
  102 CONTINUE
      !!CALL COSQ1B (N,INC,X,LENX,WSAVE,LENSAV,WORK,LENWRK,fft_sign,IER1)
      !!CALL COSQ1B (N,INC,X,LENX,WORK,LENWRK,fft_sign,IER1)
      CALL COSQ1B (N,INC,X,LENX,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINQ1B',-5)
        GO TO 300
      ENDIF
      DO 103 K=1,NS2
         KC = N-K
         XHOLD = X(1,K)
         X(1,K) = X(1,KC+1)
         X(1,KC+1) = XHOLD
  103 CONTINUE
  300 continue
      END subroutine SINQ1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINQ1F ( N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK,fft_sign, IER)
      !!SUBROUTINE SINQ1F ( N, INC, X, LENX,WORK, LENWRK,fft_sign, IER)
      SUBROUTINE SINQ1F ( N, INC, X, LENX,fft_sign, IER)
      integer     :: N, INC, LENX
      !!integer     ::  LENWRK
      integer     ::  IER
      integer     :: LENSAV
      integer           :: LJ,K,KC,M,NS2,IER1
      real(kind=rk)     :: XHOLD
      !!REAL(kind=rk)     ::  WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINQ1F', 6)
        GO TO 300
  !!    ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
  !!      IER = 2
  !!      CALL XERFFT ('SINQ1F', 8)
  !!      GO TO 300
   !!   ELSEIF (LENWRK .LT. N) THEN
   !!     IER = 3
   !!     CALL XERFFT ('SINQ1F', 10)
   !!     GO TO 300
      ENDIF
!
      IF (N .EQ. 1) RETURN
      NS2 = N/2
      DO 101 K=1,NS2
         KC = N-K
         XHOLD = X(1,K)
         X(1,K) = X(1,KC+1)
         X(1,KC+1) = XHOLD
  101 CONTINUE

      !!CALL COSQ1F (N,INC,X,LENX,WSAVE,LENSAV,WORK,LENWRK,fft_sign,IER1)
      !!CALL COSQ1F (N,INC,X,LENX,WORK,LENWRK,fft_sign,IER1)
      CALL COSQ1F (N,INC,X,LENX,fft_sign,IER1)

      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINQ1F',-5)
        GO TO 300
      ENDIF
      DO 102 K=2,N,2
         X(1,K) = -X(1,K)
  102 CONTINUE
  300 continue
      END subroutine SINQ1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINQ1I (N, WSAVE, LENSAV,fft_sign, IER)
      SUBROUTINE SINQ1I (N, fft_sign, IER)
      integer     :: N
      integer     ::  IER
      integer     :: LENSAV
      integer           :: LJ,K,KC,M,NS2,IER1
      real(kind=rk)     :: XHOLD
      real(kind=rk)     :: PIH,DT,FK
      integer     :: LENWRK
  !!    REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
!!      IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('SINQ1I', 3)
!!        GO TO 300
!!      ENDIF
!
      !!CALL COSQ1I (N, WSAVE, LENSAV,fft_sign, IER1)
      !!CALL COSQ1I (N, fft_sign, IER1)

      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(N))

      !!PIH = 2._rk*ATAN(1._rk)
      PIH = 0.5_dk*PI_long
      DT = PIH/REAL(N,kind=rk)
      FK = 0._rk
      DO 101 K=1,N
         FK = FK+1._rk
        !! WSAVE(K) = COS(FK*DT)
         fft_sign%WSAVE(K)=COS(FK*DT)
  101 CONTINUE

      !!CALL RFFT1I_orig (N, WSAVE(N+1), LNSV, IER1)
      CALL all_FFT1I(N, fft_sign)
      fft_sign%sign_type=SINQ1_sign_type

      LENWRK=N
      if(allocated(fft_sign%work) ) deallocate(fft_sign%work)
      allocate(fft_sign%work(LENWRK))

  300 continue
      END subroutine SINQ1I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINQMB (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      SUBROUTINE SINQMB (LOT, JUMP, N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer     ::  LENWRK, IER
      integer     :: LENSAV
      integer           :: LJ,K,KC,M,NS2,IER1
      real(kind=rk)     :: XHOLD
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINQMB', 6)
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('SINQMB', 8)
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('SINQMB', 10)
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('SINQMB', -1)
      ENDIF
!
      LJ = (LOT-1)*JUMP+1
      IF (N .GT. 1) GO TO 101
      DO 201 M=1,LJ,JUMP
         X(M,1) = 4._rk*X(M,1)
 201  CONTINUE
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         DO 202 M=1,LJ,JUMP
         X(M,K) = -X(M,K)
 202     CONTINUE
  102 CONTINUE
      !!CALL COSQMB (LOT,JUMP,N,INC,X,LENX,WSAVE,LENSAV,WORK,LENWRK, fft_sign,IER1)
      CALL COSQMB (LOT,JUMP,N,INC,X,LENX,WORK,LENWRK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINQMB',-5)
        GO TO 300
      ENDIF
      DO 103 K=1,NS2
         KC = N-K
         DO 203 M=1,LJ,JUMP
         XHOLD = X(M,K)
         X(M,K) = X(M,KC+1)
         X(M,KC+1) = XHOLD
 203     CONTINUE
  103 CONTINUE
  300 CONTINUE
      END subroutine SINQMB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINQMF (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      SUBROUTINE SINQMF (LOT, JUMP, N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer     ::  LENWRK, IER
      integer     :: LENSAV
      integer           :: LJ,K,KC,M,NS2,IER1
      real(kind=rk)     :: XHOLD
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINQMF', 6)
        GO TO 300
!!      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('SINQMF', 8)
!!        GO TO 300
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('SINQMF', 10)
        GO TO 300
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('SINQMF', -1)
        GO TO 300
      ENDIF
!
      IF (N .EQ. 1) RETURN
      NS2 = N/2
      LJ = (LOT-1)*JUMP+1
      DO 101 K=1,NS2
         KC = N-K
         DO 201 M=1,LJ,JUMP
         XHOLD = X(M,K)
         X(M,K) = X(M,KC+1)
         X(M,KC+1) = XHOLD
 201     CONTINUE
  101 CONTINUE
      !!CALL COSQMF (LOT,JUMP,N,INC,X,LENX,WSAVE,LENSAV,WORK,LENWRK, fft_sign,IER1)
      CALL COSQMF (LOT,JUMP,N,INC,X,LENX,WORK,LENWRK, fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINQMF',-5)
        GO TO 300
      ENDIF
      DO 102 K=2,N,2
         DO 202 M=1,LJ,JUMP
         X(M,K) = -X(M,K)
 202     CONTINUE
  102 CONTINUE
  300 CONTINUE
      END subroutine SINQMF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINQMI (N, WSAVE, LENSAV, fft_sign, IER)
      SUBROUTINE SINQMI (N, fft_sign, IER)
      integer     :: N
      integer     ::  IER
      integer     :: LENSAV
      integer           :: I,K,KC,NP1,NS2,IW1,IW2,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,T1,T2,XHOLD,FNP1S4,DT,PI
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
     !! IF (LENSAV .LT. 2*N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
     !!   IER = 2
     !!   CALL XERFFT ('SINQMI', 3)
     !!   GO TO 300
     !! ENDIF
!
      !!CALL COSQMI (N, WSAVE, LENSAV, fft_sign, IER1)
      CALL COSQMI (N, fft_sign, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINQMI',-5)
      ENDIF
  300 CONTINUE
      END subroutine SINQMI
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINT1B ( N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE SINT1B ( N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE SINT1B ( N, INC, X, LENX, fft_sign, IER)
      integer     :: N, INC
      !!integer     ::  LENWRK
      integer     ::  IER
      integer     :: LENSAV
      integer           :: I,K,KC,NP1,NS2,IW1,IW2,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,T1,T2,XHOLD,FNP1S4,DT,PI
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
      if( debug > 0 ) write(*,*) 'SINT1B start'
      IER = 0
!
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINT1B', 6)
        GO TO 100
  !!    ELSEIF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
  !!      IER = 2
  !!      CALL XERFFT ('SINT1B', 8)
  !!      GO TO 100
    !!  ELSEIF (LENWRK .LT. (2*N+2)) THEN
    !!    IER = 3
    !!    CALL XERFFT ('SINT1B', 10)
    !!    GO TO 100
      ENDIF
!
      !!CALL SINTB1(N,INC,X,WSAVE,WORK,WORK(N+2),fft_sign,IER1)
      !!CALL SINTB1(N,INC,X,WORK,WORK(N+2),fft_sign,IER1)
      CALL SINTB1(N,INC,X,fft_sign%WORK_a,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINT1B',-5)
      ENDIF
!
  100 CONTINUE
  if( debug > 0 ) write(*,*) 'SINT1B end'
      END subroutine SINT1B
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINT1F ( N, INC, X, LENX, WSAVE, LENSAV,WORK, LENWRK, fft_sign, IER)
      !!SUBROUTINE SINT1F ( N, INC, X, LENX, WORK, LENWRK, fft_sign, IER)
      SUBROUTINE SINT1F ( N, INC, X, LENX, fft_sign, IER)
      integer     :: N, INC
      !!integer     ::  LENWRK
      integer     ::  IER
      !!integer     :: LENSAV
      integer           :: I,K,KC,NP1,NS2,IW1,IW2,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,T1,T2,XHOLD,FNP1S4,DT,PI
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      !!real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      type(fft_real_sign_type) :: fft_sign
!
  if( debug > 0 ) write(*,*) 'SINT1F start'
      IER = 0
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINT1F', 6)
        GO TO 100
 !!     ELSEIF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
 !!       IER = 2
 !!       CALL XERFFT ('SINT1F', 8)
 !!       GO TO 100
  !!    ELSEIF (LENWRK .LT. (2*N+2)) THEN
  !!      IER = 3
  !!      CALL XERFFT ('SINT1F', 10)
  !!      GO TO 100
      ENDIF
!
      !!CALL SINTF1(N,INC,X,WSAVE,WORK,WORK(N+2),fft_sign,IER1)
      !!CALL SINTF1(N,INC,X,WORK,WORK(N+2),fft_sign,IER1)
      !!write(*,*) 'SINT1F: allocated(fft_sign%work_a)=',allocated(fft_sign%work_a)
      !!write(*,*) 'SINT1F: size(fft_sign%work_a)=',size(fft_sign%work_a)
      CALL SINTF1(N,INC,X,fft_sign%WORK_a,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINT1F',-5)
      ENDIF
  100 CONTINUE
  if( debug > 0 ) write(*,*) 'SINT1F end'
      END subroutine SINT1F
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINT1I (N, WSAVE, LENSAV, fft_sign, IER)
      SUBROUTINE SINT1I (N, fft_sign, IER)
      integer     :: N
      integer    ::  IER
      !!integer   :: LENSAV
      integer           :: I,K,KC,NP1,NS2,IW1,IW2,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,T1,T2,XHOLD,FNP1S4
      real(kind=rk)     :: PI,DT,FK
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
      integer                  :: size_WA
      integer     :: LENWRK
      integer                  :: NF
!
   !!   if( fft_sign%switch == 2 .and. fft_sign%N == N .and. allocated(fft_sign%FAC) .and. allocated(fft_sign%WA)) then
   !!        goto 1111
   !!   endif

      fft_sign%switch = 2

      fft_sign%N = N
      IER = 0
!
   !!   IF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
   !!     IER = 2
   !!     CALL XERFFT ('SINT1I', 3)
   !!     GO TO 300
   !!   ENDIF
!
      !!PI = 4._rk*ATAN(1._rk)
      PI = PI_long
      !!IF (N .LE. 1) RETURN
      NS2 = N/2
      NP1 = N+1

    !!  LNSV = NP1 + INT(LOG(REAL(NP1))/LOG(2._rk)) +4
    !!  size_WA=LNSV

      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(NS2))


      DT = PI/REAL(NP1,kind=rk)
      DO 101 K=1,NS2
     !!    WSAVE(K) = 2._rk*SIN(K*DT)
         fft_sign%WSAVE(K) = 2._rk*SIN(K*DT)
  101 CONTINUE

      !!CALL RFFT1I_orig (NP1, WSAVE(NS2+1), LNSV, IER1)
      CALL all_FFT1I(NP1, fft_sign)
      fft_sign%sign_type=SINT1_sign_type

      LENWRK=2*N+2
      LENWRK=N+1
      if(allocated(fft_sign%work) ) deallocate(fft_sign%work)
      allocate(fft_sign%work(LENWRK))

      if(allocated(fft_sign%work_a) ) deallocate(fft_sign%work_a)
      allocate(fft_sign%work_a(LENWRK))
      !!write(*,*) 'SINT1I: size(fft_sign%work_a)=',size(fft_sign%work_a)

!
  300 CONTINUE
  1111 continue
  if( debug > 0 ) write(*,*) 'SINT1I: end'
      END subroutine SINT1I
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINTB1(N,INC,X,WSAVE,XH,WORK,fft_sign,IER)
      !!SUBROUTINE SINTB1(N,INC,X,XH,WORK,fft_sign,IER)
      SUBROUTINE SINTB1(N,INC,X,XH,fft_sign,IER)
      INTEGER           :: LOT,JUMP,N,INC,IER
      integer           :: I,K,KC,NP1,NS2,IW1,IW2,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: SRT3S2,SSQRT3,SFNP1,XHOLD,FNP1S4,DT
      !!REAL(kind=rk)     :: WSAVE(*)
      REAL(kind=rk)     :: X(INC,*),XH(*)
      !!real(kind=rk)     :: WORK(*)    ! Kuester
      real(kind=rk)     :: T1,T2
      real(kind=dk)    :: DSUM
      type(fft_real_sign_type) :: fft_sign


      IER = 0
      IF (N-2) 200,102,103
  102 SRT3S2 = SQRT(3._rk)/2._rk
      XHOLD = SRT3S2*(X(1,1)+X(1,2))
      X(1,2) = SRT3S2*(X(1,1)-X(1,2))
      X(1,1) = XHOLD
      GO TO 200
  103 NP1 = N+1
      NS2 = N/2
      DO 104 K=1,NS2
         KC = NP1-K
         T1 = X(1,K)-X(1,KC)
         T2 = fft_sign%WSAVE(K)*(X(1,K)+X(1,KC))
         XH(K+1) = T1+T2
         XH(KC+1) = T2-T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
      XH(NS2+2) = 4._rk*X(1,NS2+1)
  124 XH(1) = 0._rk
      LNXH = NP1
      !!LNSV = NP1 + INT(LOG(REAL(NP1))/LOG(2._rk)) + 4
      LNWK = NP1
!
      !!CALL RFFT1F(NP1,1,XH,LNXH,WSAVE(NS2+1),LNSV,WORK,LNWK,fft_sign,IER1)
      !!CALL RFFT1F(NP1,1,XH,LNXH,WORK,LNWK,fft_sign,IER1)
      CALL RFFT1F(NP1,1,XH,LNXH,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINTB1',-5)
        GO TO 200
      ENDIF
!
      IF(MOD(NP1,2) .NE. 0) GO TO 30
      XH(NP1) = XH(NP1)+XH(NP1)
 30   FNP1S4 = REAL(NP1,kind=rk)/4._rk
         X(1,1) = FNP1S4*XH(1)
         DSUM = X(1,1)
      DO 105 I=3,N,2
            X(1,I-1) = FNP1S4*XH(I)
            DSUM = DSUM+FNP1S4*XH(I-1)
            X(1,I) = DSUM
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 200
         X(1,N) = FNP1S4*XH(N+1)
!
  200 CONTINUE
      END subroutine SINTB1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINTF1(N,INC,X,WSAVE,XH,WORK,fft_sign,IER)
      !!SUBROUTINE SINTF1(N,INC,X,XH,WORK,fft_sign,IER)
      SUBROUTINE SINTF1(N,INC,X,XH,fft_sign,IER)
      INTEGER           :: N,INC,IER
      integer           :: I,K,KC,NP1,NS2,IW1,IW2,IER1
      real(kind=rk)     :: SSQRT3,SFNP1,XHOLD,FNP1S4
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      !!REAL(kind=rk)     :: WSAVE(*)
      REAL(kind=rk)     :: X(INC,*),XH(*)
      !!real(kind=rk)     :: WORK(*)      ! Kuester
      real(kind=rk)     :: T1,T2
      real(kind=dk)    :: DSUM
      type(fft_real_sign_type) :: fft_sign

      IER = 0
      IF (N-2) 200,102,103
  102 SSQRT3 = 1._rk/SQRT(3._rk)
      XHOLD = SSQRT3*(X(1,1)+X(1,2))
      X(1,2) = SSQRT3*(X(1,1)-X(1,2))
      X(1,1) = XHOLD
      GO TO 200
  103 NP1 = N+1
      NS2 = N/2
      DO 104 K=1,NS2
         KC = NP1-K
         T1 = X(1,K)-X(1,KC)
         T2 = fft_sign%WSAVE(K)*(X(1,K)+X(1,KC))
         XH(K+1) = T1+T2
         XH(KC+1) = T2-T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) GO TO 124
      XH(NS2+2) = 4._rk*X(1,NS2+1)
  124 XH(1) = 0._rk
      LNXH = NP1
      LNSV = NP1 + INT(LOG(REAL(NP1))/LOG(2._rk)) + 4
      LNWK = NP1
!
      !!CALL RFFT1F(NP1,1,XH,LNXH,WSAVE(NS2+1),LNSV,WORK,LNWK,fft_sign,IER1)
      !!CALL RFFT1F(NP1,1,XH,LNXH,WORK,LNWK,fft_sign,IER1)
      CALL RFFT1F(NP1,1,XH,LNXH,fft_sign,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINTF1',-5)
        GO TO 200
      ENDIF
!
      IF(MOD(NP1,2) .NE. 0) GO TO 30
      XH(NP1) = XH(NP1)+XH(NP1)
   30 SFNP1 = 1._rk/REAL(NP1,kind=rk)
         X(1,1) = .5_rk*XH(1)
         DSUM = X(1,1)
      DO 105 I=3,N,2
            X(1,I-1) = .5_rk*XH(I)
            DSUM = DSUM+.5_rk*XH(I-1)
            X(1,I) = DSUM
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 200
      X(1,N) = .5_rk*XH(N+1)
  200 continue
      END subroutine SINTF1
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      ! Entrance
      !!SUBROUTINE SINTMB (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, DSUM, LENWRK, IER)
      SUBROUTINE SINTMB (LOT, JUMP, N, INC, X, LENX, WORK, DSUM, LENWRK, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IW1,IW2,IER1
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      real(kind=dk)    :: DSUM(*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINTMB', 6)
        GO TO 100
!!      ELSEIF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('SINTMB', 8)
!!        GO TO 100
      ELSEIF (LENWRK .LT. LOT*(2*N+4)) THEN
        IER = 3
        CALL XERFFT ('SINTMB', 10)
        GO TO 100
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('SINTMB', -1)
        GO TO 100
      ENDIF
!
      IW1 = LOT+LOT+1
      IW2 = IW1+LOT*(N+1)
      !Kuester CALL MSNTB1(LOT,JUMP,N,INC,X,WSAVE,WORK,WORK(IW1),WORK(IW2),IER1)
      !!CALL MSNTB1(LOT,JUMP,N,INC,X,WSAVE,DSUM,WORK(IW1),WORK(IW2), fft_sign, IER1)
      CALL MSNTB1(LOT,JUMP,N,INC,X,DSUM,WORK(IW1),WORK(IW2), fft_sign, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINTMB',-5)
      ENDIF
!
  100 CONTINUE
      END subroutine SINTMB
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      ! Entrance
      !!SUBROUTINE SINTMF (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV,WORK, DSUM, LENWRK, IER)
      SUBROUTINE SINTMF (LOT, JUMP, N, INC, X, LENX, WORK, DSUM, LENWRK, IER)
      integer     :: LOT, JUMP, N, INC, LENX
      integer    ::  LENWRK, IER
      !!integer   :: LENSAV
      integer           :: IW1,IW2,IER1
!!      real(kind=rk)     :: WSAVE(LENSAV)
      real(kind=rk)     :: WORK(LENWRK)
      real(kind=rk)     :: X(INC,*)
      real(kind=dk)    :: DSUM(*)
      type(fft_real_sign_type) :: fft_sign
!!      LOGICAL    XERCON
!
      IER = 0
!
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('SINTMF', 6)
        GO TO 100
!!      ELSEIF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
!!        IER = 2
!!        CALL XERFFT ('SINTMF', 8)
!!        GO TO 100
      ELSEIF (LENWRK .LT. LOT*(2*N+4)) THEN
        IER = 3
        CALL XERFFT ('SINTMF', 10)
        GO TO 100
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('SINTMF', -1)
        GO TO 100
      ENDIF
!
      IW1 = LOT+LOT+1
      IW2 = IW1+LOT*(N+1)
      !Kuester CALL MSNTF1(LOT,JUMP,N,INC,X,WSAVE,WORK,WORK(IW1),WORK(IW2),IER1)
      !!CALL MSNTF1(LOT,JUMP,N,INC,X,WSAVE,DSUM,WORK(IW1),WORK(IW2), fft_sign, IER1)
      CALL MSNTF1(LOT,JUMP,N,INC,X,DSUM,WORK(IW1),WORK(IW2), fft_sign, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINTMF',-5)
      ENDIF
  100 CONTINUE
      END subroutine SINTMF
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      !!SUBROUTINE SINTMI (N, WSAVE, LENSAV, fft_sign, IER)
      SUBROUTINE SINTMI (N, fft_sign, IER)
      integer     :: N
      integer     ::  IER
      integer     :: LENSAV
      integer           :: K,NS2,NP1,IER1
      integer           :: IDP2,LENX,LNSV,LNWK,LNXH,MODN
      real(kind=rk)     :: DT,PI
      !!REAL(kind=rk)     :: WSAVE(LENSAV)
      type(fft_real_sign_type) :: fft_sign
!
      IER = 0
!
 !!     IF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N,kind=rk))/LOG(2._rk)) +4) THEN
 !!       IER = 2
 !!       CALL XERFFT ('SINTMI', 3)
 !!       GO TO 300
 !!     ENDIF
!
      !!PI = 4._rk*ATAN(1._rk)
      PI = PI_long
      IF (N .LE. 1) RETURN
      NS2 = N/2
      NP1 = N+1

      if(allocated(fft_sign%WSAVE)) deallocate(fft_sign%WSAVE)
      allocate(fft_sign%WSAVE(NS2))

      DT = PI/REAL(NP1,kind=rk)
      DO 101 K=1,NS2
         !!WSAVE(K) = 2._rk*SIN(K*DT)
         fft_sign%WSAVE(K) = 2._rk*SIN(K*DT)
  101 CONTINUE

      !!LNSV = NP1 + INT(LOG(REAL(NP1))/LOG(2._rk)) +4
      !!CALL RFFTMI (NP1, WSAVE(NS2+1), LNSV,fft_sign, IER1)
      CALL RFFTMI (NP1,fft_sign, IER1)

      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINTMI',-5)
      ENDIF
!
  300 CONTINUE
      END subroutine SINTMI
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE TABLES (IDO,IP,WA)
      integer           :: IDO,IP
      integer           :: I,J
      real(kind=rk)     :: TPI,ARGZ,ARG1,ARG2,ARG3,ARG4
      REAL(kind=rk)     :: WA(IDO,IP-1,2)
!
      !!TPI = 8._dk*ATAN(1._dk)
      TPI=2._dk*PI_long
      ARGZ = TPI/REAL(IP,kind=rk)
      ARG1 = TPI/REAL(IDO*IP,kind=rk)
      DO 110 J=2,IP
         ARG2 = REAL(J-1,kind=rk)*ARG1
         DO 100 I=1,IDO
            ARG3 = REAL(I-1,kind=rk)*ARG2
            WA(I,J-1,1) = COS(ARG3)
            WA(I,J-1,2) = SIN(ARG3)
  100    CONTINUE
         IF (IP .LE. 5) GO TO 110
         ARG4 = REAL(J-1,kind=rk)*ARGZ
         WA(1,J-1,1) = COS(ARG4)
         WA(1,J-1,2) = SIN(ARG4)
  110 CONTINUE
      END subroutine TABLES
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      subroutine w2r(ldr,ldw,l,m,r,w)
      integer           :: ldr,ldw,l,m
      integer           :: i,j 
      real(kind=rk)     ::  r(ldr,*),w(ldw,*)
      do j=1,m
      do i=1,l
      r(i,j) = w( i,j)
      end do
      end do
      return
      end subroutine w2r
!************************************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      LOGICAL FUNCTION XERCON (INC,JUMP,N,LOT)
      INTEGER    :: INC, JUMP, N, LOT
      INTEGER    :: I, J, JNEW, LCM
!
!     Definition: positive integers INC, JUMP, N and LOT are consistent
!                                                            ----------
!     if I1*INC + J1*JUMP = I2*INC + J2*JUMP for I1,I2 < N and J1,J2
!     < LOT implies I1=I2 and J1=J2;
!
!     For multiple FFTs to execute correctly, input parameters INC,
!     JUMP, N and LOT must be consistent ... otherwise at least one
!     array element mistakenly is transformed more than once.
!
!     XERCON = .TRUE. if and only if INC, JUMP, N and LOT are
!     consistent.
!
!     ------------------------------------------------------------------
!
!     Compute I = greatest common divisor (INC, JUMP)
!
      I = INC
      J = JUMP
   10 CONTINUE
      IF (J .NE. 0) THEN
        JNEW = MOD(I,J)
        I    = J
        J    = JNEW
        GO TO 10
      ENDIF
!
! Compute LCM = least common multiple (INC, JUMP)
!
      LCM = (INC*JUMP)/I
!
! Check consistency of INC, JUMP, N, LOT
!
      IF (LCM .LE. (N-1)*INC .AND. LCM .LE. (LOT-1)*JUMP) THEN
        XERCON = .FALSE.
      ELSE
        XERCON = .TRUE.
      ENDIF
!
      END FUNCTION XERCON
!****************************************************************************************************
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
      SUBROUTINE XERFFT( SRNAME, INFO)
!
!     .. Scalar Arguments ..
      CHARACTER(len=*)  :: SRNAME
      INTEGER           :: INFO
!
!     ..
!
!  Purpose
!  =======
!
!  XERFFT  is an error handler for library FFTPACK version 5.1 routines.
!  It is called by an FFTPACK 5.1 routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER(len=*)
!          The name of the routine which called XERFFT.
!
!  INFO    (input) INTEGER
!          When a single  invalid parameter in the parameter list of
!          the calling routine has been detected, INFO is the position
!          of that parameter.  In the case when an illegal combination
!          of LOT, JUMP, N, and INC has been detected, the calling
!          subprogram calls XERFFT with INFO = -1;
!
! =====================================================================
!
!     .. Executable Statements ..
!
      IF (INFO .GE. 1) THEN
        WRITE( *, '(A,A,A,I3,A)') ' ** On entry to ', SRNAME,' parameter number ', INFO, ' had an illegal value'
      ELSEIF (INFO .EQ. -1) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,' parameters LOT, JUMP, N and INC are inconsistent'
      ELSEIF (INFO .EQ. -2) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,' parameter L is greater than LDIM'
      ELSEIF (INFO .EQ. -3) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,' parameter M is greater than MDIM'
      ELSEIF (INFO .EQ. -5) THEN
        WRITE( *, '(A,A,A,A)') ' ** Within ', SRNAME,' input error returned by lower level routine'
      ELSEIF (INFO .EQ. -6) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,' parameter LDIM is less than 2*(L/2+1)'
      ENDIF
!
      STOP
!
!     End of XERFFT
!
      END subroutine XERFFT
!************************************************************************************************************************
end module fftpack_module
