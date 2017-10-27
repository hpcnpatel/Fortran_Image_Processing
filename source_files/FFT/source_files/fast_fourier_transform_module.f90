module fast_fourier_transform_module
use constants_module
use fftpack_module
implicit none
contains
!************************************************************************************************************************
subroutine fft_complex_1D(input,output,forward_backward_switch)
      integer                                              :: forward_backward_switch
      integer                                              :: INC,LENC,LENSAV,LENWRK
      integer                                              :: IER
      integer                                              :: N
      complex(kind=rk),dimension(:)                           :: input
      complex(kind=rk),dimension(:),allocatable,intent(out)   :: output
      real(kind=rk),dimension(:),allocatable               :: C_R         
      real(kind=rk),dimension(:),allocatable               :: WSAVE,WORK  
      integer                                              :: ll
      type(fft_real_sign_type)    :: fft_sign
!
      INC=1 ! stride of data in input
      N=size(input)
      LENSAV= 2*N + INT(LOG(REAL(N))/LOG(2.)) + 4
      LENWRK=2*N
      LENC=INC*(N-1) + 1
      allocate(C_R(2*N))
      allocate(output(N))
      allocate(WSAVE(LENSAV))
      allocate(WORK(LENWRK))

      do ll=1,N
         C_R(2*ll-1)=real(input(ll),kind=rk)
         C_R(2*ll)=aimag(input(ll))
      enddo

!
! --- INITIALIZE FFT
      !!CALL CFFT1I (N,WSAVE,LENSAV,IER)
      CALL CFFT1I (N,fft_sign,IER)
      IF (IER.NE.0) THEN
         WRITE(6,*) 'fast_fourier_transform_module%fft_complex_1D: error in routine CFFT1I'
         select case (IER)
         case(2); write(*,*) 'parameter LENSAV ',LENSAV,' not big enough'
         end select
         STOP 'error in fast_fourier_transform_module%fft_complex_1D'
      END IF
!
! --- PERFORM FORWARD TRANSFORM
      if(forward_backward_switch==1) then
         !!CALL CFFT1F (N,INC,C_R,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
         !!CALL CFFT1F (N,INC,C_R,LENC,WORK,LENWRK,fft_sign,IER)
         CALL CFFT1F (N,INC,C_R,LENC,fft_sign,IER)
      elseif(forward_backward_switch==-1) then
         !!CALL CFFT1B (N,INC,C_R,LENC,WSAVE,LENSAV,WORK,LENWRK,fft_sign,IER)
         !!CALL CFFT1B (N,INC,C_R,LENC,WORK,LENWRK,fft_sign,IER)
         CALL CFFT1B (N,INC,C_R,LENC,fft_sign,IER)
      else
         WRITE(*,'(A)') 'fast_fourier_transform_module%fft_complex_1D: error  '
         WRITE(*,'(A,i0)') 'forward_backward_switch=',forward_backward_switch
         WRITE(*,'(A,i0)') 'must be -1 or +1'
         stop 'error in fast_fourier_transform_module%fft_complex_1D '
      endif
      IF (IER.NE.0) THEN
         WRITE(6,*) 'fast_fourier_transform_module%fft_complex_1D: error in routine CFFT1(BF) !'
         select case (IER)
         case(1); write(*,*) 'parameter LENC ',LENC,' not big enough'
         case(2); write(*,*) 'parameter LENSAV ',LENSAV,' not big enough'
         case(3); write(*,*) 'parameter LENWRK ',LENWRK,' not big enough'
         end select
         STOP 'error in subroutine fast_fourier_transform_module%fft_complex_1D '
      END IF

      do ll=1,N
         output(ll)=cmplx(C_R(2*ll-1),C_R(2*ll))
      enddo
!
       WRITE(*,'(A,/)') ' end fast_fourier_transform_module%fft_complex_1D '
end subroutine fft_complex_1D
!************************************************************************************************************************
subroutine fft_real_1D(input,output,forward_backward_switch)
      integer                                              :: forward_backward_switch
      integer                                              :: INC,LENR,LENSAV,LENWRK
      integer                                              :: IER
      integer                                              :: N
      real(kind=rk),dimension(:)                           :: input
      real(kind=rk),dimension(:),allocatable,intent(out)   :: output
      real(kind=rk),dimension(:),allocatable               :: WSAVE,WORK  
      type(fft_real_sign_type) :: fft_sign
!
      INC=1 ! stride of data in input
      N=size(input)
      LENSAV= N + INT(LOG(REAL(N))/LOG(2.)) + 4
      LENWRK=N
      LENR=INC*(N-1) + 1
      !!allocate(output(N))
      allocate(WSAVE(LENSAV))
      allocate(WORK(LENWRK))

      output=input

!
! --- INITIALIZE FFT
      CALL RFFT1I (N,fft_sign, IER)
      IF (IER.NE.0) THEN
         WRITE(6,*) 'fast_fourier_transform_module%fft_real_1D: error in routine RFFT1I'
         select case (IER)
         case(2); write(*,*) 'parameter LENSAV ',LENSAV,' not big enough'
         end select
         STOP ' error in fast_fourier_transform_module%fft_real_1D'
      END IF
!
! --- PERFORM FORWARD TRANSFORM
      if(forward_backward_switch==1) then
         !!CALL RFFT1F (N,INC,output,LENR,WORK,LENWRK,fft_sign,IER)
         CALL RFFT1F (N,INC,output,LENR,fft_sign,IER)
      elseif(forward_backward_switch==-1) then
         !!CALL RFFT1B (N,INC,output,LENR,WORK,LENWRK,fft_sign,IER)
         CALL RFFT1B (N,INC,output,LENR,fft_sign,IER)
      else
         WRITE(*,'(A)') 'fast_fourier_transform_module%fft_real_1D: error in subroutine fft_real_1D '
         WRITE(*,'(A,i0)') 'forward_backward_switch=',forward_backward_switch
         WRITE(*,'(A,i0)') 'must be -1 or +1'
         stop 'error in subroutine fast_fourier_transform_module%fft_real_1D '
      endif
      IF (IER.NE.0) THEN
         WRITE(6,*) 'fast_fourier_transform_module%fft_real_1D: error in routine RFFT1F !'
         select case (IER)
         case(1); write(*,*) 'parameter LENR ',LENR,' not big enough'
         case(2); write(*,*) 'parameter LENSAV ',LENSAV,' not big enough'
         case(3); write(*,*) 'parameter LENWRK ',LENWRK,' not big enough'
         end select
         STOP
      END IF
!
       WRITE(*,'(A,/)') ' end fast_fourier_transform_module%fft_real_1D '
end subroutine fft_real_1D
!************************************************************************************************************************
end module fast_fourier_transform_module
