module mod_grad_of_image

USE OMP_LIB
USE VARIABLE_MODULE
USE CONSTANTS_MODULE
USE VANDERMONDE_MODULE
USE QR_MODULE_PATCH
USE SORT_MODULE
USE GHOST
USE MOD_GNUPLOT

IMPLICIT NONE

CONTAINS
!=============================================================================
SUBROUTINE GRADIENT(ADATA,GRADATA,KSIZE,NP)

    Logical                                                   ::sing
    REAL(kind=rk8),dimension(:),allocatable                   ::cc,dd
    INTEGER                                                   ::dim_nn
    INTEGER                                                   ::dim_mm
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    REAL(KIND=RK8),ALLOCATABLE,DIMENSION(:,:)                 ::V_ARR
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    TYPE(vector),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)     ::GRADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::i,j,k
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,ALPHA
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end

!WRITE(*,'(/,"*",80("="),"*",/)')
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* calculating grandient ...'
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(3(A,I0))')'==> Kernel size for gradient:: ',KSIZE(1),'x',KSIZE(2),'x',KSIZE(3)

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

CTR=(KSIZE-1)/2

IF(KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
    CALL vandermonde_1d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 1D is being called for grad calculation'
ELSEIF(KSIZE(3) ==  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_2d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 2D is being called for grad calculation'
ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_3d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 3D is being called for grad calculation'
ELSE
    write(*,*)'CHECK the dimesions of the array being passed in MOD_GRAD_OF_IMAGE.f90'
STOP
END IF

CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        dim_mm=size(V_arr,2)
        dim_nn=size(V_arr,1)

    write(*,'(A,I0,A,I0)')'==> DIM of Vandermonde MATRIX(R,C)::',dim_nn,' X ',dim_mm

        IF(allocated(cc))deallocate(cc)
        allocate(cc(1:dim_mm))
        IF(allocated(dd))deallocate(dd)
        allocate(dd(1:dim_mm))

CALL QR_dcmp_non_square(V_arr,dim_nn,dim_mm,cc,dd,sing)

        IF(allocated(ATEMP))deallocate(ATEMP)
        allocate(ATEMP(1:dim_nn))

        IF(allocated(ALPHA))deallocate(ALPHA)
        allocate(ALPHA(1:dim_mm))

call CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

ALLOCATE(GRADATA(DSIZE(1),DSIZE(2),DSIZE(3)))

!$OMP PARALLEL private(THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,A,ATEMP,ALPHA)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length
!            write(*,'(5(A,I0))')"thread_num:",thread_num, &
!                        "; num_threads:",omp_get_num_threads(),&
!                        "; kk_start:",kk_start,&
!                        "; kk_end:", kk_end,&
!                        "; last should be:",DSIZE(3)
DO KK=KK_START,KK_END

!        if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
!            write(*,'(3(A,I0))') "thread_num:",thread_num, &
!            "; kk:",kk, " from ",kk_end
!        end if

!Do kk=1,DSIZE(3)
Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

A=0
        DO K=-CTR(3),CTR(3)
        DO J=-CTR(2),CTR(2)
        DO I=-CTR(1),CTR(1)
!        DO J=-CTR(2)+ABS(K),CTR(2)-ABS(K)
!        DO I=-CTR(1)+ABS(J+K),CTR(1)-ABS(J+K)
            A=A+1
            !ADD=ADD+PATCH(A)*DGHOST(II+I,JJ+J,KK+K)
            ATEMP(A)=DGHOST(ii+i,jj+j,kk+k)
        END DO
        END DO
        END DO

CALL QR_solv_non_square(V_arr,dim_nn,dim_mm,cc,dd,ATEMP,ALPHA)

!write(*,*)DGHOST(ii,jj,kk),ALPHA

IF(    KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
GRADATA(II,JJ,KK)%x=ALPHA(2)
GRADATA(II,JJ,KK)%y=0.0D0
GRADATA(II,JJ,KK)%z=0.0D0
ELSEIF(KSIZE(3) == 1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
GRADATA(II,JJ,KK)%x=ALPHA(3)
GRADATA(II,JJ,KK)%y=ALPHA(2)
GRADATA(II,JJ,KK)%z=0.0D0
ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
GRADATA(II,JJ,KK)%x=ALPHA(4)
GRADATA(II,JJ,KK)%y=ALPHA(3)
GRADATA(II,JJ,KK)%z=ALPHA(2)
ELSE
    write(*,*)'CHECK the dimesions of the array in mod_grad_of_image.f90'
STOP
END IF
!GRADATA(II,JJ,KK)%x=ALPHA(4)
!GRADATA(II,JJ,KK)%y=ALPHA(3)
!GRADATA(II,JJ,KK)%z=ALPHA(2)
ADATA(II,JJ,KK)=sqrt((GRADATA(II,JJ,KK)%x)**2+&
                    &(GRADATA(II,JJ,KK)%y)**2+&
                    &(GRADATA(II,JJ,KK)%z)**2)
END DO
END DO
END DO

!$OMP END PARALLEL

call CPU_TIME(cput2)

write(*,*)
write(*,'(A,F10.5)')'==> CPU TIME:',cput2-cput1
write(*,*)

WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* gradient calculated ...'
WRITE(*,'(10(" "),35("*"))')
END SUBROUTINE GRADIENT
!=============================================================================
SUBROUTINE GRADIENT_Central_difference(ADATA,GRADATA,KSIZE,NP)

    Logical                                                   ::sing
    REAL(kind=rk8),dimension(:),allocatable                   ::cc,dd
    INTEGER                                                   ::dim_nn
    INTEGER                                                   ::dim_mm
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    REAL(KIND=RK8),ALLOCATABLE,DIMENSION(:,:)                 ::V_ARR
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    TYPE(vector),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)     ::GRADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::i,j,k
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,ALPHA
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end

WRITE(*,*)'*********************************************************'
write(*,*)'GRADIENT OF AN IMAGE WITH CENTRAL DIFFERENCE'
WRITE(*,*)'*********************************************************'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)
CTR=(KSIZE-1)/2
KSIZE(1:2)=5
ALLOCATE(GRADATA(DSIZE(1),DSIZE(2),DSIZE(3)))

CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        WRITE(*,'(A,I0,A)')

call system_clock(t1)
call CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(THREAD_NUM,KK_START,KK_END,II,JJ,KK)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length

DO KK=KK_START,KK_END

!Do kk=1,DSIZE(3)
Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

!Central difference of order O(h^2)
!GRADATA(II,JJ,KK)%x=(DGHOST(II-1,JJ,KK)-DGHOST(II+1,JJ,KK))/2.0D0
!GRADATA(II,JJ,KK)%y=(DGHOST(II,JJ-1,KK)-DGHOST(II,JJ+1,KK))/2.0D0
!Central difference of order O(h^4)
GRADATA(II,JJ,KK)%x=(-DGHOST(II-2,JJ,KK)&
                    &+8*DGHOST(II+1,JJ,KK)&
                    &-8*DGHOST(II-1,JJ,KK)&
                    &+DGHOST(II+2,JJ,KK))/12.0D0
GRADATA(II,JJ,KK)%y=(-DGHOST(II,JJ-2,KK)&
                    &+8*DGHOST(II,JJ+1,KK)&
                    &-8*DGHOST(II,JJ-1,KK)&
                    &+DGHOST(II,JJ+2,KK))/12.0D0
if(DSIZE(3) > 1)GRADATA(II,JJ,KK)%z=(DGHOST(II,JJ,KK-1)-DGHOST(II,JJ,KK+1))/2.0D0
if(DSIZE(3) <= 1)GRADATA(II,JJ,KK)%z=0.0
ADATA(II,JJ,KK)=sqrt(((GRADATA(ii,jj,kk)%x)**2)+&
                    &((GRADATA(ii,jj,kk)%y)**2)+&
                    &((GRADATA(ii,jj,kk)%z)**2))
END DO
END DO
END DO

!$OMP END PARALLEL

call system_clock(t2)
call CPU_TIME(cput2)

write(*,*)'====================='
write(*,*)'SYSTEM CLOCK TIME:',t2-t1
write(*,*)'====================='
write(*,*)'====================='
write(*,*)'CPU TIME:',cput2-cput1
write(*,*)'====================='

END SUBROUTINE
!=============================================================================
SUBROUTINE DIV(ADATA,GRADATA,KSIZE,NP)

    Logical                                                   ::sing
    REAL(kind=rk8),dimension(:),allocatable                   ::cc,dd
    INTEGER                                                   ::dim_nn
    INTEGER                                                   ::dim_mm
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    REAL(KIND=RK8),ALLOCATABLE,DIMENSION(:,:)                 ::V_ARR
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    TYPE(vector),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)     ::GRADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::i,j,k
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,ALPHA
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end

WRITE(*,*)'*********************************************************'
write(*,*)'GRADIENT OF AN IMAGE'
WRITE(*,*)'*********************************************************'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)
CTR=(KSIZE-1)/2

ALLOCATE(GRADATA(DSIZE(1),DSIZE(2),DSIZE(3)))

CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        WRITE(*,'(A,I0,A)')

IF(    KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
    CALL vandermonde_1d(V_arr,KSIZE)
    write(*,*)'It is 1D data so Vandermnode 1D is being called'
ELSEIF(KSIZE(3) == 1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_2d(V_arr,KSIZE)
    write(*,*)'It is 2D data so Vandermnode 2D is being called'
ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_3d(V_arr,KSIZE)
    write(*,*)'It is 3D data so Vandermnode 3D is being called'
ELSE
    write(*,*)'CHECK the dimesions of the array being passed in MOD_QR.f90'
STOP
END IF

        dim_mm=size(V_arr,2)
        dim_nn=size(V_arr,1)

    write(*,'(A,I0,A,I0)')'DIM of Vandermonde MATRIX(R,C)::',dim_nn,' X ',dim_mm

        IF(allocated(cc))deallocate(cc)
        allocate(cc(1:dim_mm))
        IF(allocated(dd))deallocate(dd)
        allocate(dd(1:dim_mm))

CALL QR_dcmp_non_square(V_arr,dim_nn,dim_mm,cc,dd,sing)

        IF(allocated(ATEMP))deallocate(ATEMP)
        allocate(ATEMP(1:dim_nn))

        IF(allocated(ALPHA))deallocate(ALPHA)
        allocate(ALPHA(1:dim_mm))

call system_clock(t1)
call CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,A,ATEMP,ALPHA)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length
!            write(*,'(5(A,I0))')"thread_num:",thread_num, &
!                        "; num_threads:",omp_get_num_threads(),&
!                        "; kk_start:",kk_start,&
!                        "; kk_end:", kk_end,&
!                        "; last should be:",DSIZE(3)
DO KK=KK_START,KK_END

!        if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
!            write(*,'(3(A,I0))') "thread_num:",thread_num, &
!            "; kk:",kk, " from ",kk_end
!        end if

!Do kk=1,DSIZE(3)
Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

A=0
        DO K=-CTR(3),CTR(3)
        DO J=-CTR(2),CTR(2)
        DO I=-CTR(1),CTR(1)
!        DO J=-CTR(2)+ABS(K),CTR(2)-ABS(K)
!        DO I=-CTR(1)+ABS(J+K),CTR(1)-ABS(J+K)
            A=A+1
            !ADD=ADD+PATCH(A)*DGHOST(II+I,JJ+J,KK+K)
            ATEMP(A)=DGHOST(ii+i,jj+j,kk+k)
        END DO
        END DO
        END DO

CALL QR_solv_non_square(V_arr,dim_nn,dim_mm,cc,dd,ATEMP,ALPHA)

IF(    KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
GRADATA(II,JJ,KK)%x=ALPHA(3)
GRADATA(II,JJ,KK)%y=0.0D0
GRADATA(II,JJ,KK)%z=0.0D0
ELSEIF(KSIZE(3) == 1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
GRADATA(II,JJ,KK)%x=ALPHA(6)
GRADATA(II,JJ,KK)%y=ALPHA(5)
GRADATA(II,JJ,KK)%z=0.0D0
ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
GRADATA(II,JJ,KK)%x=ALPHA(4)
GRADATA(II,JJ,KK)%y=ALPHA(3)
GRADATA(II,JJ,KK)%z=ALPHA(2)
ELSE
    write(*,*)'CHECK the dimesions of the array in mod_grad_of_image.f90'
STOP
END IF
!GRADATA(II,JJ,KK)%x=ALPHA(4)
!GRADATA(II,JJ,KK)%y=ALPHA(3)
!GRADATA(II,JJ,KK)%z=ALPHA(2)
ADATA(II,JJ,KK)=(GRADATA(II,JJ,KK)%x)+&
               &(GRADATA(II,JJ,KK)%y)+&
               &(GRADATA(II,JJ,KK)%z)
END DO
END DO
END DO

!$OMP END PARALLEL

call system_clock(t2)
call CPU_TIME(cput2)

write(*,*)'====================='
write(*,*)'SYSTEM CLOCK TIME:',t2-t1
write(*,*)'====================='
write(*,*)'====================='
write(*,*)'CPU TIME:',cput2-cput1
write(*,*)'====================='

END SUBROUTINE DIV
!=============================================================================
END MODULE 
