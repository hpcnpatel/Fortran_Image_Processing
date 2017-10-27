MODULE MOD_QR

    USE OMP_LIB
    USE GHOST
    USE VARIABLE_MODULE
    USE CONSTANTS_MODULE
    USE VANDERMONDE_MODULE
    USE QR_MODULE_PATCH
    !USE QR_MODULE
    USE mod_grad_of_image

IMPLICIT NONE

CONTAINS
!=============================================================================
SUBROUTINE ADAPTIVE_QR(ADATA,KSIZE,NP)

    Logical                                                   ::sing
    REAL(kind=rk8),dimension(:),allocatable                   ::cc,dd
    INTEGER                                                   ::dim_nn
    INTEGER                                                   ::dim_mm
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    REAL(KIND=RK8),ALLOCATABLE,DIMENSION(:,:)                 ::V_ARR
    TYPE(vector),ALLOCATABLE,DIMENSION(:,:,:)                 ::DGRAD_VEC
    INTEGER                                                   ::np,tot
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGRAD
    INTEGER                                                   ::ii,jj,kk,ll
    INTEGER                                                   ::i,j,k
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,PATCH
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end
    REAL(Kind=rk8)                                            ::CV1,CV2,CV3,GV
    REAL(Kind=rk8)                                            ::min,max
    REAL(Kind=rk8)                                            ::ADD

WRITE(*,*)'========================================================================'
write(*,*)'Adaptive QR FILTER'
WRITE(*,*)'========================================================================'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

CTR=(KSIZE-1)/2


        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',size(ADATA,1),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',size(ADATA,2),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',size(ADATA,3),' *****'
        WRITE(*,'(A,3I0,A)')'***** Kernel size **********: ',KSIZE(1:3),' *****'

DGRAD=ADATA
CALL GRADIENT(DGRAD,DGRAD_VEC,KSIZE,NP)
CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        min=minval(DGRAD)
        max=maxval(DGRAD)
        CV1=50.0
        CV2=25.0
        CV3=10.0
        !CV1=max/10.0
        !CV2=cv1/5.0
        !CV3=cv2/2.0

        write(*,*)'min ::',min
        write(*,*)'CV3 ::',CV3
        write(*,*)'CV2 ::',CV2
        write(*,*)'CV1 ::',CV1
        write(*,*)'max ::',max


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
        allocate(ATEMP(dim_nn))

        ATEMP=0.0
        ATEMP(1:dim_mm)=V_arr((dim_nn+1)/2,1:dim_mm)

        IF(allocated(PATCH))deallocate(PATCH)
        allocate(PATCH(1:dim_nn))

CALL QR_solv_non_square_FOR_PATCH(V_arr,dim_nn,dim_mm,cc,dd,ATEMP,PATCH)

!write(*,*)'patch is'
!kk=0
!do k=1,DSIZE(3)
!do j=1,DSIZE(2)
!    write(*,*)patch(1+kk:ksize(1)+kk)
!    kk=kk+ksize(1)
!    enddo
!    enddo
!===========================================================

DO ll=1,4
IF(ll == 1)GV=max
IF(ll == 2)GV=CV1
IF(ll == 3)GV=CV2
IF(ll == 4)GV=CV3

call cpu_time(cput1)
call system_clock(t1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(TOT,THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,A,ADD)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length
            write(*,'(5(A,I0))')"thread_num:",thread_num, &
                        "; num_threads:",omp_get_num_threads(),&
                        "; kk_start:",kk_start,&
                        "; kk_end:", kk_end,&
                        "; last should be:",DSIZE(3)
TOT=0
DO KK=KK_START,KK_END

        if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
            write(*,'(3(A,I0))') "thread_num:",thread_num, &
            "; kk:",kk, " from ",kk_end
        end if

!Do kk=1,DSIZE(3)
Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

    IF(DGRAD(ii,jj,kk) .lt. GV)then
    TOT=TOT+1

ADD=0.0
A=0
        DO K=-CTR(3),CTR(3)
        DO J=-CTR(2),CTR(2)
        DO I=-CTR(1),CTR(1)
!        DO J=-CTR(2)+ABS(K),CTR(2)-ABS(K)
!        DO I=-CTR(1)+ABS(J+K),CTR(1)-ABS(J+K)
            A=A+1
            ADD=ADD+PATCH(A)*DGHOST(II+I,JJ+J,KK+K)
        END DO
        END DO
        END DO
ADATA(II,JJ,KK)=ADD!+((1-CT_GRAD(II,JJ,KK)/MAX)*(ADD-CT_DATA(II,JJ,KK)))
    ENDIF

END DO
END DO
END DO

write(*,'(A,F7.3,A,I0,A,I0,A,F5.2)')'Scalar value of pixels less then ', GV ,&
&' out of ',SIZE(ADATA),' are ',TOT,' in percentage :: ',TOT*100.0/SIZE(ADATA)
write(*,*)'======================================================================='
!$OMP END PARALLEL

DO kk=1,DSIZE(3)
DO jj=1,DSIZE(2)
DO ii=1,DSIZE(1)
DGHOST(ii,jj,kk)=ADATA(ii,jj,kk)
END DO
END DO
END DO


call system_clock(t2)
call CPU_TIME(cput2)

write(*,*)'====================='
write(*,*)'SYSTEM CLOCK TIME:',t2-t1
write(*,*)'====================='
write(*,*)'====================='
write(*,*)'CPU TIME:',cput2-cput1
write(*,*)'====================='

ENDdo 

END SUBROUTINE ADAPTIVE_QR
!=============================================================================
SUBROUTINE QR(ADATA,KSIZE,NP)

    Logical                                                   ::sing
    REAL(kind=rk8),dimension(:),allocatable                   ::cc,dd
    INTEGER                                                   ::dim_nn
    INTEGER                                                   ::dim_mm
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    REAL(KIND=RK8),ALLOCATABLE,DIMENSION(:,:)                 ::V_ARR
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::i,j,k
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,PATCH
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end
    REAL(Kind=rk8)                                            ::ADD

WRITE(*,*)'========================================================================'
write(*,*)'QR FILTER'
WRITE(*,*)'========================================================================'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

CTR=(KSIZE-1)/2

CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',size(ADATA,1),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',size(ADATA,2),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',size(ADATA,3),' *****'
        WRITE(*,'(A,3I0,A)')'***** Kernel size **********: ',KSIZE(1:3),' *****'

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
        allocate(ATEMP(dim_nn))

        ATEMP=0.0
        ATEMP(1:dim_mm)=V_arr((dim_nn+1)/2,1:dim_mm)

        IF(allocated(PATCH))deallocate(PATCH)
        allocate(PATCH(1:dim_nn))

CALL QR_solv_non_square_FOR_PATCH(V_arr,dim_nn,dim_mm,cc,dd,ATEMP,PATCH)

call system_clock(t1)
call CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,A,ADD)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length
            write(*,'(5(A,I0))')"thread_num:",thread_num, &
                        "; num_threads:",omp_get_num_threads(),&
                        "; kk_start:",kk_start,&
                        "; kk_end:", kk_end,&
                        "; last should be:",DSIZE(3)
DO KK=KK_START,KK_END

        if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
            write(*,'(3(A,I0))') "thread_num:",thread_num, &
            "; kk:",kk, " from ",kk_end
        end if

!Do kk=1,DSIZE(3)
Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

ADD=0.0
A=0
        DO K=-CTR(3),CTR(3)
        DO J=-CTR(2),CTR(2)
        DO I=-CTR(1),CTR(1)
!        DO J=-CTR(2)+ABS(K),CTR(2)-ABS(K)
!        DO I=-CTR(1)+ABS(J+K),CTR(1)-ABS(J+K)
            A=A+1
            ADD=ADD+PATCH(A)*DGHOST(II+I,JJ+J,KK+K)
        END DO
        END DO
        END DO
ADATA(II,JJ,KK)=ADD!+((1-CT_GRAD(II,JJ,KK)/MAX)*(ADD-CT_DATA(II,JJ,KK)))

END DO
END DO
END DO

DO kk=1,DSIZE(3)
DO jj=1,DSIZE(2)
DO ii=1,DSIZE(1)
DGHOST(ii,jj,kk)=ADATA(ii,jj,kk)
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

END SUBROUTINE QR
!=============================================================================
END MODULE MOD_QR
