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
    INTEGER,DIMENSION(3)                                      ::GRAD_KSIZE
    REAL(KIND=RK8),ALLOCATABLE,DIMENSION(:,:)                 ::V_ARR
    TYPE(vector),ALLOCATABLE,DIMENSION(:,:,:)                 ::DGRAD_VEC
    INTEGER                                                   ::np,tot
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGRAD
    INTEGER                                                   ::ii,jj,kk,ll
    INTEGER                                                   ::i,j,k
    INTEGER                                                   ::a
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,PATCH
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:)                      ::PATCH3
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end
    REAL(Kind=rk8)                                            ::CV1,CV2,CV3,GV
    REAL(Kind=rk8)                                            ::min,max
    REAL(Kind=rk8)                                            ::ADD
    character(len=256)::str

WRITE(*,'(/,"*",80("="),"*",/)')
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* Adaptive QR FILTER begin ...'
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(3(A,I0))')'==> Kernel size :: ',KSIZE(1),'x',KSIZE(2),'x',KSIZE(3)

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

CTR=(KSIZE-1)/2

IF(KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
    CALL vandermonde_1d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 1D is being called'
    GRAD_KSIZE=(/3,1,1/)
ELSEIF(KSIZE(3) == 1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_2d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 2D is being called'
    GRAD_KSIZE=(/3,3,1/)
ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_3d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 3D is being called'
    GRAD_KSIZE=(/3,3,3/)
ELSE
    write(*,*)'CHECK the dimesions of the array being passed in MOD_QR.f90'
STOP
END IF

WRITE(*,'(3(A,I0))')'==> Kernel size for gradient :: ',GRAD_KSIZE(1),'x',GRAD_KSIZE(2),'x',GRAD_KSIZE(3)

        DGRAD=ADATA

        CALL GRADIENT(DGRAD,DGRAD_VEC,GRAD_KSIZE,NP)

        CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        min=minval(DGRAD)
        max=maxval(DGRAD)
        CV1=50.0
        CV2=25.0
        CV3=10.0
        !CV1=max/10.0
        !CV2=cv1/5.0
        !CV3=cv2/2.0

        write(*,*)'min of gradient ::',min
        write(*,*)'CV3             ::',CV3
        write(*,*)'CV2             ::',CV2
        write(*,*)'CV1             ::',CV1
        write(*,*)'max of gradient ::',max

        dim_mm=size(V_arr,2)
        dim_nn=size(V_arr,1)

WRITE(*,'(A,I0,A,I0)')'==> DIM of Vandermonde MATRIX(R,C)::',dim_nn,'x',dim_mm

        IF(allocated(cc))deallocate(cc)
        allocate(cc(1:dim_mm))
        IF(allocated(dd))deallocate(dd)
        allocate(dd(1:dim_mm))


OPEN(unit=300,file='V_arr.dat')
Do j=1,dim_nn
write(300,*)real(v_arr(j,1:dim_mm),4)
Enddo
write(300,*)'======================================='


CALL QR_dcmp_non_square(V_arr,dim_nn,dim_mm,cc,dd,sing)

Do j=1,dim_nn
write(300,*)real(v_arr(j,1:dim_mm),4)
Enddo
CLOSE(300)
        IF(allocated(ATEMP))deallocate(ATEMP)
        allocate(ATEMP(dim_nn))

        ATEMP=0.0
        ATEMP(1:dim_mm)=V_arr((dim_nn+1)/2,1:dim_mm)
        !ATEMP(dim_mm)=V_arr(8,6)
        !ATEMP(dim_mm-1)=V_arr(8,5)



!OPEN(unit=200,file='array.dat',position='append',status='old')
OPEN(unit=200,file='array.dat',position='append')
write(200,*)dim_nn,dim_mm
write(200,*)real(ATEMP(1:dim_mm),4),sum(ATEMP(1:dim_mm))
close(200)

        IF(allocated(PATCH))deallocate(PATCH)
        allocate(PATCH(1:dim_nn))


CALL QR_solv_non_square_FOR_PATCH(V_arr,dim_nn,dim_mm,cc,dd,ATEMP,PATCH)

write(str,'(A,I0,A)')'patch_',ksize(1),'.csv'
OPEN(unit=100,file=trim(str))
write(100,'(A)')'i,j,k,scalar'
kk=1
do k=-ctr(3),ctr(3)
do j=-ctr(2),ctr(2)
do i=-ctr(1),ctr(1)
write(100,'(3(I0,A),F10.6)')i,',',j,',',k,',',patch(kk)
kk=kk+1
enddo
enddo
enddo
close(100)

STOP
!write(*,*)'i,k'
!write(*,*)"i,k,'=',j,i,'*',j,k"
write(*,*)"j,k,'=',j,k,'-',j,i,'*',i,k"
DO k=1,6
DO i=1,k-1
DO j=1,9
write(*,*)j,k,'=',j,k,'-',j,i,'*',i,k
ENDDO
write(*,*)'--------------------------------'
ENDDO
ENDDO



DO ll=1,1
IF(ll == 1)GV=max
IF(ll == 2)GV=CV1
IF(ll == 3)GV=CV2
IF(ll == 4)GV=CV3

call CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,A,ADD,TOT)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length
!            write(*,'(5(A,I0))')"thread_num:",thread_num, &
!                        "; num_threads:",omp_get_num_threads(),&
!                        "; kk_start:",kk_start,&
!                        "; kk_end:", kk_end,&
!                        "; last should be:",DSIZE(3)
TOT=0
DO KK=KK_START,KK_END

!        if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
!            write(*,'(3(A,I0))') "thread_num:",thread_num, &
!            "; kk:",kk, " from ",kk_end
!        end if

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

write(*,'(A,F7.3,A,I0,A,I0,A,F6.2)')'==> &
                &number of pixels less then ',GV,' out of ',SIZE(ADATA),' :: ',TOT,'&
                &. In percentage :: ',TOT*100.0/SIZE(ADATA)

!$OMP END PARALLEL

DO kk=1,DSIZE(3)
DO jj=1,DSIZE(2)
DO ii=1,DSIZE(1)
DGHOST(ii,jj,kk)=ADATA(ii,jj,kk)
END DO
END DO
END DO

call CPU_TIME(cput2)

write(*,*)
write(*,'(A,F10.5)')'==> CPU TIME:',cput2-cput1
write(*,*)

ENDDO 

WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* Adaptive QR FILTER end ...'
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(/,"*",80("="),"*",/)')
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
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::ATEMP,PATCH
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end
    REAL(Kind=rk8)                                            ::ADD

WRITE(*,'(/,"*",80("="),"*",/)')
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* QR FILTER begin ...'
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(3(A,I0))')'==> Kernel size :: ',KSIZE(1),'x',KSIZE(2),'x',KSIZE(3)

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

CTR=(KSIZE-1)/2

IF(KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
    CALL vandermonde_1d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 1D is being called'

ELSEIF(KSIZE(3) == 1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_2d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 2D is being called'

ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    CALL vandermonde_3d(V_arr,KSIZE)
    WRITE(*,'(A)')'==> Vandermnode 3D is being called'

ELSE
    write(*,*)'CHECK the dimesions of the array being passed in MOD_QR.f90'
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
        allocate(ATEMP(dim_nn))

        ATEMP=0.0
        !ATEMP(1:dim_mm)=V_arr((dim_nn+1)/2,1:dim_mm)
        ATEMP(1)=1.0

        IF(allocated(PATCH))deallocate(PATCH)
        allocate(PATCH(1:dim_nn))

CALL QR_solv_non_square_FOR_PATCH(V_arr,dim_nn,dim_mm,cc,dd,ATEMP,PATCH)

!write(*,*)'patch is'
!write(*,*)patch
!kk=0
!do i=1,KSIZE(1)
!    write(*,'(5(F10.5))')patch(1+kk:ksize(1)+kk)
!    kk=kk+ksize(1)
!    enddo


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


!$OMP END PARALLEL
DO kk=1,DSIZE(3)
DO jj=1,DSIZE(2)
DO ii=1,DSIZE(1)
DGHOST(ii,jj,kk)=ADATA(ii,jj,kk)
END DO
END DO
END DO

call CPU_TIME(cput2)

write(*,*)
write(*,'(A,F10.5)')'==> CPU TIME:',cput2-cput1
write(*,*)

WRITE(*,'(10(" "),A)')'* QR FILTER end ...'
WRITE(*,'(/,"*",80("="),"*",/)')
END SUBROUTINE QR
!=============================================================================
END MODULE MOD_QR
