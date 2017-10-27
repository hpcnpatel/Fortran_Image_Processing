MODULE mod_gauss_filter

use constants_module
use variable_module
use OMP_LIB
use ghost

IMPLICIT NONE

CONTAINS
!=============================================================================
SUBROUTINE GAUSS(ADATA,KSIZE,np)
 
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    INTEGER                                                   ::KSHAPE
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    Real(kind=rk8),dimension(:,:,:),allocatable               ::G_ARR
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::DTEMP
    INTEGER                                                  ::thread_num
    INTEGER                                                  ::num_threads
    INTEGER                                                  ::local_length
    INTEGER                                                  ::gl_length
    INTEGER                                                  ::kk_start
    INTEGER                                                  ::kk_end
    INTEGER                                                   ::i,j,k

WRITE(*,*)'========================================================================'
WRITE(*,*)'Gauss Filter'
WRITE(*,*)'========================================================================'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',size(ADATA,1),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',size(ADATA,2),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',size(ADATA,3),' *****'

     KSHAPE=KSIZE(1)*KSIZE(2)*KSIZE(3)
     CTR=(KSIZE-1)/2

        CALL ghost_region(ADATA,DGHOST,KSIZE)

        CALL GAUSS_PATCH(G_ARR,KSIZE)

    IF(ALLOCATED(DTEMP))DEALLOCATE(DTEMP)
    ALLOCATE(DTEMP(KSHAPE))

CALL SYSTEM_CLOCK(t1)
CALL CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(DTEMP,THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,A)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length
            write(*,'(5(A,I0))')"thread_num:",thread_num, &
                        "; num_threads:",omp_get_num_threads(),&
                        "; kk_start:",kk_start,&
                        "; kk_end:", kk_end,&
                        "; last should be:",DSIZE(3)
!DO kk=1,DSIZE(3)
DO KK=KK_START,KK_END

        if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
            write(*,'(3(A,I0))') "thread_num:",thread_num, &
            "; kk:",kk, " from ",kk_end
        end if

Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

        a=0
        do k=-CTR(3),CTR(3)
        do j=-CTR(2),CTR(2)
        do i=-CTR(1),CTR(1)
!         do j=-CTR(2)+abs(k),CTR(2)-abs(k)
!         do i=-CTR(1)+abs(j+k),CTR(1)-abs(j+k)
            a=a+1
            DTEMP(a)=DGHOST(ii+i,jj+j,kk+k)*G_arr(i,j,k)
        enddo
        enddo
        enddo
        ADATA(ii,jj,kk)=SUM(DTEMP)!/a

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
subroutine gauss_patch(G_arr,KSIZE)

INTEGER,DIMENSION(3)                                      ::CTR
integer                                                   ::X,Y,Z,i
integer                                                   ::XX,YY,ZZ
real(kind=rk8),allocatable,dimension(:,:,:),INTENT(OUT)   ::G_arr
integer,dimension(3),INTENT(INOUT)                        ::KSIZE
integer,dimension(3)                                      ::len
real(kind=rk8),dimension(3)                               ::sigma
real(kind=rk8)                                            ::cc,rr
real(kind=rk8) ::pow
!=====================

LEN=KSIZE

CTR=(LEN-1)/2

DO i=1,3
if(len(i)==1)then
    sigma(i)=0.00
elseif(len(i)==3)then
    sigma(i)=0.25
elseif(len(i)==5)then
    sigma(i)=0.75
elseif(len(i)==7)then
    sigma(i)=1.0
elseif(len(i)==9)then
    sigma(i)=1.25
elseif(len(i)==11)then
    sigma(i)=1.75
elseif(len(i)==13)then
    sigma(i)=2.0
elseif(len(i)==15)then
    sigma(i)=2.25
elseif(len(i)==17)then
    sigma(i)=2.75
elseif(len(i)==19)then
    sigma(i)=3.0
elseif(len(i)==21)then
    sigma(i)=3.25
elseif(len(i)==23)then
    sigma(i)=3.75
elseif(len(i)==25)then
    sigma(i)=4.0
else
    WRITE(*,*)'VALUE OF PATCH SIZE IS OUT OF RANGE, CHOOSE BETWEEN 1 AND 25'
    STOP
endif
END DO

        allocate(G_arr(-CTR(1):CTR(1),&
                       -CTR(2):CTR(2),&
                       -CTR(3):CTR(3)))
G_arr=0.0
        IF(KSIZE(3) == 1) pow=1.0_8 ! if 2D
        IF(KSIZE(2) == 1 .and. KSIZE(3)==1) pow=0.5_8 !if 1D
        IF(KSIZE(1) > 1 .and. KSIZE(2)>1 .and. KSIZE(3) > 1) pow=3_8/2_8 !if 3D

        cc=(1/((2*pi)**(pow))*sigma(1)**2)

write(*,*)'Gaussian deviation is::',sigma(1:3)
write(*,*)'cc is::',cc,'CTR is',CTR

!Do Z=-CTR(3),CTR(3)
Do Y=-CTR(2),CTR(2)
Do X=-CTR(1),CTR(1)

        rr=(X**2/2*sigma(1)**2)+(Y**2/2*sigma(2)**2)!+(Z**2/2*sigma(3)**2)
!write(*,*)X,Y,Z,'rr',rr,X**2/2*sigma(1)**2,Y**2/2*sigma(2)**2,Z**2/2*sigma(3)**2
        G_arr(X,Y,0)= cc*exp(-rr)  

END DO
write(*,*)G_arr(:,Y,0)
END DO
!END DO

G_arr=G_arr/sum(G_arr)

end subroutine gauss_patch
!=============================================================================
end module
