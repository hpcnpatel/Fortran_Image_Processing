MODULE mod_median

    USE mod_grad_of_image
    USE ghost
    USE constants_module
    USE variable_module
    USE OMP_LIB
    USE sort_module
    USE mod_gnuplot


IMPLICIT NONE

CONTAINS
!=============================================================================
SUBROUTINE MEDIAN(ADATA,KSIZE,np)

    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    INTEGER                                                   ::KSHAPE
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    INTEGER                                                   ::ct_pt
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::DTEMP
    INTEGER,ALLOCATABLE,DIMENSION(:)                          ::SORT_INDEX
    INTEGER                                                  ::thread_num
    INTEGER                                                  ::num_threads
    INTEGER                                                  ::local_length
    INTEGER                                                  ::gl_length
    INTEGER                                                  ::kk_start
    INTEGER                                                  ::kk_end
    INTEGER                                                   ::i,j,k

!WRITE(*,*)'========================================================================'
!WRITE(*,*)'MEDIAN module.....'
!WRITE(*,*)'========================================================================'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',size(ADATA,1),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',size(ADATA,2),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',size(ADATA,3),' *****'
        WRITE(*,'(A,3I0,A)')'***** Kernel size **********: ',KSIZE(1:3),' *****'

        CALL ghost_region(ADATA,DGHOST,KSIZE)

     KSHAPE=KSIZE(1)*KSIZE(2)*KSIZE(3)
     CTR=(KSIZE-1)/2

    IF(ALLOCATED(DTEMP))DEALLOCATE(DTEMP)
    ALLOCATE(DTEMP(KSHAPE))
    IF(ALLOCATED(SORT_INDEX))DEALLOCATE(SORT_INDEX)
    ALLOCATE(SORT_INDEX(KSHAPE))

CALL SYSTEM_CLOCK(t1)
CALL CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(DTEMP,SORT_INDEX,THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,CT_PT,A)

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
            DTEMP(a)=DGHOST(ii+i,jj+j,kk+k)
        enddo
        enddo
        enddo
        call SORTRX_real8(DTEMP(1:a),SORT_INDEX)
        CT_PT=(a+1)/2
        !ADATA(ii,jj,kk)=ADATA(ii,jj,kk)+&
        !&((1-(DGRAD(ii,jj,kk)/max))*(DTEMP(SORT_INDEX(CT_PT))-ADATA(ii,jj,kk)))
        ADATA(ii,jj,kk)=DTEMP(SORT_INDEX(CT_PT))

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
SUBROUTINE ADAPTIVE_MEDIAN(ADATA,KSIZE,np)

    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    INTEGER,DIMENSION(3)                                      ::GRAD_KSIZE
    INTEGER                                                   ::KSHAPE
    INTEGER                                                   ::np,tot
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    TYPE(vector),ALLOCATABLE,DIMENSION(:,:,:)                 ::DGRAD_VEC
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGRAD
    INTEGER                                                   ::ii,jj,kk,ll
    INTEGER                                                   ::a
    REAL                                                      ::cput1,cput2
    INTEGER                                                   ::ct_pt
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::DTEMP
    INTEGER,ALLOCATABLE,DIMENSION(:)                          ::SORT_INDEX
    INTEGER                                                   ::thread_num
    INTEGER                                                   ::num_threads
    INTEGER                                                   ::local_length
    INTEGER                                                   ::gl_length
    INTEGER                                                   ::kk_start
    INTEGER                                                   ::kk_end
    REAL(Kind=rk8) ::CV1,CV2,CV3,CV4,CV5
    REAL(Kind=rk8)                                            ::min,max
    INTEGER                                                   ::i,j,k
    REAL(Kind=rk8)                                            ::GV!gauge value

WRITE(*,'(/,"*",80("="),"*",/)')
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* Adaptive Median FILTER begin ...'
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(3(A,I0))')'==> Kernel size :: ',KSIZE(1),'x',KSIZE(2),'x',KSIZE(3)

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)


IF(KSIZE(3) == 1 .AND. KSIZE(2) == 1 .AND. KSIZE(1) > 1)then
    GRAD_KSIZE=(/3,1,1/)
ELSEIF(KSIZE(3) == 1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    GRAD_KSIZE=(/3,3,1/)
ELSEIF(KSIZE(3) >  1 .AND. KSIZE(2) > 1  .AND. KSIZE(1) > 1)then
    GRAD_KSIZE=(/3,3,3/)
ELSE
    write(*,*)'CHECK the dimesions of the array being passed in MOD_MEDIAN.f90'
STOP
END IF

WRITE(*,'(3(A,I0))')'==> Kernel size for gradient :: ',GRAD_KSIZE(1),'x',GRAD_KSIZE(2),'x',GRAD_KSIZE(3)

        DGRAD=ADATA

        CALL GRADIENT(DGRAD,DGRAD_VEC,GRAD_KSIZE,NP)

        CALL GHOST_REGION(ADATA,DGHOST,KSIZE)

        min=minval(DGRAD)
        max=maxval(DGRAD)

        CV1=100.0
        CV2=50.0
        CV3=10.0
        !CV1=50.0
        !CV2=25.0
        !CV3=10.0

        write(*,*)'min of gradient ::',min
        write(*,*)'CV3             ::',CV3
        write(*,*)'CV2             ::',CV2
        write(*,*)'CV1             ::',CV1
        write(*,*)'max of gradient ::',max

     KSHAPE=KSIZE(1)*KSIZE(2)*KSIZE(3)
     CTR=(KSIZE-1)/2

    IF(ALLOCATED(DTEMP))DEALLOCATE(DTEMP)
    ALLOCATE(DTEMP(KSHAPE))
    IF(ALLOCATED(SORT_INDEX))DEALLOCATE(SORT_INDEX)
    ALLOCATE(SORT_INDEX(KSHAPE))

DO ll=1,4
IF(ll == 1)GV=max
IF(ll == 2)GV=CV1
IF(ll == 3)GV=CV2
IF(ll == 4)GV=CV3

CALL CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(TOT,DTEMP,SORT_INDEX,THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,CT_PT,A)

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

!       if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
!           write(*,'(3(A,I0))') "thread_num:",thread_num, &
!           "; kk:",kk, " from ",kk_end
!       end if

Do jj=1,DSIZE(2)
Do ii=1,DSIZE(1)

    !IF(DGRAD(ii,jj,kk) .gt. CV1)write(*,*)DGRAD_VEC(ii,jj,kk)
    IF(DGRAD(ii,jj,kk) .lt. GV)then

        TOT=TOT+1
!goto 100
        a=0
         do k=-CTR(3),CTR(3)
         do j=-CTR(2),CTR(2)
         do i=-CTR(1),CTR(1)
!        do j=-CTR(2)+abs(k),CTR(2)-abs(k)
!        do i=-CTR(1)+abs(j+k),CTR(1)-abs(j+k)
            a=a+1
            DTEMP(a)=DGHOST(ii+i,jj+j,kk+k)
        enddo
        enddo
        enddo
        call SORTRX_real8(DTEMP(1:a),SORT_INDEX)
        CT_PT=(a+1)/2
        !ADATA(ii,jj,kk)=ADATA(ii,jj,kk)+&
        !&((1-(DGRAD(ii,jj,kk)/max))*(DTEMP(SORT_INDEX(CT_PT))-ADATA(ii,jj,kk)))
        ADATA(ii,jj,kk)=DTEMP(SORT_INDEX(CT_PT))
!100 continue
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

END DO

WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A)')'* Adaptive Median FILTER end ...'
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(/,"*",80("="),"*",/)')

END SUBROUTINE
!=============================================================================
SUBROUTINE MEAN(ADATA,KSIZE,np)
 
    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    INTEGER                                                   ::KSHAPE
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
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
WRITE(*,*)'MEAN FILTER'
WRITE(*,*)'========================================================================'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',size(ADATA,1),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',size(ADATA,2),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',size(ADATA,3),' *****'
        WRITE(*,'(A,3I0,A)')'***** Kernel size **********: ',KSIZE(1:3),' *****'

        CALL ghost_region(ADATA,DGHOST,KSIZE)

     KSHAPE=KSIZE(1)*KSIZE(2)*KSIZE(3)
     CTR=(KSIZE-1)/2

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
            DTEMP(a)=DGHOST(ii+i,jj+j,kk+k)
        enddo
        enddo
        enddo
        ADATA(ii,jj,kk)=SUM(DTEMP)/a

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
SUBROUTINE TRY_MEDIAN(ADATA,KSIZE,np)

    INTEGER,DIMENSION(3)                                      ::KSIZE,CTR,DSIZE
    INTEGER                                                   ::KSHAPE
    INTEGER                                                   ::np
    REAL(rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)      ::ADATA
    Real(kind=rk8),dimension(:,:,:),allocatable               ::DGHOST
    INTEGER                                                   ::ii,jj,kk
    INTEGER                                                   ::a
    INTEGER                                                   ::t1,t2
    REAL                                                      ::cput1,cput2
    INTEGER                                                   ::ct_pt
    REAL(rk8),ALLOCATABLE,DIMENSION(:)                        ::DTEMP
    INTEGER,ALLOCATABLE,DIMENSION(:)                          ::SORT_INDEX
    INTEGER                                                  ::thread_num
    INTEGER                                                  ::num_threads
    INTEGER                                                  ::local_length
    INTEGER                                                  ::gl_length
    INTEGER                                                  ::kk_start
    INTEGER                                                  ::kk_end
    INTEGER                                                   ::i,j,k

!WRITE(*,*)'========================================================================'
!WRITE(*,*)'try MEDIAN module.....'
!WRITE(*,*)'========================================================================'

DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',size(ADATA,1),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',size(ADATA,2),' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',size(ADATA,3),' *****'
        WRITE(*,'(A,3I0,A)')'***** Kernel size **********: ',KSIZE(1:3),' *****'

        CALL ghost_region(ADATA,DGHOST,KSIZE)

     KSHAPE=KSIZE(1)*KSIZE(2)*KSIZE(3)
     CTR=(KSIZE-1)/2

    IF(ALLOCATED(DTEMP))DEALLOCATE(DTEMP)
    ALLOCATE(DTEMP(KSHAPE))
    IF(ALLOCATED(SORT_INDEX))DEALLOCATE(SORT_INDEX)
    ALLOCATE(SORT_INDEX(KSHAPE))

CALL SYSTEM_CLOCK(t1)
CALL CPU_TIME(cput1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = DSIZE(3)
        local_length = gl_length/num_threads

!$OMP PARALLEL private(DTEMP,SORT_INDEX,THREAD_NUM,KK_START,KK_END,II,JJ,I,J,K,KK,CT_PT,A)

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

            a=a+1
            DTEMP(a)=DGHOST(ii+i,jj+j,kk+k)

        enddo
        enddo
        enddo
        call SORTRX_real8(DTEMP(1:a),SORT_INDEX)
        CT_PT=(a+1)/2
        !ADATA(ii,jj,kk)=ADATA(ii,jj,kk)+&
        !&((1-(DGRAD(ii,jj,kk)/max))*(DTEMP(SORT_INDEX(CT_PT))-ADATA(ii,jj,kk)))
        ADATA(ii,jj,kk)=DTEMP(SORT_INDEX(CT_PT))

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
END MODULE
