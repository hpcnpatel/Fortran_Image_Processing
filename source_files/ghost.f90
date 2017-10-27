Module ghost

      USE constants_module
      USE variable_module

      IMPLICIT NONE
CONTAINS
!=============================================================================
! To include ghost points outside the border of image data. 
!=============================================================================
SUBROUTINE ghost_region(ADATA,DGHOST,KSIZE)

INTEGER,dimension(3)                                      ::KSIZE,CTR,DSIZE
REAL(kind=rk8),dimension(:,:,:),ALLOCATABLE,INTENT(INOUT) ::ADATA
REAL(kind=rk8),dimension(:,:,:),ALLOCATABLE,INTENT(OUT)   ::DGHOST
INTEGER                                                   ::i,j,k
INTEGER                                                   ::ii,jj,kk


DSIZE(1)=size(ADATA,1)
DSIZE(2)=size(ADATA,2)
DSIZE(3)=size(ADATA,3)

CTR=(KSIZE-1)/2

        ALLOCATE(DGHOST(1-CTR(1):DSIZE(1)+CTR(1),&
                      &1-CTR(2):DSIZE(2)+CTR(2),&
                      &1-CTR(3):DSIZE(3)+CTR(3)))

DO kk=1,DSIZE(3)
DO jj=1,DSIZE(2)
DO ii=1,DSIZE(1)

DGHOST(ii,jj,kk)=ADATA(ii,jj,kk)

IF(CTR(1)==0)exit

    IF(ii==1)then
        DO i=1,CTR(1)
            DGHOST(ii-i,jj,kk)=ADATA(ii+i,jj,kk)
        ENDDO
    ELSEIF(ii==DSIZE(1))then
        DO i=1,CTR(1)
            DGHOST(ii+i,jj,kk)=ADATA(ii-i,jj,kk)
        ENDDO
    ENDIF

END DO
END DO
END DO

DO kk=1,DSIZE(3)
DO jj=1,DSIZE(2)
DO ii=1-CTR(1),DSIZE(1)+CTR(1)

IF(CTR(1)==0)exit
IF(CTR(2)==0)exit

    IF(jj==1)then
        DO j=1,CTR(2)
            DGHOST(ii,jj-j,kk)=DGHOST(ii,jj+j,kk)
        ENDDO
    ELSEIF(jj==DSIZE(2))then
        DO j=1,CTR(2)
            DGHOST(ii,jj+j,kk)=DGHOST(ii,jj-j,kk)
        ENDDO
    ENDIF

END DO
END DO
END DO

DO kk=1,DSIZE(3)
DO jj=1-CTR(2),DSIZE(2)+CTR(2)
DO ii=1-CTR(1),DSIZE(1)+CTR(1)

IF(CTR(3)==0)exit
IF(CTR(2)==0)exit
IF(CTR(1)==0)exit

    IF(kk==1)then
        DO k=1,CTR(3)
            DGHOST(ii,jj,kk-k)=DGHOST(ii,jj,kk+k)
        ENDDO
    ELSEIF(kk==DSIZE(3))then
        DO k=1,CTR(3)
            DGHOST(ii,jj,kk+k)=DGHOST(ii,jj,kk-k)
        ENDDO
    ENDIF

END DO
END DO
END DO

END SUBROUTINE
!==================================================================
END MODULE 
