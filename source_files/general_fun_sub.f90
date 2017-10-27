MODULE GENERAL_FUN_SUB

IMPLICIT NONE

CONTAINS
!===============================================
FUNCTION provide_new_unit()result(new_unit)

INTEGER :: new_unit
LOGICAL :: log_openornot
LOGICAL :: log_flag
INTEGER :: ii

ii=200
DO

    INQUIRE(UNIT=ii,OPENED=log_openornot)

    IF(.NOT. log_openornot)THEN
    new_unit=ii
    log_flag=.TRUE.
    ELSE
    ii=ii+1
    log_flag=.FALSE.
    ENDIF

IF(log_flag)EXIT
END DO

END FUNCTION
!===============================================
END MODULE
!===============================================
