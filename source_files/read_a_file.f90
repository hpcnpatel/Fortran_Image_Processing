MODULE READ_A_FILE

USE CONSTANTS_MODULE
USE GENERAL_FUN_SUB

IMPLICIT NONE

CONTAINS

!================================================
        SUBROUTINE read_no_of_lines_in_file(file_name,count)

        CHARACTER(len=*),INTENT(IN)      ::file_name
        INTEGER,INTENT(OUT)              ::count
        INTEGER                          ::ios  
        INTEGER                          ::unit

    unit=provide_new_unit()

    OPEN(unit,FILE=trim(file_name),ACTION='read',STATUS='old',IOSTAT=ios)

IF(ios >= 1)WRITE(*,*)'Error while reading file',file_name
IF(ios >= 1)STOP

        count=0
        DO WHILE(ios == 0)
        !write(*,*)count
        read(200,'(A)',IOSTAT=ios)
        count=count+1
        END DO

        count=count-1
    CLOSE(200)

        END SUBROUTINE
!================================================
        SUBROUTINE read_content_in_file(file_name,count,content_arr)

        CHARACTER(len=*  ),INTENT(IN)                          ::file_name
        INTEGER,INTENT(OUT)                                    ::count
        INTEGER                                                ::ios
        INTEGER                                                ::i
        INTEGER                                                ::unit
        CHARACTER(c_len),allocatable,dimension(:),INTENT(OUT)::content_arr

    unit=provide_new_unit()

    OPEN(unit,FILE=trim(file_name),ACTION='read',STATUS='old',IOSTAT=ios)
IF(ios >= 1)WRITE(*,*)'Error while reading file',file_name
IF(ios >= 1)STOP
        count=0
        DO WHILE(ios == 0)
        read(unit,'(A)',IOSTAT=ios)
            count=count+1
        END DO
        count=count-1
    CLOSE(unit)

        ALLOCATE(content_arr(1:count))

        content_arr=''

    unit=provide_new_unit()

    OPEN(unit,FILE=trim(file_name),ACTION='read',STATUS='old',IOSTAT=ios)
        DO i=1,count
            READ(unit,'(A)',IOSTAT=ios)content_arr(i)
        END DO
    CLOSE(unit)

        END SUBROUTINE
!================================================

END MODULE 
