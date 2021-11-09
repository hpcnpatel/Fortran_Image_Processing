PROGRAM IMAGE_PROCESS

    USE OMP_LIB
    USE UTILITIES_MODULE; USE VARIABLE_MODULE; USE CONSTANTS_MODULE
    USE VTKREADER; USE VTKWRITTER; USE RAWREADER; USE binaryreader
    USE READ_A_FILE
    USE VANDERMONDE_MODULE
    USE GHOST; USE INTERPOLATE_CTDATA
    USE MOD_QR; USE MOD_MEDIAN; USE MOD_grad_of_image
    USE MOD_GAUSS_FILTER; !USE SNAKE
    USE ISO_SURFACE
    USE SORT_MODULE
    USE MOD_GNUPLOT
    USE FAST_FOURIER_TRANSFORM_MODULE

IMPLICIT NONE

!=========================================================
!variables and types
!=========================================================
TYPE(TCT)                                                 ::CT
REAL(KIND=RK8),DIMENSION(:,:,:),ALLOCATABLE               ::CT_DATA
REAL(KIND=RK8),DIMENSION(:,:,:),ALLOCATABLE               ::CT_NEW
REAL(KIND=RK8),DIMENSION(:,:),ALLOCATABLE               ::V_ARR
REAL(KIND=RK8),DIMENSION(:),ALLOCATABLE                   ::TEMP,TEMP1
INTEGER,ALLOCATABLE,DIMENSION(:)                          ::SORT_INDEX
INTEGER                                                   ::count
CHARACTER(C_LEN),ALLOCATABLE,DIMENSION(:)                 ::Abuffer
CHARACTER(C_LEN)                                          ::buffer
CHARACTER(C_LEN)                                          ::CONTROL_PARAMETER_FILE
INTEGER                                                   ::i,j,k,a
TYPE(vector),DIMENSION(:,:,:),ALLOCATABLE                 ::GRAD_VEC
!=========================================================
!variables and types from CONTROL FILE
!=========================================================
INTEGER,DIMENSION(3)                                      ::KSIZE
INTEGER,DIMENSION(3)                                      ::DSIZE
INTEGER                                                   ::MED_LEN
INTEGER                                                   ::QR_LEN
INTEGER                                                   ::GAUSS_LEN
INTEGER                                                   ::np
CHARACTER(C_LEN)                                          ::FILE_NAME
CHARACTER(C_LEN)                                          ::FNAME,FPATH,FEXT
!=========================================================

CALL GET_COMMAND_ARGUMENT(1,CONTROL_PARAMETER_FILE)
CALL read_content_in_file(CONTROL_PARAMETER_FILE,count,ABUFFER)

DO i=1,count

BUFFER=ADJUSTL(ABUFFER(i))

IF(BUFFER(1:9) == 'kernel_X=')read(BUFFER(10:),*)KSIZE(1)
IF(BUFFER(1:9) == 'kernel_Y=')read(BUFFER(10:),*)KSIZE(2)
IF(BUFFER(1:9) == 'kernel_Z=')read(BUFFER(10:),*)KSIZE(3)
IF(BUFFER(1:6) == 'fpath='   )read(BUFFER(7: ),'(A)')FPATH
IF(BUFFER(1:6) == 'fname='   )read(BUFFER(7: ),'(A)')FNAME
IF(BUFFER(1:5) == 'fext='    )read(BUFFER(6: ),'(A)')FEXT
IF(BUFFER(1:3) == 'np='      )read(BUFFER(4: ),*)np

END DO


FILE_NAME=trim(fpath)//trim(fname)//'.'//trim(fext)

WRITE(*,'(/,"*",80("="),"*",/)')
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A,A)')'* Reading File ::',trim(FILE_NAME)
WRITE(*,'(10(" "),35("*"),/)')

!==========================================================
IF (trim(fext) == 'vtk') THEN
    CALL READ_VTK_FILE(CT,CT_DATA,FILE_NAME)
    WRITE(*,*)  'vtk'
ELSE IF (trim(fext) == 'raw') THEN
    CALL READ_RAW_FILE(CT,CT_DATA,FILE_NAME)
    WRITE(*,*)  'RAW'
ELSE
    WRITE(*,*)  'The program can read RAW or VTK files. Please specify either type of file name'
   STOP
END IF
    WRITE(FILE_NAME,'(A,A,A)')'gray_value_',trim(FNAME),'.vtk'
    CALL WRITE_VTK(CT,CT_DATA,FILE_NAME)
!==============================================================
    CALL IN_DATA_3D_Z_ZERO(CT,CT_DATA,CT_NEW,FILE_NAME)
!    CALL IN_DATA_3D_Z_one(CT,CT_DATA,CT_NEW,file_name)
!    CALL IN_DATA_3D_Z_two(CT,CT_DATA,CT_NEW,file_name)

    CALL ADAPTIVE_MEDIAN(CT_NEW,KSIZE,np)
    WRITE(FILE_NAME,'(A,I0,A,A)')'median_',KSIZE(1),trim(FNAME),'.vtk'
    CALL WRITE_VTK(CT,CT_NEW,FILE_NAME)

!! Other filters (subroutines) that can be called.    

!    CALL ADAPTIVE_QR(CT_NEW,KSIZE,np)
!    CALL GRADIENT(CT_NEW,GRAD_VEC,KSIZE,np)
!    CALL MEAN(CT_DATA,KSIZE,np)
!    CALL MEDIAN(CT_NEW,(/9,9,9/),np)

!==============================================================
END PROGRAM IMAGE_PROCESS
!==============================================================
