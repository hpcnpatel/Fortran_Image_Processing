PROGRAM IMAGE_PROCESS
!
!
!
!
!
!

    USE OMP_LIB
    USE UTILITIES_MODULE; USE VARIABLE_MODULE; USE CONSTANTS_MODULE
    USE VTKREADER; USE VTKWRITTER; USE RAWREADER; USE binaryreader
    USE READ_A_FILE
    USE VANDERMONDE_MODULE
    USE GHOST; USE INTERPOLATE_CTDATA
    USE MOD_QR; USE MOD_MEDIAN; USE MOD_grad_of_image
    !USE FILTERS; 
    USE MOD_GAUSS_FILTER; !USE SNAKE
    USE ISO_SURFACE
    USE SORT_MODULE
    use MOD_GNUPLOT
    use fast_fourier_transform_module

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
!FILE_NAME=trim(ADJUSTL(FILE_NAME))

WRITE(*,'(/,"*",80("="),"*",/)')
WRITE(*,'(10(" "),35("*"))')
WRITE(*,'(10(" "),A,A)')'* Reading File ::',trim(FILE_NAME)
WRITE(*,'(10(" "),35("*"),/)')

!==========================================================
CALL READ_VTK_FILE(CT,CT_DATA,FILE_NAME)
    WRITE(FILE_NAME,'(A,A,A)')'gray_value_',trim(FNAME),'.vtk'
    CALL WRITE_VTK(CT,CT_DATA,FILE_NAME)
!FILE_NAME='/home/hpcnpate/HLRS/pixel_files/torax/torax_10f/pixel'
!CALL READ_RAW_FILE(CT,CT_DATA,FILE_NAME)
!CALL READ_binary_file(dim1d,DATA,FILE_NAME)
!==============================================================
!     CALL IN_DATA_3D_Z_ZERO(CT,CT_DATA,CT_NEW,FILE_NAME)
!    CALL IN_DATA_3D_Z_one(CT,CT_DATA,CT_NEW,file_name)

!CT_NEW=CT_DATA
!    CALL ADAPTIVE_MEDIAN(CT_NEW,KSIZE,np)
!    WRITE(FILE_NAME,'(A,I0,A,A)')'adap_median_',KSIZE(1),trim(FNAME),'.vtk'
!    CALL WRITE_VTK(CT,CT_NEW,FILE_NAME)
!CT_NEW=CT_DATA
!    CALL ADAPTIVE_QR(CT_NEW,KSIZE,np)
!    WRITE(FILE_NAME,'(A,I0,A,A)')'adap_qr_',KSIZE(1),trim(FNAME),'.vtk'
!    CALL WRITE_VTK(CT,CT_NEW,FILE_NAME)
!    KSIZE=(/3,3,1/)
CT_NEW=CT_DATA
    CALL GRADIENT(CT_NEW,GRAD_VEC,KSIZE,np)
    WRITE(FILE_NAME,'(A,I0,A,A)')'grad_',KSIZE(1),trim(FNAME),'.vtk'
    CALL WRITE_VTK(CT,CT_NEW,FILE_NAME)
!    CALL GAUSS(CT_DATA,KSIZE,np)
!    CALL MEAN(CT_NEW,KSIZE,np)

!DSIZE(1)=size(CT_DATA,1)
!!DSIZE(2)=size(CT_DATA,2)
!DSIZE(3)=size(CT_DATA,3)

!ALLOCATE(TEMP(DSIZE(1)*DSIZE(2)*DSIZE(3)))
!ALLOCATE(TEMP1(DSIZE(1)*DSIZE(2)*DSIZE(3)))
!a=0
!DO k=1,DSIZE(3)
!DO j=1,DSIZE(2)
!DO i=1,DSIZE(1)
!a=a+1
!temp(a)=ct_data(i,j,k)
!END DO
!END DO
!END DO

!CALL fft_real_1D(TEMP,TEMP1,-1) 

!a=0
!DO k=1,DSIZE(3)
!DO j=1,DSIZE(2)
!DO i=1,DSIZE(1)
!a=a+1
!ct_data(i,j,k)=temp1(a)
!END DO
!END DO
!END DO
!    WRITE(FILE_NAME,'(A,A,A)')'FFT_',trim(FNAME),'.vtk'
!    CALL WRITE_VTK(CT,CT_DATA,FILE_NAME)
!CALL gnu_initial_script()
!CALL vector_plot_array_real(GRAD_VEC,filename=fname)
!CALL surface_plot_array_real(ct_data,filename='data')
!CALL surface_plot_array_real(ct_data,filename='data_grad')


!==============================================================
END PROGRAM IMAGE_PROCESS
!==============================================================
!TEMP=RESHAPE(CT_DATA,(/size(CT_DATA)/))
!CALL SORTRX_real8(TEMP,SORT_INDEX)
!TEMP=TEMP(SORT_INDEX)
!    CALL MEDIAN(CT_NEW,KSIZE,np)
!    WRITE(FILE_NAME,'(A,I0,A,A,A)')'median_',KSIZE(1),'_',trim(FNAME),'.vtk'
!    CALL ADAPTIVE_MEDIAN(CT_NEW,KSIZE,np)
!    WRITE(FILE_NAME,'(A,I0,A,A,A)')'adaptive_median_',KSIZE(1),'_',trim(FNAME),'.vtk'
!    CALL QR(CT_NEW,KSIZE,np)
