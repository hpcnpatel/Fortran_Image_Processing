module binaryreader

    USE CONSTANTS_MODULE
    USE VARIABLE_MODULE 

  IMPLICIT NONE

 CONTAINS

  Subroutine READ_binary_file(dim1d,data,file_name)

    Integer,Parameter            :: un_hd  = 32 !** Unit number CT-header-file
    Integer,Parameter            :: un_raw = 132!** Unit number CT-binary data file
    Integer,Parameter            :: un_lf  = 232!** Unit number log-file

    Character(c_len),INTENT(INOUT):: file_name ! Name of CT-header-file

    Integer                                                   ::io_stat,dim1d
CHARACTER(C_LEN)                                              ::data_file
    REAL(kind=rk8),allocatable,dimension(:),INTENT(OUT)       ::data
    REAL(kind=4),allocatable,dimension(:)       ::temp
    Integer                                                   ::ii

    !==========================================================================
WRITE(*,*)'    | =================================|'
WRITE(*,*)'    |*.dat file will be read.          |'
WRITE(*,*)'    | This is binary file with floats  |'
WRITE(*,*)'    | =================================|'

!data_file=trim(file_name)//'.dat'
data_file=trim(file_name)

!!!!!write(*,*)'***** RAW file is being read'
write(*,*)'***** name of data file:: ',trim(file_name),' ****'

    !==========================================================================
    dim1d=25000
    ALLOCATE(data(1:dim1d))
    ALLOCATE(temp(1:dim1d))

    Open(UNIT=un_raw, FILE=Trim(file_name), STATUS='OLD', ACTION='READ', ACCESS='STREAM') 

!    Open(UNIT=un_hd, FILE=Trim(file_name), STATUS='OLD', ACTION='READ', &
!         IOSTAT=io_stat)
!    Do While (io_stat == 0)
                Do ii=1,dim1d
                Read(un_raw) temp(ii)
!                no_data = no_data + 1
                END DO
!    End Do
!    Close(UNIT=un_hd, STATUS='KEEP') 
    Close(UNIT=un_raw) 
data=real(temp,8)

  End Subroutine READ_binary_file
  !****************************************************************************
END MODULE binaryreader
