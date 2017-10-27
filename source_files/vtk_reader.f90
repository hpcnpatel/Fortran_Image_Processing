MODULE vtkreader

    USE CONSTANTS_MODULE 
    USE VARIABLE_MODULE 
    USE READ_A_FILE 

IMPLICIT NONE

CONTAINS

!====================================================================
!====================================================================
SUBROUTINE READ_VTK_FILE(ct,arr,file_name)

CHARACTER(c_len)                                          :: dataset_TYPE
CHARACTER(c_len),INTENT(INOUT)                            ::file_name
CHARACTER(c_len),allocatable,dimension(:)                 ::content_arr
INTEGER(kind=rk4)                                         ::count
Type(tCt)                                                 :: ct
CHARACTER(c_len)                                          ::temp
INTEGER                                                   ::i
LOGICAL                                                   ::EXT
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(INOUT) ::arr

WRITE(*,'(10(" "),A)')'* Opening .vtk file to read content in it ...'

INQUIRE(FILE=file_name,EXIST=EXT)

IF(EXT .eqv. .FALSE.) WRITE(*,*)'==== NO FILE WITH A FILENAME, "',trim(file_name),'" found! ===='
IF(.NOT. EXT) STOP

CALL read_no_of_lines_in_file(file_name,count)

!=======================
!read vtk file for determining dataset type
CALL read_content_in_file(file_name,count,content_arr)

    DO i=1,count

        temp=trim(content_arr(i))

        IF(temp(1:7) == 'DATASET')then
        READ(temp(8:),*)dataset_TYPE
        END IF

    ENDDO
        dataset_TYPE=trim(ADJUSTL(dataset_TYPE))
        write(*,*)
        WRITE(*,'(10(" "),A,A)')'==> DATASET TYPE of VTK file is :: ',dataset_TYPE
!=======================

    IF (dataset_TYPE == 'STRUCTURED_POINTS')THEN
        CALL HEADER_STRUCTURED_POINTS(ct,arr,content_arr,file_name)
!    ELSEIF (DATASET == 'STRUCTURED_GRID')THEN
!        CALL HEADER_STRUCTURED_GRID()
!    ELSEIF (DATASET == 'UNSTRUCTURED_GRID')THEN
!        CALL HEADER_UNSTRUCTURED_GRID()
    ELSEIF (DATASET_TYPE == 'POLYDATA')THEN
    WRITE(*,*)'POLYDATA vtk file formate not supported. Update the VTK polydata reader'
    ! update how to extract dimesions of the structured data set from the
    ! polydata file format. 
    STOP
!        CALL HEADER_POLYDATA(ct,arr,content_arr,file_name)
!    ELSEIF (DATASET == 'RECTILINEAR_GRID')THEN
!        CALL HEADER_RECTILINEAR_GRID()
!    ELSEIF (DATASET == 'FIELD')THEN
!        CALL HEADER_FIELD()
    ELSE
        WRITE(*,*)'==== DATASET TYPE IS NOT FOULD IN VTK FILE, &
            &CHECK VTK FILE AND VTK_READER.F90 FILE ====' 
    ENDIF

END SUBROUTINE READ_VTK_FILE
!====================================================================
!====================================================================
SUBROUTINE HEADER_POLYDATA(ct,arr,content_arr,file_name)

CHARACTER(c_len)                                        :: data_type
CHARACTER(c_len),INTENT(INOUT)                          ::file_name
CHARACTER(c_len),allocatable,dimension(:),INTENT(IN)    ::content_arr
INTEGER(kind=rk4)                                       ::count
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::temp
INTEGER                                                 ::tot_points
INTEGER                                                 ::i
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(INOUT) ::arr

count=size(content_arr)

    DO i=1,count

        temp=trim(content_arr(i))

    IF(temp(1:6) == 'POINTS')THEN
        READ(temp(7:),*)tot_points,data_type
        WRITE(*,*)'*******************************************'
        WRITE(*,*)'**** tot no of points(size of all dimesions) ',tot_points
        WRITE(*,*)'**** datatype is :: ',data_type
    END IF

    ENDDO

!=======================
data_type=trim(data_type)

If (data_type == 'short')then
    !Int kind_type=2
    CALL READ_CONTENT_VTK_OF_FILE_INT_2(ct,arr,file_name)
ELSEIF(data_type == 'int')then
    !Int kind_type=4
    CALL READ_CONTENT_VTK_OF_FILE_INT_4(ct,arr,file_name)
ELSEIF(data_type == 'long')then
    !Int kind_type=8
    CALL READ_CONTENT_VTK_OF_FILE_INT_8(ct,arr,file_name)
ELSEIF(data_type == 'float')then
    !Real kind_type=4
    CALL READ_CONTENT_VTK_OF_FILE_REAL_4(ct,arr,file_name)
ELSEIF(data_type == 'double')then
    !Real kind_type=8
    CALL READ_CONTENT_VTK_OF_FILE_REAL_8(ct,arr,file_name)
ELSE
    WRITE(*,*)'CHECK YOUR VTK FILE AND VTK_READER.F90 file there is something wrong with datatype definition'
ENDIF

END SUBROUTINE HEADER_POLYDATA
!====================================================================
SUBROUTINE HEADER_STRUCTURED_POINTS(ct,arr,content_arr,file_name)

CHARACTER(c_len)                                        :: data_name
CHARACTER(c_len)                                        :: data_type
CHARACTER(c_len),INTENT(INOUT)                          ::file_name
CHARACTER(c_len),allocatable,dimension(:),INTENT(IN)    ::content_arr
INTEGER(kind=rk4)                                       ::count
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::temp
INTEGER                                                 ::i
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(INOUT) ::arr

count=size(content_arr)

    DO i=1,count

        temp=trim(content_arr(i))

    IF(temp(1:10) == 'DIMENSIONS')THEN
        READ(temp(11:),*)ct%pixel_x,ct%pixel_y,ct%slices
        WRITE(*,'(10(" "),A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(10(" "),A,I0,A)')'***** Num. of Pixel in y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(10(" "),A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
    END IF

    IF(temp(1:7) == 'SPACING')THEN
        READ(temp(8:),*)ct%dx,ct%dy,ct%dz
        WRITE(*,'(10(" "),A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(10(" "),A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(10(" "),A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
    END IF

    IF(temp(1:6) == 'ORIGIN')THEN
        READ(temp(7:),*)ct%offset_x,ct%offset_y,ct%offset_z
        WRITE(*,'(10(" "),A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(10(" "),A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(10(" "),A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'
    END IF

    ! We will find what attribute is used for data type format. 
    IF(temp(1:7) == 'SCALARS')THEN
        READ(temp(8:),*)data_name, data_type
        WRITE(*,'(10(" "),A,A,A)')'***** data name is ::',trim(data_name),' *****'
        WRITE(*,'(10(" "),A,A,A)')'***** data type is ::',trim(data_type),' *****'
    END IF

END DO

ct%x%a=1
ct%y%a=1
ct%z%a=1

ct%x%e=ct%pixel_x
ct%y%e=ct%pixel_y
ct%z%e=ct%slices

!=======================
data_type=trim(data_type)

If (data_type == 'short')then
    !Int kind_type=2
    CALL READ_CONTENT_VTK_OF_FILE_INT_2(ct,arr,file_name)
ELSEIF(data_type == 'int')then
    !Int kind_type=4
    CALL READ_CONTENT_VTK_OF_FILE_INT_4(ct,arr,file_name)
ELSEIF(data_type == 'long')then
    !Int kind_type=8
    CALL READ_CONTENT_VTK_OF_FILE_INT_8(ct,arr,file_name)
ELSEIF(data_type == 'float')then
    !Real kind_type=4
    CALL READ_CONTENT_VTK_OF_FILE_REAL_4(ct,arr,file_name)
ELSEIF(data_type == 'double')then
    !Real kind_type=8
    CALL READ_CONTENT_VTK_OF_FILE_REAL_8(ct,arr,file_name)
ELSE
    WRITE(*,*)' ==== CHECK YOUR VTK FILE AND VTK_READER.F90 file there is something wrong with datatype definition ===='
ENDIF

END SUBROUTINE HEADER_STRUCTURED_POINTS
!====================================================================
SUBROUTINE READ_CONTENT_VTK_OF_FILE_INT_2(ct,array,file_name)

INTEGER,PARAMETER                                       ::rk_var=2
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::file_name
INTEGER(kind=rk_var),allocatable,dimension(:,:,:)       ::tarray
INTEGER(kind=rk4)                                       ::fs,hs
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(OUT) ::array

    ALLOCATE(array(ct%pixel_x,ct%pixel_y,ct%slices))
    array=0.0
    ALLOCATE(tarray(ct%pixel_x,ct%pixel_y,ct%slices))
    tarray=0

    INQUIRE(file=file_name,size=fs)
    hs=(fs-1)-(ct%pixel_x*ct%pixel_y*ct%slices*rk_var)

WRITE(*,'(/,10(" "),A,I0)')'==> Header size ::',hs
WRITE(*,'(10(" "),A,I0)')'==> File size ::',fs

OPEN(12,FILE=file_name,ACTION='READ',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
    READ(12,pos=hs+1)tarray
CLOSE(12)

array=real(tarray,8)

END SUBROUTINE
!====================================================================
SUBROUTINE READ_CONTENT_VTK_OF_FILE_INT_4(ct,array,file_name)

INTEGER,PARAMETER                                       ::rk_var=4
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::file_name
INTEGER(kind=rk_var),allocatable,dimension(:,:,:)       ::tarray
INTEGER(kind=rk4)                                       ::fs,hs
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(OUT) ::array

    ALLOCATE(array(ct%pixel_x,ct%pixel_y,ct%slices))
    array=0.0
    ALLOCATE(tarray(ct%pixel_x,ct%pixel_y,ct%slices))
    tarray=0

    INQUIRE(file=file_name,size=fs)
    hs=(fs-1)-(ct%pixel_x*ct%pixel_y*ct%slices*rk_var)

WRITE(*,'(/,10(" "),A,I0)')'==> Header size ::',hs
WRITE(*,'(10(" "),A,I0)')'==> File size ::',fs

OPEN(12,FILE=file_name,ACTION='READ',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
    READ(12,pos=hs+1)tarray
CLOSE(12)

array=real(tarray,8)

END SUBROUTINE
!====================================================================
SUBROUTINE READ_CONTENT_VTK_OF_FILE_INT_8(ct,array,file_name)

INTEGER,PARAMETER                                       ::rk_var=8
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::file_name
INTEGER(kind=rk_var),allocatable,dimension(:,:,:)       ::tarray
INTEGER(kind=rk4)                                       ::fs,hs
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(OUT) ::array

    ALLOCATE(array(ct%pixel_x,ct%pixel_y,ct%slices))
    array=0.0
    ALLOCATE(tarray(ct%pixel_x,ct%pixel_y,ct%slices))
    tarray=0

    INQUIRE(file=file_name,size=fs)
    hs=(fs-1)-(ct%pixel_x*ct%pixel_y*ct%slices*rk_var)

WRITE(*,'(/,10(" "),A,I0)')'==> Header size ::',hs
WRITE(*,'(10(" "),A,I0)')'==> File size ::',fs

OPEN(12,FILE=file_name,ACTION='READ',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
    READ(12,pos=hs+1)tarray
CLOSE(12)

array=real(tarray,8)

END SUBROUTINE
!====================================================================
SUBROUTINE READ_CONTENT_VTK_OF_FILE_REAL_4(ct,array,file_name)

INTEGER,PARAMETER                                       ::rk_var=4
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::file_name
REAL(kind=rk_var),allocatable,dimension(:,:,:)          ::tarray
INTEGER(kind=rk4)                                       ::fs,hs
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(OUT) ::array

    ALLOCATE(array(ct%pixel_x,ct%pixel_y,ct%slices))
    array=0.0
    ALLOCATE(tarray(ct%pixel_x,ct%pixel_y,ct%slices))
    tarray=0.0

    INQUIRE(file=file_name,size=fs)
    hs=(fs-1)-(ct%pixel_x*ct%pixel_y*ct%slices*rk_var)

WRITE(*,'(/,10(" "),A,I0)')'==> Header size ::',hs
WRITE(*,'(10(" "),A,I0)')'==> File size ::',fs

OPEN(12,FILE=file_name,ACTION='READ',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
    READ(12,pos=hs+1)tarray
CLOSE(12)

array=real(tarray,8)

END SUBROUTINE
!====================================================================
SUBROUTINE READ_CONTENT_VTK_OF_FILE_REAL_8(ct,array,file_name)

INTEGER,PARAMETER                                       ::rk_var=8
Type(tCt)                                               :: ct
CHARACTER(c_len)                                        ::file_name
REAL(kind=rk_var),allocatable,dimension(:,:,:)          ::tarray
INTEGER(kind=rk4)                                       ::fs,hs
REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(OUT) ::array

    ALLOCATE(array(ct%pixel_x,ct%pixel_y,ct%slices))
    array=0.0
    ALLOCATE(tarray(ct%pixel_x,ct%pixel_y,ct%slices))
    tarray=0.0

    INQUIRE(file=file_name,size=fs)
    hs=(fs-1)-(ct%pixel_x*ct%pixel_y*ct%slices*rk_var)

WRITE(*,'(/,10(" "),A,I0)')'==> Header size ::',hs
WRITE(*,'(10(" "),A,I0)')'==> File size ::',fs

OPEN(12,FILE=file_name,ACTION='READ',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
    READ(12,pos=hs+1)tarray
CLOSE(12)

array=real(tarray,8)

END SUBROUTINE
!====================================================================
END MODULE vtkreader
