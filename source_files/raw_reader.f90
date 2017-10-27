module rawreader

    USE CONSTANTS_MODULE
    USE VARIABLE_MODULE 

  IMPLICIT NONE

 CONTAINS

  !****************************************************************************
  !*                                                                         **
  !* Subroutine to load header information of CT-data file                   **
  !*                                                                         **
  !* last edited : on  18.08.2014 by :Nisarg                                 **
  !*                                                                         **
  !****************************************************************************
  Subroutine READ_RAW_FILE(ct,ct_data,file_name)

    Type(tCT), Intent(Out)       :: ct
    Character(c_len)             :: tmp
    Character(Len=*),Parameter   :: fmt_line  = '(A256)'
    Character(Len=*),Parameter   :: fmt_A8    = '(A8)'
    Character(Len=*),Parameter   :: fmt_A10   = '(A10)'
    Character(Len=*),Parameter   :: fmt_A12   = '(A12)'
    Character(Len=*),Parameter   :: fmt_A13   = '(A13)'
    Character(Len=*),Parameter   :: fmt_A18   = '(A18)'
    Character(Len=*),Parameter   :: fmt_AF12  = '(A,F12.6)'
    Character(Len=*),Parameter   :: fmt_AI    = '(A,I6)'
    Character(Len=*),Parameter   :: fmt_2AI   = '(2(A,I6))'
    Character(Len=*),Parameter   :: fmt_I12A  = '(I12,A)'
    Character(Len=*),Parameter   :: fmt_shift = '(A,F12.6)'
    Character(Len=*),Parameter   :: fmt_field_lim = '(2(A,I5),A)'
    Character(Len=*),PARAMETER   :: fmt_str   = '(A)'
    Character(Len=*),PARAMETER   :: fmt_sep   = "(/,80('='),/)"

    Integer,Parameter            :: un_hd  = 32 !** Unit number CT-header-file
    Integer,Parameter            :: un_raw = 132!** Unit number CT-binary data file
    Integer,Parameter            :: un_lf  = 232!** Unit number log-file

    Character(c_len),INTENT(INOUT):: file_name ! Name of CT-header-file
    Character(c_len)             :: data_file   ! Name of CT-data-file

    Integer                                                   ::io_stat
    REAL(kind=rk8),allocatable,dimension(:,:,:),INTENT(INOUT) ::ct_data
    Integer                                                   ::ii, jj, kk
    Integer                                                   ::no_data = 0

    !==========================================================================
WRITE(*,*)'    | =================================|'
WRITE(*,*)'    |*.raw file will be read. To read  |'
WRITE(*,*)'    | *.raw files two files with binary|'
WRITE(*,*)'    | and header info is required      |'
WRITE(*,*)'    | =================================|'

data_file=trim(file_name)//'.bin'
file_name=trim(file_name)//'.dat'

write(*,*)'***** name of header file:: ',trim(file_name),' ****'

    !==========================================================================
    Open(UNIT=un_hd, FILE=Trim(file_name), STATUS='OLD', ACTION='READ', &
         IOSTAT=io_stat)

    Do While (io_stat == 0)

       Read(un_hd, fmt_line, IOSTAT=io_stat) tmp

       If (tmp(1:23) == '# Voxel und Dimensionen') Then

          Read(un_hd, *, IOSTAT=io_stat) tmp
          Read(un_hd, *, IOSTAT=io_stat) tmp

          Read(un_hd, fmt_A8, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%slices

          Read(un_hd, fmt_A10, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%dx

          Read(un_hd, fmt_A10, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%dy

          Read(un_hd, fmt_A10, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%dz

          ct%dz=abs(ct%dz)

       Else If (tmp(1:31) == '# Daten des Fensterausschnitts:') Then

          Read(un_hd, fmt_A13, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%x%a

          Read(un_hd, fmt_A13, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%y%a

          Read(un_hd, fmt_A12, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%pixel_x

          Read(un_hd, fmt_A13, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%pixel_y  

       Else If (tmp(1:32) == '# Offset zu Scanner Koordinaten:') Then
          
          Read(un_hd, fmt_A18, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%offset_x

          Read(un_hd, fmt_A18, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%offset_y

          Read(un_hd, fmt_A18, IOSTAT=io_stat, ADVANCE='No') tmp
          Read(un_hd, *, IOSTAT=io_stat) ct%offset_z

       End If

    End Do

    Close(UNIT=un_hd, STATUS='KEEP') 

    ct%x%a=1
    ct%y%a=1
    ct%z%a=1

    ct%x%e=ct%pixel_x
    ct%y%e=ct%pixel_y
    ct%z%e=ct%slices


        WRITE(*,*)'*******************************************'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
        WRITE(*,*)'*******************************************'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
        WRITE(*,*)'*******************************************'
        WRITE(*,'(A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'
        WRITE(*,*)'*******************************************'

    ALLOCATE(ct_data(ct%pixel_x,ct%pixel_y,ct%slices))

write(*,*)'***** RAW file is being read, name of data file:: ',trim(data_file),' ****'

    !** Open binary input file 
    
      open(un_raw,file=data_file,action='read',status='old',access='stream')
           
       Do kk = ct%z%a, ct%z%e
          Do jj = ct%y%a, ct%y%e
             Do ii = ct%x%a, ct%x%e
                
                Read(un_raw) ct_data(ii,jj,kk)
                no_data = no_data + 1

             End Do
          End Do
       End Do
   
    Close(un_raw)

!WRITE(*,*)lbound(ct_data),ubound(ct_data)

  End Subroutine READ_RAW_FILE
  !****************************************************************************
END MODULE rawreader
