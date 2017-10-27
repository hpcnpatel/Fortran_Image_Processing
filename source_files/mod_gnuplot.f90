module mod_gnuplot

USE constants_module
USE variable_module
USE GENERAL_FUN_SUB

IMPLICIT NONE

!=============================================================
!interface write_data
!   module procedure write_data_real_1D
!   module procedure write_data_real_2D
!end interface
!=============================================================

CONTAINS
!=====================================================================
subroutine gnu_initial_script(Xlabel,Ylabel,title)
!suppliment data on which the Xdata is to be plotted if present or upper&lower bonds will be considered. 

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
character(len=*)             ,optional               :: title
character(len=*)             ,optional               :: Xlabel
character(len=*)             ,optional               :: Ylabel
!=====================================================================
integer                                          :: unit_int
!=====================================================================

unit_int=provide_new_unit()

OPEN(unit_int,FILE='initial_gnu.gplt') 

!=====================================================================
!GNU INITIAL PLOT SCRIPT WILL BE WRITTEN
!=====================================================================

           WRITE(unit_int,'(a)') 'set style data lines'
           !WRITE(unit_int,'(a)') 'set multiplot' 
           WRITE(unit_int,'(a)') 'set autoscale' 

        IF(present(title))THEN
           WRITE(unit_int,'(a)') 'set title "',trim(title),'"'
        ELSE
           WRITE(unit_int,'(a)') 'unset title'
        ENDIF

        IF(present(Xlabel))THEN
           WRITE(unit_int,'(a)') 'set xlabel "',trim(Xlabel),'"'
        ELSE
           WRITE(unit_int,'(a)') 'unset xlabel'
        ENDIF

           !WRITE(unit_int,'(A)')'unset grid'
           !WRITE(unit_int,'(a)') 'unset title'
           !WRITE(unit_int,'(a)') 'unset xlabel'
           WRITE(unit_int,'(A)')'unset grid'
           WRITE(unit_int,'(A)')'set palette defined (\'
           WRITE(unit_int,'(A)')'0 "#000090",1 "#000fff",\'
           WRITE(unit_int,'(A)')'2 "#0090ff",3 "#0fffee",\'
           WRITE(unit_int,'(A)')'4 "#90ff70",5 "#ffee00",\'
           WRITE(unit_int,'(A)')'6 "#ff7000",7 "#ee0000",\'
           WRITE(unit_int,'(A)')'8 "#7f0000")'

close(unit_int)

end subroutine gnu_initial_script
!=====================================================================
subroutine plot_array_real(data,Ydata,filename)
!suppliment data on which the Xdata is to be plotted if present or upper&lower bonds will be considered. 

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
real(kind=rk),dimension(:)            ,allocatable   :: data
character(len=*)                                     :: filename
real(kind=rk),dimension(:)   ,optional,allocatable   :: Ydata
!=====================================================================
!TYPE local
!=====================================================================
character(256)                                       :: fname
integer                                              :: unit_int
real(kind=rk)                                        :: Xdata1,Xdata2
real(kind=rk)                                        :: Ydata1,Ydata2
!=====================================================================

            fname=trim(filename)//'.dat'

CALL write_data_real_1d(data,Ydata,trim(fname))

            Xdata1=lbound(data,1)
            Xdata2=ubound(data,1)

        IF(present(Ydata))then
            Ydata1=lbound(data,1)
            Ydata2=ubound(data,1)
        ELSE
            Ydata1=minval(data,1)
            IF(Ydata1 .le. 0.001 .and. Ydata1 .ge. -0.001)Ydata1= 0.0
            Ydata2=maxval(data,1)
        ENDIF

unit_int=provide_new_unit()

OPEN(unit_int,FILE=trim(filename)//'.gplt') 

  !! call delete_file(trim(filename)//'.gnu_script')
  ! call append_plot(trim(filename)//'.gnu_script')

!=====================================================================
!GNU PLOT SCRIPT WILL BE WRITTEN
!=====================================================================
   WRITE(unit_int,'(A,ES13.5,A,ES13.5,A)')'set xrange[',Xdata1,':',Xdata2,']'
   WRITE(unit_int,'(A,ES13.5,A,ES13.5,A)')'set yrange[',Ydata1,':',Ydata2,']'
   WRITE(unit_int,'(A)') 'load "initial_gnu.gplt"' 
   WRITE(unit_int,'(A)') 'plot \' 
   WRITE(unit_int,'(A,A,A)') '"',trim(fname),'" using 0:1 lw 2, \'
   WRITE(unit_int,'(A)')'#'
   WRITE(unit_int,'(a,a,a)') 'pause -1 '

close(unit_int)

end subroutine plot_array_real
!=====================================================================
subroutine write_data_real_1d(Xdata,Ydata,filename)

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
real(kind=rk),dimension(:),allocatable               :: Xdata
real(kind=rk),dimension(:),allocatable,optional      :: Ydata
character(len=*)                                     :: filename
!=====================================================================
!TYPE local
!=====================================================================
integer                                              :: nn
integer                                              :: unit_int
!=====================================================================

unit_int=provide_new_unit()

open(newunit=unit_int,file=filename,STATUS='REPLACE',ACTION='WRITE') 

!call delete_file(filename)

if(present(Ydata)) then
   do nn=1,size(Ydata)
      write(unit_int,'(2e15.7e2)',advance='yes') Xdata(nn),Ydata(nn)
   enddo
else
   do nn=1,size(Xdata)
      write(unit_int,'(e15.7e2)',advance='yes') Xdata(nn)
   enddo
endif

close(unit_int)

end subroutine write_data_real_1D
!=====================================================================
subroutine vector_plot_array_real(data,filename)
!suppliment data on which the data is to be plotted if present or upper&lower bonds will be considered. 

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
type(vector),dimension(:,:,:)         ,allocatable   :: data
character(len=*)                                     :: filename
!=====================================================================
!TYPE local
!=====================================================================
character(256)                                       :: fname
integer                                              :: unit_int
real(kind=rk)                                        :: Xdata1,Xdata2
real(kind=rk)                                        :: Ydata1,Ydata2
!=====================================================================

            fname=trim(filename)//'.dat'

CALL vector_write_data_real_2d(data,trim(fname))

            Xdata1=lbound(data,1)
            Xdata2=ubound(data,1)
            Ydata1=lbound(data,2)
            Ydata2=ubound(data,2)

unit_int=provide_new_unit()

OPEN(unit_int,FILE=trim(filename)//'.gplt') 

  !! call delete_file(trim(filename)//'.gnu_script')
  ! call append_plot(trim(filename)//'.gnu_script')

!=====================================================================
!GNU PLOT SCRIPT WILL BE WRITTEN
!=====================================================================
   WRITE(unit_int,'(A,ES13.5,A,ES13.5,A)')'set xrange[',Xdata1,':',Xdata2,']'
   WRITE(unit_int,'(A,ES13.5,A,ES13.5,A)')'set yrange[',Ydata1,':',Ydata2,']'
   WRITE(unit_int,'(A)') 'load "initial_gnu.gplt"' 
   WRITE(unit_int,'(A)') 'plot \' 
   WRITE(unit_int,'(A,A,A)') '"',trim(fname),'" using 1:2:3:4 with vectors &
   &head size 0.1,20,60 filled lw 2, \'
   WRITE(unit_int,'(A)')'#'
   WRITE(unit_int,'(a,a,a)') 'pause -1 '

close(unit_int)

end subroutine vector_plot_array_real
!=====================================================================
subroutine vector_write_data_real_2d(data,filename)

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
TYPE(vector),dimension(:,:,:),allocatable            :: data
character(len=*)                                     :: filename
!=====================================================================
!TYPE local
!=====================================================================
integer                                              :: nn,mm
integer                                              :: unit_int
real(KIND=RK8)                                       ::maxx,maxy,minx,miny
!=====================================================================

unit_int=provide_new_unit()

open(newunit=unit_int,file=filename,STATUS='REPLACE',ACTION='WRITE') 

!call delete_file(filename)

maxx=maxval(data%x)
minx=minval(data%x)
maxy=maxval(data%y)
miny=minval(data%y)

write(*,*)'MAX and MIN in X',maxx,minx
write(*,*)'MAY and MIN in Y',maxy,miny

   do mm=1,SIZE(data,2)
   do nn=1,SIZE(data,1)
      !write(unit_int,'(e15.7e2)',advance='yes') nn,mm,data(nn,mm,1)%x,data(nn,mm,1)%y
      !write(unit_int,*) nn,mm,data(nn,mm,1)%x/maxx,data(nn,mm,1)%y/maxy
      write(unit_int,*) nn,mm,data(nn,mm,1)%x/maxx,data(nn,mm,1)%y/maxy
   enddo
   enddo

close(unit_int)

end subroutine vector_write_data_real_2D
!=====================================================================
subroutine surface_plot_array_real(data,filename)
!suppliment data on which the data is to be plotted if present or upper&lower bonds will be considered. 

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
Real(kind=rk8),dimension(:,:,:)       ,allocatable   :: data
character(len=*)                                     :: filename
!=====================================================================
!TYPE local
!=====================================================================
character(256)                                       :: fname
integer                                              :: unit_int
real(kind=rk)                                        :: Xdata1,Xdata2
real(kind=rk)                                        :: Ydata1,Ydata2
!=====================================================================

            fname=trim(filename)//'.dat'

CALL surface_write_data_real_2d(data,trim(fname))

            Xdata1=lbound(data,1)
            Xdata2=ubound(data,1)
            Ydata1=lbound(data,2)
            Ydata2=ubound(data,2)

unit_int=provide_new_unit()

OPEN(unit_int,FILE=trim(filename)//'.gplt') 

  !! call delete_file(trim(filename)//'.gnu_script')
  ! call append_plot(trim(filename)//'.gnu_script')

!=====================================================================
!GNU PLOT SCRIPT WILL BE WRITTEN
!=====================================================================
   WRITE(unit_int,'(A,ES13.5,A,ES13.5,A)')'set xrange[',Xdata1,':',Xdata2,']'
   WRITE(unit_int,'(A,ES13.5,A,ES13.5,A)')'set yrange[',Ydata1,':',Ydata2,']'
   WRITE(unit_int,'(A)') 'load "initial_gnu.gplt"' 
   WRITE(unit_int,'(A)') 'plot \' 
   WRITE(unit_int,'(A,A,A)') '"',trim(fname),'" using 1:2:3 with image, \'
   WRITE(unit_int,'(A)')'#'
   WRITE(unit_int,'(a,a,a)') 'pause -1 '

close(unit_int)

end subroutine
!=====================================================================
subroutine surface_write_data_real_2d(data,filename)

!=====================================================================
!TYPE incoming and outgoing
!=====================================================================
real(kind=rk8),dimension(:,:,:),allocatable            :: data
character(len=*)                                     :: filename
!=====================================================================
!TYPE local
!=====================================================================
integer                                              :: nn,mm
integer                                              :: unit_int
!=====================================================================

unit_int=provide_new_unit()

open(newunit=unit_int,file=filename,STATUS='REPLACE',ACTION='WRITE') 

   do mm=1,SIZE(data,2)
   do nn=1,SIZE(data,1)
      write(unit_int,*) nn,mm,data(nn,mm,1)
   enddo
   enddo

close(unit_int)

end subroutine 
!=====================================================================
END MODULE
!PROGRAM test_gnu
!
!USE mod_gnuplot
!USE constants_module
!IMPLICIT NONE
!
!
!real(kind=rk),dimension(:),allocatable               :: Xdata
!
!Xdata=(/100,1000,2000,300,500,700,400,1100,1400,1800,1100,1300,600,800,1900/)
!
!CALL plot_array_real(Xdata,filename='data')
!
!
!END PROGRAM test_gnu

