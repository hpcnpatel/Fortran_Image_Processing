module gnuplot_module
use constants_module
use utilities_module
use sort_module
implicit none

integer                            :: lt_undirected_small=14
integer                            :: lt_undirected=-1
integer                            :: lt_directed_small=2
integer                            :: lt_directed=1 
!!real(kind=rk)                      :: skip_edges_threshould=0.09
!!real(kind=rk)                      :: skip_edges_threshould=0.0
real(kind=rk)                      :: skip_edges_threshould=0.050 
!!real(kind=rk)                      :: linewidth_scale=10.
!!real(kind=rk)                      :: linewidth_scale=1.
!!real(kind=rk)                      :: linewidth_scale=13.
real(kind=rk)                      :: linewidth_scale=4.

type valued_egde_type
  integer                                  :: number
  integer,allocatable,dimension(:,:)       :: index
  integer,allocatable,dimension(:)         :: type
  real(kind=rk),allocatable,dimension(:)   :: value
end type valued_egde_type

type curve_type
   !!character(len=:),allocatable              :: name
   character(len=64)                          :: name=''
   real(kind=rk),dimension(:,:),allocatable  :: data
end type curve_type

type curve_set_type
   integer                                    :: lb=1,ub=0 
   type(curve_type),dimension(:),allocatable  :: curve
end type curve_set_type

interface allocate
    module procedure allocate_valued_egde_type
end interface allocate

interface plot_array_allocatable
   module procedure plot_array_allocatable_real
   module procedure plot_array_allocatable_complex
   module procedure plot_array_allocatable_2d
end interface

interface plot_array
   module procedure plot_array_real
end interface

interface write_data
   module procedure write_data_real_1D
   module procedure write_data_real_2D
end interface

interface put_data
   module procedure put_data_curve_set_1d
   module procedure put_data_curve_set_2d
end interface put_data

interface deallocate
   module procedure deallocate_curve_set
end interface deallocate

contains
!****************************************************
subroutine write_gnu_parameter_file(xmin,ymin,xmax,ymax,zmin,zmax &
                                 & ,x_label,y_label,z_label &
                                 & ,colorvec &
                                 & ,coord_type &
                                 & )

real(kind=rk)                              :: xmin,xmax,ymin,ymax
real(kind=rk),optional                     :: zmin,zmax
character(len=*),optional                  :: z_label
character(len=*)                           :: x_label,y_label
integer                                    :: parameter_unit
character(len=*),optional                  :: coord_type
integer,optional                           :: colorvec


    open(newunit=parameter_unit,file='parameter.gnu_script')

!    if(present(colorvec)) write(parameter_unit,'(a,i0)') 'colored_vector=',colorvec

    write(parameter_unit,'(a,e10.3)') 'xmin=',xmin
    write(parameter_unit,'(a,e10.3)') 'xmax=',xmax
    write(parameter_unit,'(a,e10.3)') 'ymin=',ymin
    write(parameter_unit,'(a,e10.3)') 'ymax=',ymax
    if(present(zmin))then
    write(parameter_unit,'(a,e10.3)') 'zmin=',zmin
    else
          write(parameter_unit,'(a)') 'zmin=-1.'
    endif
    if(present(zmax))then
    write(parameter_unit,'(a,e10.3)') 'zmax=',zmax
    else
          write(parameter_unit,'(a)') 'zmax=+1.'
    endif
    if(present(coord_type)) write(parameter_unit,'(3a)')&
                         & 'coordinate_type="',trim(coord_type),'"'
    write(parameter_unit,'(3a)') 'x_label="',trim(x_label),'"'
    write(parameter_unit,'(3a)') 'y_label="',trim(y_label),'"'
    if(present(z_label)) write(parameter_unit,'(3a)')&
                               & 'z_label="',trim(z_label),'"'

    close(parameter_unit)
end subroutine write_gnu_parameter_file
!****************************************************
subroutine allocate_valued_egde_type(valued_egde,number)
integer                :: number
type(valued_egde_type) :: valued_egde

         if(allocated(valued_egde%index)) deallocate(valued_egde%index)  
         if(allocated(valued_egde%type))  deallocate(valued_egde%type)   
         if(allocated(valued_egde%value)) deallocate(valued_egde%value) 

          valued_egde%number=number

         allocate(valued_egde%index(valued_egde%number,2))
         allocate(valued_egde%type (valued_egde%number))
         allocate(valued_egde%value(valued_egde%number))

end subroutine allocate_valued_egde_type
!****************************************************
subroutine allocate_curve_set(curve_set,lb)
type(curve_set_type)                   :: curve_set
integer                                :: lb,ub

              curve_set%lb=lb
              curve_set%ub=0
              !!allocate(curve_set%curve(curve_set%lb:curve_set%ub))
              allocate(curve_set%curve(curve_set%lb:10))

end subroutine allocate_curve_set
!****************************************************
subroutine deallocate_curve_set(curve_set)
type(curve_set_type)                   :: curve_set
integer                                :: lb,ub
integer                                :: ll    

              curve_set%lb=lb
              curve_set%ub=0
              do ll=lb,ub
                 deallocate(curve_set%curve(ll)%data)
              enddo
                 deallocate(curve_set%curve)
              curve_set%lb=0
              curve_set%ub=-1

end subroutine deallocate_curve_set
!****************************************************
subroutine put_data_curve_set_1d(string,curve_set,data1,data2,name)
character(len=*)                           :: string
type(curve_set_type)                       :: curve_set
real(kind=rk),dimension(:)                 :: data1,data2
character(len=*)                           :: name
integer                                    :: lb,ub
type(curve_type),dimension(:),allocatable  :: curve
logical                                    :: is_alloc
integer                                    :: ll

if(lbound(data1,1)/=lbound(data2,1))then
        write(*,'(*(a))') 'put_data_curve_set_1d: data1 and data2 do not conform at ',trim(string)
        write(*,'(*(a,i0))') 'lbound(data1,1)=',lbound(data1,1),' ubound(data1,1)=',ubound(data1,1)
        write(*,'(*(a,i0))') 'lbound(data2,1)=',lbound(data2,1),' ubound(data2,1)=',ubound(data2,1)
        stop 'put_data_curve_set_1d: data1 and data2 do not conform'
endif
if(ubound(data1,1)/=ubound(data2,1))then
        write(*,'(*(a))') 'put_data_curve_set_1d: data1 and data2 do not conform at ',trim(string)
        write(*,'(*(a,i0))') 'lbound(data1,1)=',lbound(data1,1),' ubound(data1,1)=',ubound(data1,1)
        write(*,'(*(a,i0))') 'lbound(data2,1)=',lbound(data2,1),' ubound(data2,1)=',ubound(data2,1)
        stop 'put_data_curve_set_1d: data1 and data2 do not conform'
endif

   lb=curve_set%lb 
   ub=curve_set%ub 

   is_alloc=allocated(curve_set%curve) 
   if(is_alloc) then
      if(ub+1<=ubound(curve_set%curve,1) ) then
         ub=ub+1
      else
         allocate(curve(lb:ub+10)) 
         do ll=lb,ub
            curve(ll)=curve_set%curve(ll)
         enddo
         call move_alloc(curve,curve_set%curve)
         ub=ub+1
      endif
   else
      call allocate_curve_set(curve_set,1)
      lb=lbound(curve_set%curve,1)
      ub=lb
   endif

      curve_set%curve(ub)%name=name
      allocate(curve_set%curve(ub)%data(lbound(data1,1):ubound(data1,1),2))
      curve_set%curve(ub)%data(:,1)=data1(:)
      curve_set%curve(ub)%data(:,2)=data2(:)
      curve_set%lb=lb
      curve_set%ub=ub

end subroutine put_data_curve_set_1d
!****************************************************
subroutine put_data_curve_set_2d(string,curve_set,data,name)
character(len=*)                           :: string
type(curve_set_type)                       :: curve_set
real(kind=rk),dimension(:,:)               :: data
character(len=*)                           :: name
integer                                    :: lb,ub
type(curve_type),dimension(:),allocatable  :: curve
logical                                    :: is_alloc
integer                                    :: ll

   lb=curve_set%lb 
   ub=curve_set%ub 

   is_alloc=allocated(curve_set%curve) 
   if(is_alloc) then
      if(ub+1<=ubound(curve_set%curve,1) ) then
         ub=ub+1
      else
         allocate(curve(lb:ub+10)) 
         do ll=lb,ub
            curve(ll)=curve_set%curve(ll)
         enddo
         call move_alloc(curve,curve_set%curve)
         ub=ub+1
      endif
   else
      call allocate_curve_set(curve_set,1)
      lb=lbound(curve_set%curve,1)
      ub=lb
   endif

      curve_set%curve(ub)%name=name
      curve_set%curve(ub)%data=data
      curve_set%lb=lb
      curve_set%ub=ub

end subroutine put_data_curve_set_2d
!****************************************************
subroutine write_data_real_1D(routine,xdata,xdescriptor,data,lb,ub,title,labels,descriptors,filename)
character(len=*)                       :: routine
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=:),allocatable           :: descr
real(kind=rk),dimension(:)             :: data
real(kind=rk),dimension(:),optional    :: xdata
integer                                :: mm,nn
integer                                :: unit
integer                                :: number_of_data
integer                                :: lb,ub
character(len=:),allocatable           :: data_filename


!      write(*,'(*(a))') 'write_data_real_1D: ',trim(filename),' will be opened'

open(newunit=unit,file=filename) 
call delete_file(filename)

!      write(unit,'(a,a)') '# written by write_data_real_1D called from ',trim(routine)

!if(present(xdescriptor) .or. present(descriptors)) then
!      write(unit,'(a)',advance='no') '# ' 
!endif
!if(present(xdescriptor)) then
!      write(unit,'(a14,a1)',advance='no')  trim(xdescriptor),' '
!endif
!if(present(descriptors)) then
!   do mm=1,size(descriptors)
!      write(unit,'(a15)',advance='no') trim(descriptors(mm))
!   enddo
!endif
!if(present(xdescriptor) .or. present(descriptors)) then
!      write(unit,'(a)',advance='yes') '' 
!endif
!
!   if(present(title))  write(unit,'(*(a,a))',advance='yes')  '# title= "',trim(title),'"'
!
!   if(present(labels)) write(unit,'(*(a,a))',advance='yes')  '# xlabel= "',trim(labels(1)),'" ylabel= "',trim(labels(2)),'"'


if(present(xdata)) then
        number_of_data=2
   write(unit,'(*(a,i0))',advance='yes')  '# lb= ',lb,' ub= ',ub,' number_of_data= ',number_of_data
   write(unit,'(*(a,i0))',advance='yes')  '# end_of_header '
   do nn=1,size(data)
      write(unit,'(2e15.7e2)',advance='yes') xdata(nn),data(nn)
   enddo
else
        number_of_data=1
   write(unit,'(*(a,i0))',advance='yes')  '# lb= ',lb,' ub= ',ub,' number_of_data= ',number_of_data
   write(unit,'(*(a,i0))',advance='yes')  '# end_of_header '
   do nn=1,size(data)
      write(unit,'(2e15.7e2)',advance='yes') data(nn)
   enddo
endif
   close(unit)
      write(*,'(*(a))') 'write_data_real_1D: ',trim(filename),' has been closed'
end subroutine write_data_real_1D
!****************************************************
subroutine write_data_real_2D(routine,xdata,xdescriptor,data,lb,ub,title,labels,descriptors,filename)
character(len=*)                       :: routine
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=:),allocatable           :: descr
real(kind=rk),dimension(:,:)           :: data
real(kind=rk),dimension(:),optional    :: xdata
integer                                :: mm,nn
integer                                :: unit
integer                                :: number_of_data
integer                                :: lb,ub
integer                                :: size2 
character(len=:),allocatable           :: data_filename


      write(*,'(*(a))') 'write_data_real_2D: ',trim(filename),' will be opened'

open(newunit=unit,file=filename) 
call delete_file(filename)

      write(unit,'(a,a)') '# written by write_data_real_2D called from ',trim(routine)

if(present(xdescriptor) .or. present(descriptors)) then
      write(unit,'(a)',advance='no') '# ' 
endif

if(present(xdescriptor)) then
      write(unit,'(a14,a1)',advance='no')  trim(xdescriptor),' '
endif

if(present(descriptors)) then
   do mm=1,size(descriptors)
      write(unit,'(a15)',advance='no') trim(descriptors(mm))
   enddo
endif
if(present(xdescriptor) .or. present(descriptors)) then
      write(unit,'(a)',advance='yes') '' 
endif

   if(present(title))  write(unit,'(*(a,a))',advance='yes')  '# title= "',trim(title),'"'

   if(present(labels)) write(unit,'(*(a,a))',advance='yes')  '# xlabel= "',trim(labels(1)),'" ylabel= "',trim(labels(2)),'"'


      size2=size(data,2)

if(present(xdata)) then
        number_of_data=1+size2
   write(unit,'(*(a,i0))',advance='yes')  '# lb= ',lb,' ub= ',ub,' number_of_data= ',number_of_data
   write(unit,'(*(a,i0))',advance='yes')  '# end_of_header '
   do nn=1,size(data,1)
      write(unit,'(*(e15.7e2))',advance='yes') xdata(nn),(data(nn,mm),mm=1,size2)
   enddo
else
        number_of_data=size2
   write(unit,'(*(a,i0))',advance='yes')  '# lb= ',lb,' ub= ',ub,' number_of_data= ',number_of_data
   write(unit,'(*(a,i0))',advance='yes')  '# end_of_header '
   do nn=1,size(data,1)
      write(unit,'(*(e15.7e2))',advance='yes') (data(nn,mm),mm=1,size2)
   enddo
endif
   close(unit)
      write(*,'(*(a))') 'write_data_real_2D: ',trim(filename),' has been closed'
end subroutine write_data_real_2D
!****************************************************
subroutine write_curve_set(routine,curve_set,title,labels,descriptors,filename)

character(len=*)                       :: routine
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=:),allocatable           :: descr
type(curve_set_type)                   :: curve_set
integer                                :: ll,mm,nn
integer                                :: unit
integer                                :: number_of_data
integer                                :: lb,ub
integer                                :: size2 
character(len=:),allocatable           :: data_filename


      write(*,'(*(a))') 'write_curve_set: ',trim(filename),' will be opened'

open(newunit=unit,file=filename) 
call delete_file(filename)

      write(unit,'(a,a)') '# written by write_curve_set called from ',trim(routine)



if(present(descriptors)) then
      write(unit,'(a)',advance='no') '# ' 
   do mm=1,size(descriptors)
      write(unit,'(a15)',advance='no') trim(descriptors(mm))
   enddo
      write(unit,'(a)',advance='yes') '' 
endif

   if(present(title))  write(unit,'(*(a,a))',advance='yes')  '# title= "',trim(title),'"'

   if(present(labels)) write(unit,'(*(a,a))',advance='yes')  '# xlabel= "',trim(labels(1)),'" ylabel= "',trim(labels(2)),'"'



   write(unit,'(*(a,i0))',advance='yes')  '# curve_set_lower_bound= ',curve_set%lb,' curve_set_upper_bound= ',curve_set%ub          
   write(unit,'(*(a,i0))',advance='yes') '# end_of_header '

   do ll=curve_set%lb,curve_set%ub
         lb=lbound(curve_set%curve(ll)%data,1)
         ub=ubound(curve_set%curve(ll)%data,1)
         number_of_data=size(curve_set%curve(ll)%data,2)
         write(unit,'(4(a,i0),1(a,a),a)',advance='yes')  '# ll= ',ll,' lb= ',lb,' ub= ',ub,' number_of_data= ',number_of_data &
                                                     & ,' property= "',trim(curve_set%curve(ll)%name),'"'
         write(unit,'(*(a,i0))',advance='yes') '# end_of_header '
      do nn=lb,ub
         write(unit,'(*(e15.7e2))',advance='yes') (curve_set%curve(ll)%data(nn,mm),mm=1,number_of_data)
      enddo
         write(unit,'(a)',advance='yes') ''
         write(unit,'(a)',advance='yes') ''
   enddo
   close(unit)
      write(*,'(*(a))') 'write_curve_set: ',trim(filename),' has been closed'

end subroutine write_curve_set
!****************************************************
subroutine read_curve_set(curve_set,title,filename)

    type(curve_set_type):: curve_set
    character(len=*)    :: title
    character(len=*)    :: filename
    character(len=256)  :: string
    integer             :: unit
    integer             :: lb,ub
    integer             :: ind  
    integer             :: ll,mm,nn  
    integer             :: number_of_data
    integer             :: number_of_curves
    real(kind=dk)       :: value
    character(len=64),dimension(:),allocatable     :: keys
    character(len=64),dimension(:),allocatable     :: param
    logical,dimension(:),allocatable               :: defined
    integer             :: status  
    integer             :: start,end_i

    write(*,'(a,a)') 'read_curve_set: opening file ',trim(filename)
open(newunit=unit,file=filename,status='old')

allocate(keys(7),param(7),defined(7))
title=''
keys(1)='title='
keys(2)='curve_set_lower_bound='
keys(3)='curve_set_upper_bound='
keys(4)='lb='
keys(5)='ub='
keys(6)='number_of_data='
keys(7)='property='

!call read_header_curve_set(keys,param,defined,unit)

if(defined(1)) title=trim(param(1))
if(defined(2)) read(param(2),*) curve_set%lb 
if(defined(3)) read(param(3),*) curve_set%ub
!!if(defined(4)) read(param(4),*) lb
!!if(defined(5)) read(param(5),*) ub
!!if(defined(6)) read(param(6),*) number_of_data
!!if(defined(7)) curve_set%curve(ll)%name=trim(param(7))

allocate(curve_set%curve(curve_set%lb:curve_set%ub))
do ll=curve_set%lb,curve_set%ub


!call read_header_curve_set(keys,param,defined,unit)

if(defined(1)) title=trim(param(1))
!!if(defined(2)) read(param(2),*) curve_set%lb 
!!if(defined(3)) read(param(3),*) curve_set%ub
if(defined(4)) read(param(4),*) lb
if(defined(5)) read(param(5),*) ub
if(defined(6)) read(param(6),*) number_of_data
if(defined(7)) curve_set%curve(ll)%name=trim(param(7))

allocate(curve_set%curve(ll)%data(lb:ub,number_of_data))

      do nn=lb,ub
         read(unit,'(a)',end_i=1000) string
         read(string,*) (curve_set%curve(ll)%data(nn,mm),mm=1,number_of_data)
      enddo

enddo

1000 continue
close(unit)

end subroutine read_curve_set
!***************************************************************************
subroutine plot_array_real(data,filename,x_axis,xdescriptor,xrange,yrange,title,labels,descriptors,load_name,append_name,no_legend)

!**************************************************
! TYPES BEING CALLED IN SEND OUT
!**************************************************
real(kind=rk),dimension(:),allocatable :: data
character(len=*)                       :: filename
!**************************************************
! LOCAL TYPES
!**************************************************
character(len=:),allocatable           :: data_filename
integer                                :: lb,ub
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=*),optional              :: load_name
character(len=*),optional              :: append_name
character(len=:),allocatable           :: descr
real(kind=rk),dimension(:),optional    :: x_axis
real(kind=rk),dimension(:),optional    :: xrange
real(kind=rk),dimension(:),optional    :: yrange
integer                                :: size_descr
integer                                :: mm,nn
integer                                :: unit,gnu_unit
logical,optional                       :: no_legend
logical                                :: loc_no_legend
character(len=:),allocatable           :: title_string
character(len=:),allocatable           :: append_string

!set xtics
!plot 'datafile' using 3:4:xticlabels(1) with linespoints

!if(present(no_legend) ) then
!        loc_no_legend=no_legend
!else
!        loc_no_legend=.false.
!endif
!if(present(descriptors)) then
!   size_descr=size(descriptors)
!endif

if(SIZE(DATA)>0) then

data_filename=trim(filename)//'.dat' 

      lb=lbound(data,1)
      ub=ubound(data,1)

call write_data('plot_array_real',x_axis,xdescriptor,data,lb,ub,title,labels,descriptors,data_filename)

call init_gnuplot_script(filename,title,gnu_unit)

!   if(present(descriptors)) then
!      descr=trim(descriptors(1))
!   else
!      descr=trim(filename)//'_'//trim(string_of(1))
!   endif
!
!   if(loc_no_legend) then
!        title_string=' notitle '
!   else
!        title_string=' title "'//descr//'" ' 
!   endif

!   if(present(append_name)) then
!        append_string=', '//trim(append_name)
!   else
!        append_string=' '                     
!   endif

   write(gnu_unit,'(*(a))',advance='yes') 'set multiplot ' 
if(present(x_axis)) then
   write(gnu_unit,'(*(a))',advance='no') 'plot "',data_filename,'" using 1:2 ', title_string,append_string
else
   write(gnu_unit,'(*(a))',advance='no') 'plot "',data_filename,'" using 0:1 ', title_string,append_string
endif
   write(gnu_unit,'(*(a))',advance='yes') ''

  if(present(load_name)) then
   write(gnu_unit,'(*(a))',advance='yes') 'load "',trim(load_name)//'.gnu_plot','"'
  endif

   write(gnu_unit,'(*(a))',advance='yes') 'unset multiplot ' 
   write(gnu_unit,'(*(a))') 'pause -1 "',trim(filename),'"'

     call write_gnuplot_trailer(gnu_unit,filename)

close(gnu_unit)

endif

end subroutine plot_array_real
!***************************************************************************
subroutine plot_array_allocatable_complex(data,circle_area,circle_scale,xrange,yrange &
                          & ,filename,title,labels,descriptors,xdescriptor,load_name &
                          & ,line_type,line_or_points,no_legend)
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),optional              :: load_name
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=:),allocatable           :: descr
complex(kind=rk),dimension(:),allocatable :: data
complex(kind=rk),dimension(:),allocatable,optional :: circle_area
real(kind=rk),optional                 :: circle_scale
real(kind=rk)                          :: loc_circle_scale
real(kind=rk),dimension(:,:),allocatable :: real_data
real(kind=rk),dimension(:),optional    :: xrange
real(kind=rk),dimension(:),optional    :: yrange
integer                                :: size_1,size_descr
integer                                :: mm,nn
integer                                :: unit,gnu_unit,gnu_unit_plot
integer                                :: unit_end
integer,optional                       :: line_type
integer                                :: loc_line_type
character(len=*),optional              :: line_or_points
character(len=20)                      :: loc_line_or_points
character(len=:),allocatable           :: mylinetype
character(len=:),allocatable           :: data_filename
integer                                :: lb,ub
logical,optional                       :: no_legend
logical                                :: loc_no_legend
character(len=:),allocatable           :: title_string
character(len=:),allocatable           :: line_type_string
character(len=:),allocatable           :: point_type_string

!set xtics
!plot 'datafile' using 3:4:xticlabels(1) with linespoints

if(present(no_legend) ) then
        loc_no_legend=no_legend
else
        loc_no_legend=.false.
endif
if(present(circle_scale) ) then
        loc_circle_scale=circle_scale
else
        loc_circle_scale=1._dk
endif

if(present(line_or_points) ) then
        loc_line_or_points=line_or_points
      point_type_string=' with '//trim(loc_line_or_points)//' ps mypointsize '
else
        loc_line_or_points='lines'
      point_type_string=' with lines '
endif

   mylinetype='mylinetype_'//trim(filename)
if(present(line_type) ) then
        loc_line_type=line_type
        line_type_string=' lt '//mylinetype
else
        loc_line_type=1
        line_type_string=' '
endif

size_1=size(data)
if(present(descriptors)) then
   size_descr=size(descriptors)
endif

if(size_1>0) then


   open(newunit=gnu_unit_plot,file=trim(filename)//'.gnu_plot') 
   call delete_file(trim(filename)//'.gnu_plot')

data_filename=trim(filename)//'.dat' 


lb=lbound(data,1)
ub=ubound(data,1)

if(present(circle_area)) then
   allocate(real_data(lb:ub,3))
   real_data(:,1)=real(data(:),kind=dk)
   real_data(:,2)=aimag(data(:))
   real_data(:,3)=abs(circle_area)
else
   allocate(real_data(lb:ub,2))
   real_data(:,1)=real(data(:),kind=dk)
   real_data(:,2)=aimag(data(:))
endif

call write_data('plot_array_allocatable_complex',xdescriptor=xdescriptor,data=real_data,lb=lb,ub=ub,title=title,labels=labels &
               & ,descriptors=descriptors,filename=data_filename)


call init_gnuplot_script(filename,title,gnu_unit)

   write(gnu_unit,'(*(a))') 'set style data lines '
   if(present(descriptors)) then
      descr=trim(descriptors(1))
   else
      !!descr=trim(filename)//'_'//trim(string_of(1))
      descr=trim(filename)
   endif

if(present(circle_area)) then
   write(gnu_unit,'(a,f0.3)') 'circle_scale=',circle_scale
   write(gnu_unit,'(*(a))') 'set style fill transparent solid my_transparency noborder'
endif


   write(gnu_unit,'(a)') 'key_x=key_x_start'
   write(gnu_unit,'(a)') 'key_y=key_y_start'
   !!write(gnu_unit,'(a,f0.3)') 'mypointsize=',0.2
   write(gnu_unit,'(*(a))',advance='yes') 'set multiplot ' 
   write(gnu_unit,'(*(a))',advance='yes') 'load "',trim(filename)//'.gnu_plot','"' 
   write(gnu_unit_plot,'(a,a,i0)') mylinetype,'=',loc_line_type
   write(gnu_unit_plot,'(*(a))') 'test_case="',trim(filename),'"' 

   write(gnu_unit_plot,'(*(a))') 'key_x=key_x+key_dx'
   write(gnu_unit_plot,'(*(a))') 'key_y=key_y+key_dy'
   write(gnu_unit_plot,'(*(a))') 'set key at character key_x,key_y'


   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

   if(present(circle_area)) then
      write(gnu_unit_plot,'(*(a))',advance='yes') 'plot ',& 
      & '"',data_filename,'" using 1:2:(circle_scale*sqrt($3)) ',title_string,' with circles'


      write(gnu_unit_plot,'(*(a))',advance='yes') 'plot "',data_filename &
                           & ,'" using 1:2 notitle ',point_type_string,line_type_string
   else
      write(gnu_unit_plot,'(*(a))',advance='yes') &
       & 'plot "',data_filename,'" using 1:2 ',title_string,point_type_string,line_type_string
   endif

  if(present(load_name)) then
   write(gnu_unit,'(*(a)     )',advance='yes') 'load "',trim(load_name)//'.gnu_plot','"'
  endif

   write(gnu_unit,'(*(a)     )',advance='yes') 'unset multiplot ' 
   write(gnu_unit,'(*(a) )') 'pause -1 "',trim(filename),'"'

     call write_gnuplot_trailer(gnu_unit,filename)

close(gnu_unit_plot) 
close(gnu_unit)

endif

end subroutine plot_array_allocatable_complex
!***************************************************************************
!***************************************************************************
subroutine plot_array_allocatable_2d(data,x_axis,xdescriptor,xrange,yrange,filename,title,labels,descriptors,no_legend)
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=:),allocatable           :: descr
real(kind=rk),dimension(:,:),allocatable :: data
real(kind=rk),dimension(:),optional    :: x_axis
real(kind=rk),dimension(:),allocatable :: numbers
real(kind=rk),dimension(:),optional    :: xrange
real(kind=rk),dimension(:),optional    :: yrange
integer                                :: size_1,size_2,size_descr
integer                                :: mm,nn
integer                                :: unit,gnu_unit
character(len=:),allocatable          :: close
integer                                :: lb,ub
character(len=:),allocatable           :: data_filename
logical,optional                       :: no_legend
logical                                :: loc_no_legend
character(len=:),allocatable           :: title_string

!set xtics
!plot 'datafile' using 3:4:xticlabels(1) with linespoints

size_1=size(data,1)
size_2=size(data,2)

!if(present(descriptors)) then
!   size_descr=size(descriptors)
!endif
!if(present(no_legend) ) then
!        loc_no_legend=no_legend
!else
!        loc_no_legend=.false.
!endif

if(size_1>0) then
data_filename=trim(filename)//'.dat' 


      lb=lbound(data,1)
      ub=ubound(data,1)

      if(present(x_axis) ) then
         call write_data('plot_array_allocatable_2d',x_axis,xdescriptor,data,lb,ub,title,labels,descriptors,data_filename)
      else
         allocate(numbers(lb:ub))
         numbers(lb:)=[(mm,mm=lb,ub)]
         descr='current_number'
         call write_data('plot_array_allocatable_2d',numbers,descr,data,lb,ub,title,labels,descriptors,data_filename)
      endif

call init_gnuplot_script(filename,title,gnu_unit)
   if(present(descriptors)) then
      descr=trim(descriptors(1))
   else
      descr=trim(filename)//'_'//trim(string_of(1))
   endif

   write(gnu_unit,'(a,a,a,a,a)',advance='yes') 'plot \'
   do mm=1,min(size_2,10)
   if(.not. present(descriptors)) then
      descr=trim(filename)//'_'//trim(string_of(mm))
   else
      descr=descriptors(mm)
   endif


   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

   close=',\'
   if(mm==min(size_2,10)) close=''
   write(gnu_unit,'(1(a,a),2(a,i0),a,a)',advance='yes') &
                           &  ' "',data_filename&
                           & ,'" using ',1 &
                           & ,':',mm+1&
                           & , title_string,close
   enddo
   write(gnu_unit,'(a,a,a)') 'pause -1 "',trim(filename),'"'

   do mm=2,size_2
   if(.not. present(descriptors)) then
      descr=trim(filename)//'_'//trim(string_of(1))
   else
      descr=descriptors(1)
   endif
   write(gnu_unit,'(a,a,a,a,a)',advance='yes') 'plot \'
   close=',\'
   write(gnu_unit,'(1(a,a),2(a,i0),a,a)',advance='yes') &
                           &  ' "',data_filename&
                           & ,'" using ',1 &
                           & ,':',1+1&
                           & , title_string,close
   if(.not. present(descriptors)) then
      descr=trim(filename)//'_'//trim(string_of(mm))
   else
      descr=descriptors(mm)
   endif
    close=''
   write(gnu_unit,'(1(a,a),2(a,i0),a,a)',advance='yes') &
                           &  ' "',data_filename&
                           & ,'" using ',1 &
                           & ,':',mm+1&
                           & , title_string,close
   write(gnu_unit,'(a,a,i0,a)') 'pause -1 "',trim(filename),mm,'"'
   enddo

   !!write(gnu_unit,'(a,a,a,i0)',advance='yes') ''
   !!write(gnu_unit,'(a,a,a)') 'pause -1 "',trim(filename),'"'

     call write_gnuplot_trailer(gnu_unit,filename)

close(gnu_unit)

endif

end subroutine plot_array_allocatable_2d
!***************************************************************************
subroutine plot_curve_set(curve_set,xrange,yrange,filename,title,labels,descriptors,line_type,line_or_points,no_legend)
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=:),allocatable          :: descr
type(curve_set_type)                   :: curve_set
real(kind=rk),dimension(:),allocatable :: numbers
real(kind=rk),dimension(:),optional    :: xrange
real(kind=rk),dimension(:),optional    :: yrange
integer                                :: size_1,size_2,size_descr
integer                                :: mm,nn
integer                                :: unit,gnu_unit,gnu_unit_plot
character(len=:),allocatable          :: close
integer                                :: lb,ub
character(len=:),allocatable           :: data_filename
integer,optional                       :: line_type
integer                                :: loc_line_type
character(len=*),optional              :: line_or_points
character(len=20)                      :: loc_line_or_points
logical,optional                       :: no_legend
logical                                :: loc_no_legend
character(len=:),allocatable           :: title_string
character(len=:),allocatable           :: line_type_string

!set xtics
!plot 'datafile' using 3:4:xticlabels(1) with linespoints



data_filename=trim(filename)//'.dat' 



if(size(curve_set%curve)==0) then
        write(*,'(*(a))') 'plot_curve_set: no data found for "',trim(filename),'"'
        goto 1000
endif
call write_curve_set('plot_curve_set',curve_set,title,labels,descriptors,data_filename)

if(present(line_or_points) ) then
        loc_line_or_points=line_or_points
else
        loc_line_or_points='lines'
endif

if(present(line_type) ) then
        loc_line_type=line_type
        line_type_string=' lt '//string_of(line_type)
else
        loc_line_type=1
        line_type_string=''
endif

if(present(no_legend) ) then
        loc_no_legend=no_legend
else
        loc_no_legend=.false.
endif

call init_gnuplot_script(filename,title,gnu_unit)

   open(newunit=gnu_unit_plot,file=trim(filename)//'.gnu_plot') 
   call delete_file(trim(filename)//'.gnu_plot')


   write(gnu_unit,'(a)') 'key_x=key_x_start'
   write(gnu_unit,'(a)') 'key_y=key_y_start'

   write(gnu_unit_plot,'(*(a))') 'key_x=key_x+key_dx'
   write(gnu_unit_plot,'(*(a))') 'key_y=key_y+key_dy'
   write(gnu_unit_plot,'(*(a))') 'set key at character key_x,key_y'

     size_2=size(curve_set%curve)
     lb=curve_set%lb
     ub=curve_set%ub

     nn=2
   write(gnu_unit_plot,'(a,a,a,a,a)',advance='yes') 'plot \'

   !!do mm=1,min(size_2,10)
   !!do mm=1,size_2
   do mm=lb,ub
   if(.not. present(descriptors)) then
      !!descr=trim(filename)//'_'//trim(string_of(mm))
      descr=trim(curve_set%curve(mm)%name) 
   else
      descr=descriptors(mm-lb+1)
   endif

   close=',\'
   if(mm==ub) close=''

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

   write(gnu_unit_plot,'(1(a,a),3(a,i0),1(a),1(a,a),1(a),1(a))',advance='yes') &
                           &  ' "',data_filename&
                           & ,'" using ',1  ,':',nn  ,' index ',mm-lb &
                           & ,title_string &
                           & ,' with ',trim(loc_line_or_points)  &
                           & ,line_type_string &
                           & ,close
   enddo
   write(gnu_unit,'(a,a,a)') 'load "',trim(filename)//'.gnu_plot','"'
   write(gnu_unit,'(a,a,a)') 'pause -1 "',trim(filename),'"'



     call write_gnuplot_trailer(gnu_unit,filename)

close(gnu_unit)
close(gnu_unit_plot)

1000 continue

end subroutine plot_curve_set
!***************************************************************************
subroutine plot_function(data,x_axis,filename,title,labels,descriptors,xdescriptor,no_legend)
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=:),allocatable           :: descr
real(kind=rk),dimension(:,:)           :: data
real(kind=rk),dimension(:),optional    :: x_axis
integer                                :: size_1,size_2,size_descr
integer                                :: mm,nn
integer                                :: unit,gnu_unit
integer                                :: lb,ub        
character(len=:),allocatable           :: data_filename
logical,optional                       :: no_legend
logical                                :: loc_no_legend
character(len=:),allocatable           :: title_string


if(present(no_legend) ) then
        loc_no_legend=no_legend
else
        loc_no_legend=.false.
endif
size_1=size(data,1)
size_2=size(data,2)
if(present(descriptors)) then
   size_descr=size(descriptors)
endif
if(present(no_legend) ) then
        loc_no_legend=no_legend
else
        loc_no_legend=.false.
endif

if(size_1*size_2>0) then
data_filename=trim(filename)//'.dat' 

      lb=lbound(data,1)
      ub=ubound(data,1)

call write_data('plot_function',x_axis,xdescriptor,data,lb,ub,title,labels,descriptors,data_filename)


call init_gnuplot_script(filename,title,gnu_unit)

if(size_2==1) then
if(present(descriptors)) then
   descr=trim(descriptors(1))
else
   descr=trim(filename)//'_'//trim(string_of(1))
endif

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

   write(gnu_unit,'(a,a,a,a,a)',advance='no') 'plot "',data_filename,'" using 0:1 ', title_string
elseif(size_2==2) then
if(present(descriptors)) then
   descr=trim(descriptors(2))
else
   descr=trim(filename)//'_'//trim(string_of(2))
endif

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

   write(gnu_unit,'(a,a,a,a,a)',advance='no') 'plot "',data_filename,'" using 1:2 ', title_string
else

   write(gnu_unit,'(a,a,a,i0)') 'plot \'
do mm=1,size_2
   if(mm>1) write(gnu_unit,'(a,a,a,i0)',advance='no') ','
if(present(descriptors)) then
   if(mm<=size_descr) then
       descr=trim(descriptors(mm))
   else
       descr=trim(descriptors(size_descr))//'_'//trim(string_of(mm-size_descr+1))
   endif
else
        descr=trim(filename)//'_'//trim(string_of(mm-size_descr+1))
endif

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

   write(gnu_unit,'(a,a,a,i0,a)',advance='no') '"',data_filename,'" using 1:',mm, title_string
   if(mm<size_2) write(gnu_unit,'(a,a,a,i0)',advance='yes') '\'
enddo
endif
   write(gnu_unit,'(a,a,a,i0)',advance='yes') ''
   write(gnu_unit,'(a,a,a)') 'pause -1 "',trim(filename),'"'

     call write_gnuplot_trailer(gnu_unit,filename)

close(gnu_unit)

endif

end subroutine plot_function
!***************************************************************************
subroutine plot_curve(data,filename,title,labels,descriptors,xdescriptor,no_legend)
character(len=*)                       :: filename
character(len=*),optional              :: title
character(len=*),dimension(:),optional :: labels
character(len=*),dimension(:),optional :: descriptors
character(len=*),optional              :: xdescriptor
character(len=:),allocatable           :: descr
real(kind=rk),dimension(:,:)           :: data
integer                                :: size_1,size_2,size_descr
integer                                :: mm,nn
integer                                :: mm1,mm2
integer                                :: unit,gnu_unit
character(len=:),allocatable           :: data_filename
logical,optional                       :: no_legend
logical                                :: loc_no_legend
character(len=:),allocatable           :: title_string


if(present(no_legend) ) then
        loc_no_legend=no_legend
else
        loc_no_legend=.false.
endif
size_1=size(data,1)
size_2=size(data,2)
if(present(descriptors)) then
   size_descr=size(descriptors)
else
   size_descr=0
endif

if(size_1*size_2>0) then
data_filename=trim(filename)//'.dat' 
open(newunit=unit,file=data_filename) 
call delete_file(data_filename)

if(present(descriptors)) then
      write(unit,'(a1)',advance='no') '#'
   do mm=1,size_descr
      write(unit,'(a15)',advance='no') trim(descriptors(mm))
   enddo
   do mm=size_descr+1,size_2
      write(unit,'(a15)',advance='no') trim(descriptors(size_descr))//'_'//trim(string_of(mm-size_descr+1))
   enddo
      write(unit,'(a)',advance='yes') '' 
endif

mm1=1; mm2=2
!!do nn=1,size_1
!!      write(unit,'(2e15.7e2)',advance='yes') data(nn,mm1),data(nn,mm2)
!!enddo
do nn=1,size_1
      write(unit,'(*(e15.7e2))',advance='yes') (data(nn,mm),mm=1,size_2)
enddo
   !!   write(unit,'(2e15.7e2)',advance='yes') 


call init_gnuplot_script(filename,title,gnu_unit)
                                                       
if(size_2==1) then
   if(present(descriptors)) then
      descr=trim(descriptors(1))
   else
      descr=trim(filename)//'_'//trim(string_of(1))
   endif

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

      write(gnu_unit,'(a,a,a,a,a)',advance='no') 'plot "',data_filename,'" using 0:1 ', title_string
elseif(size_2==2) then
   if(present(descriptors)) then
      descr=trim(descriptors(2))
   else
      descr=trim(filename)//'_'//trim(string_of(2))
   endif

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

      write(gnu_unit,'(*(a))',advance='no') 'plot "',data_filename,'" using 1:2 ', title_string
   else

      write(gnu_unit,'(a,a,a,i0)') 'plot \'
   do mm=1,size_2
      if(mm>1) write(gnu_unit,'(a,a,a,i0)',advance='no') ','
   if(present(descriptors)) then
      if(mm<=size_descr) then
          descr=trim(descriptors(mm))
      else
          descr=trim(descriptors(size_descr))//'_'//trim(string_of(mm-size_descr+1))
      endif
   else
          descr=trim(filename)//'_'//trim(string_of(mm-size_descr+1))
   endif

   if(loc_no_legend) then
        title_string=' notitle '
   else
        title_string=' title "'//descr//'" ' 
   endif

      write(gnu_unit,'(a,a,a,i0,a)',advance='no') '"',data_filename,'" using 1:',mm,title_string
      if(mm<size_2) write(gnu_unit,'(a,a,a,i0)',advance='yes') '\'
   enddo
endif
   write(gnu_unit,'(a,a,a,i0)',advance='yes') ''
   write(gnu_unit,'(a,a,a)') 'pause -1 "',trim(filename),'"'

     call write_gnuplot_trailer(gnu_unit,filename)

close(gnu_unit)
close(unit)

endif

end subroutine plot_curve
!***************************************************************************
subroutine init_gnuplot_script(filename,title,unit)
character(len=*)                       :: filename
character(len=*),optional              :: title
integer                                :: unit
   open(newunit=unit,file=trim(filename)//'.gnu_script') 
  !! call delete_file(trim(filename)//'.gnu_script')
  ! call append_plot(trim(filename)//'.gnu_script')

   write(unit,'(a)') 'set style data lines'

   write(unit,'(a)') 'unset title '
   write(unit,'(a)') 'unset xlabel '
   write(unit,'(a)') 'unset ylabel '
   write(unit,'(a)') 'load "parameter.gnu_script"'
   write(unit,'(a)') 'set autoscale '
   write(unit,'(a)') 'set xrange[xmin:xmax]'
   write(unit,'(a)') 'set yrange[ymin:ymax]'
   write(unit,'(a)') 'set zrange[zmin:zmax]'
!   write(unit,'(3a)') '#set log x '
!   write(unit,'(3a)') '#set log y '
!   if(present(title)) write(unit,'(3a)') 'if(plot_title==1) set title test_case." ',trim(title),'"'
   write(unit,'(3a)') 'set xlabel x_label'
   write(unit,'(3a)') 'set ylabel y_label'
   write(unit,'(3a)') 'unset arrow '
   write(unit,'(3a)') 'unset label '

end subroutine init_gnuplot_script
!***************************************************************************
subroutine plot_incrementing_labels(node_coord,label_number,filename,gnu_script_unit)

real(kind=rk)          :: node_coord(:,:)
integer,dimension(:)   :: label_number         
integer                :: unset_label_unit
integer                :: last
integer                :: gnu_script_unit
integer                :: nn
integer                :: max_size
character(len=2)       :: memo='nd'
character(len=*)       :: filename
integer,allocatable,dimension(:) :: color
integer,allocatable,dimension(:) :: enumeration

   last=size(label_number)
   max_size=max(last,maxval(label_number))
   allocate(color(max_size))
   color=0

   allocate(enumeration(max_size))
   enumeration=0
  do nn=1,last 
   !enumeration(label_number(nn))=nn
   enumeration(nn)=nn
  enddo

    call plot_node_labels(node_coord,filename,label_number&
                    & ,enumeration=enumeration,memo=memo,color=color,gnu_script_unit=gnu_script_unit)

     unset_label_unit=50
     !!open(unset_label_unit,file='unset_labels.node_unlabel')
     open(unset_label_unit,file=trim(filename)//'.node_unlabel')
            call delete_file(trim(filename)//'.node_unlabel')


     write(unset_label_unit,'(a,a,a)') 'load "',trim(filename),'.node_label"'      
  do nn=1,last 
     !!write(unset_label_unit,'(a,i0)') 'unset label !',label_number(nn)
     write(unset_label_unit,'(a,i0,a)') 'set label ',label_number(nn),' front tc lt -1'
     write(unset_label_unit,'(a,i0)') 'replot '
     write(unset_label_unit,'(a,a,i0,a)') 'pause -1 "',trim(memo),(nn),'"'
  enddo
     close(unset_label_unit)

end subroutine plot_incrementing_labels
!***************************************************************************
subroutine plot_node_labels(node_coord,filename,index,color,enumeration,memo,implicit,gnu_script_unit,extension)
!        memo: characters the labels are starting with
!       index: the nodes which are labelled if implicit==0
!       index: the nodes for which index(nodes)==implicit 
! enumeration: a new enumeration for the labels
!
real(kind=rk)              :: node_coord(:,:)
integer,optional           :: index(:)
integer,optional           :: color(:)
integer,optional           :: enumeration(:)
character(len=*)           :: filename
integer                    :: kk1,kk2,kk3
integer                    :: node,new_node
integer                    :: gnu_label
real(kind=rk)              :: x1,y1,z1                 
integer                    :: max_nodes
character(len=*),optional  :: memo
character(len=8)           :: loc_memo
integer,optional           :: implicit
logical                    :: loc_implicit
integer,optional           :: gnu_script_unit
character(len=*),optional  :: extension
logical                    :: found
integer                    :: num_node
logical                    :: ball=.false.    ! Nummern an den Knoten
!!logical                    :: ball=.true.    ! kleine Kugeln an den Knoten
character(len=1),parameter :: ball_sign="l"
character(len=*),parameter :: ball_font='Wingdings,14'


            gnu_label=20
            if(present(extension)) then
            open(gnu_label,file=trim(filename)//'.'//trim(extension))
          !!  write(0,*) 'node labels will be written to "',trim(filename)//'.'//trim(extension),'"'
            call delete_file(trim(filename)//'.'//trim(extension))
            if(present(gnu_script_unit)) then
               write(gnu_script_unit,'(5a)') 'if( plot_'//trim(extension)//'==1) load "',trim(filename)//'.'//trim(extension),'"'
            endif
            else
            open(gnu_label,file=trim(filename)//'.node_label')
          !!  write(0,*) 'node labels will be written to "',trim(filename)//'.node_label','"'
            call delete_file(trim(filename)//'.node_label')
            if(present(gnu_script_unit)) then
               write(gnu_script_unit,'(5a)') 'if( plot_node_label==1) load "',trim(filename)//'.node_label','"'
            endif
            endif

            max_nodes=size(node_coord,1)

            found=.false.
     if(present(index)) then
            if(size(index)==0) then
                   goto 1000 
            endif
           !! call print('plot_node_labels: index for '//trim(filename),index)
    endif

            !!write(*,'(a,i0)') 'The extent of these nodes is ',max_nodes

            kk1=1
            kk2=2
            kk3=3
            if(size(node_coord,2)<3 ) then
            kk3=2
            endif

                loc_implicit=.false.
        if(present(implicit)) then
                loc_implicit=implicit/=0
        endif

                loc_memo='n'
        if(present(color)) then
                loc_memo='c'
        endif
        if(present(memo)) then
                loc_memo=memo
        endif

!        Syntax of gnuplot set label:
!        set label {<tag>} {"<label text>"} {at <position>}
!        {left | center | right}
!        {norotate | rotate {by <degrees>}}
!        {font "<name>{,<size>}"}
!        {noenhanced}
!        {front | back}
!        {textcolor <colorspec>}
!        {point <pointstyle> | nopoint}
!        {offset <offset>}

        ! the node names
index_if:     if(present(index)) then
         !!    write(0,*) 'index'
             write(*,*) 'plot_node_labels: index'
implicit_if:      if(.not. loc_implicit) then
         !!    write(0,*) 'not implicit'
             write(*,*) 'plot_node_labels: not implicit '
        do node=1,size(index)   
            if(index(node) <=0) then
                 !!write(*,'(a,i0,3a)') 'plot_node_labels(a): warning, index is 0 for node=',node,' and "',trim(filename),'"'
                 cycle
            endif

            new_node=index(node)

            x1=node_coord(new_node,kk1)
            y1=node_coord(new_node,kk2)
            z1=node_coord(new_node,kk3)

            num_node=new_node
            !!num_node=color(new_node)
            if(present(enumeration)) num_node=enumeration(new_node)
            !!if(present(enumeration)) num_node=enumeration(node)

            if(.not.ball) then
               write(gnu_label,'(2(a,i0),a)',advance='no') 'set label ',new_node,' " '//trim(loc_memo),num_node,'"'
               write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
            else
               write(gnu_label,'(1(a,i0),a,a,a)',advance='no') 'set label ',new_node,' "',ball_sign,'"'
               write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
               write(gnu_label,'(a,a,a)',advance='no') '  center font "',ball_font,'" '
            endif
            if(present(color)) then
               write(gnu_label,'(a,i0)',advance='yes') ' tc lt ',color(new_node)
            else
               write(gnu_label,'(a,i0)',advance='yes') ''
            endif
            found=.true.

        enddo
      else implicit_if 
           !!  write(0,*) 'implicit'
             write(*,*) 'plot_node_labels: implicit'
           write(*,'(*(a,i0))') 'implicit ',implicit
          !! call print('plot_node_labels: index',index)
        do node=1,max_nodes   
            if(index(node) <=0) then
                 !!write(*,'(a,i0,3a)') 'plot_node_labels(b): warning, index is 0 for node=',node,' and "',trim(filename),'"'
                 cycle
            endif
          if(index(node)==implicit) then
            x1=node_coord(node,kk1)
            y1=node_coord(node,kk2)
            z1=node_coord(node,kk3)
            if(present(color)) then
               if(.not.ball) then
                  write(gnu_label,'(2(a,i0),a)',advance='no') 'set label ',node,' " '//trim(loc_memo),color(node) ,'"'
                  write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
               else
                  write(gnu_label,'(1(a,i0),a,a,a)',advance='no') 'set label ',node,' "',ball_sign,'"'
                  write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
                  write(gnu_label,'(a,a,a)',advance='no') '  center font "',ball_font,'" '
               endif
               write(gnu_label,'(a,i0)',advance='yes') ' tc lt ',color(node)
            else
               write(gnu_label,'(2(a,i0),a)',advance='no') 'set label ',node,' " '//trim(loc_memo),node ,'"'
               write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
               write(gnu_label,'(a,i0)',advance='yes') ''
            endif
            found=.true.
          endif
        enddo
      endif implicit_if 
     else index_if
          !!   write(0,*) 'not index'
        do node=1,max_nodes   
            x1=node_coord(node,kk1)
            y1=node_coord(node,kk2)
            z1=node_coord(node,kk3)
            if(present(color)) then
               if(.not.ball) then
                  write(gnu_label,'(2(a,i0),a)',advance='no') 'set label ',node,' "'//trim(loc_memo),color(node) ,'"'
                  write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
               else
                  write(gnu_label,'(1(a,i0),a,a,a)',advance='no') 'set label ',node,' "',ball_sign,'"'
                  write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
                  write(gnu_label,'(a,a,a)',advance='no') '  center font "',ball_font,'" '
               endif
               write(gnu_label,'(a,i0)',advance='yes') ' tc lt ',color(node)

            else
               write(gnu_label,'(2(a,i0),a)',advance='no') 'set label ',node,' "'//trim(loc_memo),node,'"'
               write(gnu_label,'(3(a,e10.3))',advance='no') ' at ',x1,',',y1,',',z1
               write(gnu_label,'(a,i0)',advance='yes') ''
            endif
            found=.true.
        enddo
     endif index_if

     1000 continue
        if(.not.found) then
            !!    write(0,'(3a)') 'warning, no nodes found in plot_node_labels; file ',trim(filename)//'.node_label',' will be void'
               write(gnu_label,'(a)',advance='yes') ''
        endif

            close(gnu_label)

end subroutine plot_node_labels
!***************************************************************************
subroutine plot_graph(node_coord,edges,filename)
real(kind=rk)              :: node_coord(:,:)
integer,dimension(:,:)     :: edges
character(len=*)           :: filename
integer                    :: edge   
integer                    :: gnu_unit
integer                                :: nn
integer                    :: max_edges,max_nodes
integer                    :: jmax

            open(newunit=gnu_unit,file=trim(filename)//'.gnu_dat')
            call delete_file(trim(filename)//'.gnu_dat')

            max_edges=size(edges,1)
            write(*,'(a,i0,a,a,a)') 'plot_graph: we found ',max_edges,' edges for "',trim(filename),'"'

            max_nodes=size(node_coord,1)
            jmax=size(node_coord,2)

        do edge=1,max_edges

            ! all edges
            do nn=1,jmax
               write(gnu_unit,'(e15.7e2)',advance='no') node_coord(edges(edge,1),nn)
            enddo
            write(gnu_unit,'(a)') ''
            do nn=1,jmax
               write(gnu_unit,'(e15.7e2)',advance='no') node_coord(edges(edge,2),nn)
            enddo
            !hier dicke logar rausschreiben
            write(gnu_unit,'(a)') ''
            write(gnu_unit,'(a)') ''
            write(gnu_unit,'(a)') ''

        enddo

        if(max_edges==0) then
                write(0,*) 'warning, no edges; artificial coordinates in ',trim(filename)//'.gnu_dat'
                ! for the case that nothing is present we write a single point
               write(gnu_unit,'(e15.7e2)',advance='yes') node_coord(1,1)
               write(gnu_unit,'(e15.7e2)',advance='yes') node_coord(1,1)
        endif

            close(gnu_unit)

end subroutine plot_graph
!***************************************************************************
subroutine plot_edges(string,nodes,edges,edge_value,region,filename,line_type,reverse_arrows,append,extra)
character(len=*)                   :: string
real(kind=rk)                      :: nodes(:,:)
character(len=*)                   :: filename
logical,optional                   :: append
logical                            :: loc_append
logical,optional                   :: extra
logical                            :: loc_extra
integer,dimension(:,:)             :: edges
real(kind=rk),dimension(:)         :: edge_value
integer,dimension(:),optional      :: region
integer,dimension(:),allocatable,optional      :: line_type
logical,optional                   :: reverse_arrows
integer                            :: undir_arrows_unit 
integer                            :: dir_arrows_unit 
integer,dimension(:),allocatable   :: sort_index
integer,dimension(:),allocatable   :: directed
integer,dimension(:),allocatable   :: undirected
integer                            :: nn_edge 
integer                            :: max_temp
integer                            :: max_num
integer,dimension(:),allocatable   :: key
integer                            :: i1,i2
real(kind=rk)                      :: mue,vue       
integer                            :: mm1,mm2 
integer                            :: row,ind
integer                            :: end_directed 
integer                            :: end_undirected 
integer                            :: n1,n2
real(kind=rk)                      :: x1,y1,z1,x2,y2,z2
real(kind=rk)                      :: xlabel,ylabel,diff_label
logical                            :: has_directed
logical                            :: has_directed_small
logical                            :: has_undirected
logical                            :: has_undirected_small
integer,allocatable,dimension(:)   :: all_lt
integer                            :: ii
integer                            :: mm,ll
integer                            :: min_lt,max_lt
integer                            :: nmax
character(len=8)                   :: arrow_dir
!!character(len=8)                   :: arrow_dir='backhead'
!!character(len=8)                   :: arrow_dir='head'
!!character(len=8)                   :: arrow_dir='heads'
!!character(len=8)                   :: arrow_dir='nohead'
!real(kind=rk)                      :: arrow_head_size=0.15

            loc_extra=.false.
            if(present(extra)) loc_extra=extra
            loc_append=.false.
            if(present(append)) loc_append=append

          !!  write(*,*) 'will open file ',trim(filename),' at ',trim(string)
            if(loc_append) then
               open(newunit=undir_arrows_unit,file=trim(filename)//'.segm',access='append')
               open(newunit=dir_arrows_unit,file=trim(filename)//'.arrows',access='append')
            else
               if(loc_extra) then
                  open(newunit=undir_arrows_unit,file=trim(filename)//'.extra_segm')
                  open(newunit=dir_arrows_unit,file=trim(filename)//'.extra_arrows')
                  call delete_file(trim(filename)//'.extra_segm')
                  call delete_file(trim(filename)//'.extra_arrows')
               else
                  open(newunit=undir_arrows_unit,file=trim(filename)//'.segm')
                  open(newunit=dir_arrows_unit,file=trim(filename)//'.arrows')
                  call delete_file(trim(filename)//'.segm')
                  call delete_file(trim(filename)//'.arrows')
               endif
            endif


          !!  call print('edges for '//trim(filename),edges)

            if(present(line_type)) then
                 min_lt=minval(line_type)
                 min_lt=1
                 max_lt=maxval(line_type)
                 allocate(all_lt(min_lt:max_lt+1))
                 all_lt=[(mm,mm=min_lt,max_lt)]
            endif

                      arrow_dir='head'
            if(present(reverse_arrows)) then
                    if(reverse_arrows) then
                      arrow_dir='backhead'
                    endif
            endif

            nmax=size(nodes,2)
          !!  write(*,*) '„„„ ',trim(filename),' ',nmax
      nn_edge=size(edges,1)

    !!  write(0,'(a,i0,a,a)') ' will write ',nn_edge,' edges to file ',trim(filename)//'.arrows'
    !!  write(*,'(a,i0,a,a)') ' will write ',nn_edge,' edges to file ',trim(filename)//'.arrows'

      max_temp=nn_edge
      allocate(key(max_temp))

      max_num=maxval(edges(:,1))+1

      do mm=1,max_temp
            row=edges(mm,1); ind=edges(mm,2)   
            i1=min(row,ind); i2=max(row,ind)
            key(mm)=i1*max_num+i2-1  ! is unique for both edges connecting the node row and the node ind
      enddo

      call sort_up(key,sort_index)
      if(max_temp>size(sort_index)) then
              write(*,*) 'plot_edges: sort_index too small'
              stop 'plot_edges: sort_index too small'
      endif

      allocate(undirected(max_temp))
      allocate(directed(max_temp))
      !!call print('sorted edges for '//trim(filename),edges(sort_index(1:max_temp),:))

       ! the edges are divided in those appearing one time (directed) and those appearing two times (undirected)
         ll=1
         mm1=0
         mm2=0
      do 
         !!write(*,*) ll,max_temp
         if(ll  >  max_temp) exit
         if(ll  == max_temp) then    ! must be a single edge
                 mm2=mm2+1
                 directed(mm2)=sort_index(ll)
                 exit
         endif

         if(key(sort_index(ll))==key(sort_index(ll+1)) ) then
            if(ll<max_temp-1) then
            if(key(sort_index(ll))==key(sort_index(ll+2)) ) then
                    write(0,*) 'warning: index appears three times for ',trim(filename)
                    stop 'plot_edges: index appears three times'
            endif
            endif
         !!write(*,'(a,i0,x,i0,x,i0)') 'a:',ll,sort_index(ll),key(sort_index(ll))
                 mm1=mm1+1
                 ! waehle die Kante mit dem kleineren Wert
                 if(abs(edge_value(sort_index(ll))) <=abs(edge_value(sort_index(ll+1))) ) then
                   undirected(mm1)=sort_index(ll)
                 else
                   undirected(mm1)=sort_index(ll+1)
                 endif
                 ll=ll+2
         else
         !!write(*,'(a,i0,x,i0,x,i0)') 'b:',ll,sort_index(ll),key(sort_index(ll))
                 mm2=mm2+1
                 directed(mm2)=sort_index(ll)
                 ll=ll+1
         endif

      enddo

      if(2*mm1+mm2 /= max_temp) then
              write(0,'(3(a,i0))') 'mm1=',mm1,' mm2=',mm2,' max_temp=',max_temp
              stop 'plot_edges: error, mm1+mm2 /= max_temp'
      endif


      end_undirected=mm1
      end_directed=mm2
   !!   call print('undirected edges',edges(undirected(1:end_undirected),:))

    !!  if(present(region)) call print(' plot_edges: region',region)

 !!   write(*,'(2(a,i0),a,a)') 'we have ',end_directed,' directed edges and ',end_undirected &
 !!                            & ,' undirected edges for ',trim(filename)//'.arrows'
 !!   write(0,'(2(a,i0),a,a)') 'we have ',end_directed,' directed edges and ',end_undirected &
 !!                            & ,' undirected edges for ',trim(filename)//'.arrows'

has_directed=.false.
has_directed_small=.false.
has_undirected=.false.
has_undirected_small=.false.

            xlabel=maxval(nodes(:,1))
            ylabel=maxval(nodes(:,2))
            diff_label=0.05*(ylabel-minval(nodes(:,2)))
            xlabel=xlabel+diff_label

                 mue=0.97
                 vue=0.03

                 mue=0.99
                 vue=0.01  

                 mue=1.00
                 vue=0.00  
       do mm=1,end_undirected

          ll=undirected(mm)

    !!      write(*,'(2(a,i0))') 'undirected(',mm,')=',undirected(mm)
         if(present(region)) then
            if(region(edges(ll,1)) == 0  ) then
                    i1=2; i2=1
            elseif(region(edges(ll,2)) == 0  ) then
                    i1=1; i2=2
            elseif(region(edges(ll,1)) <  region(edges(ll,2)) ) then
                    i1=1; i2=2
            else
                    i1=2; i2=1
            endif
         else
         !!        write(0,*) 'no unidirected edge direction information given for ',trim(filename) 
         !!        stop 'no unidirected egde direction information'
            if(edges(ll,1) >= edges(ll,2) ) then
                    i1=1; i2=2
            else
                    i1=2; i2=1
            endif
         endif

!         Syntax of gnuplot set arrow:
!         set arrow {<tag>} {from <position>} {to|rto <position>}
!         { {arrowstyle | as <arrow_style>}
!            | { {nohead | head | backhead | heads}
!                {size <length>,<angle>{,<backangle>}}
!                {filled | empty | nofilled}
!                {front | back}
!                { {linestyle | ls <line_style>}
!                  | {linetype | lt <line_type>}
!                    {linewidth (ll,1) >= edges(ll,2)<line_width} 
!                } 
!              } 
!         }

         if(abs(edge_value(ll))<=skip_edges_threshould) then
         !!   write(*,'(a,i0,a,i0)') '# small edges_value for undirected edge ',edges(ll,i1),'->',edges(ll,i2)
            write(undir_arrows_unit,'(a,i0,a,i0)') '# small edges_value for undirected edge ',edges(ll,i1),'->',edges(ll,i2)
            !!write(undir_arrows_unit,'(a)',advance='no') '# '
         else
         !!   write(*,'(a,i0,a,i0)') '# undirected edge ',edges(ll,i1),'<->',edges(ll,i2)
            write(undir_arrows_unit,'(a,i0,a,i0)') '# undirected edge ',edges(ll,i1),'<->',edges(ll,i2)
            !!write(undir_arrows_unit,'(a)') '# edges_value for undirected edge'
         endif
            write(undir_arrows_unit,'(a    )',advance='no') 'set arrow from '                 
            if( nmax >= 3 )  then
            n1=edges(ll,i1); n2=edges(ll,i2)
            x1=nodes(n1,1); x2=nodes(n2,1)
            y1=nodes(n1,2); y2=nodes(n2,2)
            z1=nodes(n1,3); z2=nodes(n2,3)

            write(undir_arrows_unit,'(e12.5)',advance='no') x1+vue*(x2-x1)
            write(undir_arrows_unit,'(a    )',advance='no') ', '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') y1+vue*(y2-y1)
            write(undir_arrows_unit,'(a    )',advance='no') ', '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') z1+vue*(z2-z1)
            write(undir_arrows_unit,'(a    )',advance='no') ' to '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') x1+mue*(x2-x1)
            write(undir_arrows_unit,'(a    )',advance='no') ', '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') y1+mue*(y2-y1)
            write(undir_arrows_unit,'(a    )',advance='no') ', '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') z1+mue*(z2-z1)
            else
            n1=edges(ll,i1); n2=edges(ll,i2)
            x1=nodes(n1,1); x2=nodes(n2,1)
            y1=nodes(n1,2); y2=nodes(n2,2)

            write(undir_arrows_unit,'(e12.5)',advance='no') x1+vue*(x2-x1)
            write(undir_arrows_unit,'(a    )',advance='no') ', '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') y1+vue*(y2-y1)
            write(undir_arrows_unit,'(a    )',advance='no') ' to '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') x1+mue*(x2-x1)
            write(undir_arrows_unit,'(a    )',advance='no') ', '                 
            write(undir_arrows_unit,'(e12.5)',advance='no') y1+mue*(y2-y1)
            endif

         if(abs(edge_value(ll))<=skip_edges_threshould) then
          !!  write(*,'(a,i0,a,e10.2)') ' nohead lt ',lt_undirected_small,' lw ',linewidth_scale*edge_value(ll)
            write(undir_arrows_unit,'(a,i0,a,e10.2)') ' nohead lt ',lt_undirected_small,' lw ',linewidth_scale*edge_value(ll)
            has_undirected_small=.true.
         else
            has_undirected=.true.
        !!    write(*,'(a,i0,a,e10.2)') ' nohead lt ',lt_undirected,' lw ',linewidth_scale*edge_value(ll)
            write(undir_arrows_unit,'(a,i0,a,e10.2)') ' nohead lt ',lt_undirected,' lw ',linewidth_scale*edge_value(ll)
         endif
       enddo

 if(has_undirected) then
 write(undir_arrows_unit,'(2(a,e10.2),a,i0)')&
          &  ' set label "undirected"       at ',xlabel,',',ylabel-0*diff_label,' left tc lt ',lt_undirected
 endif
 if(has_undirected_small) then
 write(undir_arrows_unit,'(2(a,e10.2),a,i0)')&
          &  ' set label "undirected_small" at ',xlabel,',',ylabel-1*diff_label,' left tc lt ',lt_undirected_small
 endif
            close(undir_arrows_unit)

                 mue=0.97
                 vue=0.03

                 mue=0.97
                 vue=0.10  ! head

                 mue=0.99
                 vue=0.01  ! head
       do mm=1,end_directed

          ll=directed(mm)

         if(present(region)) then
            if(region(edges(ll,1)) == 0  ) then
                    i1=2; i2=1
            elseif(region(edges(ll,2)) == 0  ) then
                    i1=1; i2=2
            elseif(region(edges(ll,1)) <  region(edges(ll,2)) ) then
                    i1=1; i2=2
            else
                    i1=2; i2=1
            endif
         else
                 !!write(0,*) 'no unidirected edge direction information given for ',trim(filename) 
               !!  stop 'no unidirected egde direction information'
        !!    if(edges(ll,1) >= edges(ll,2) ) then
        !!            i1=1; i2=2
        !!    else
        !!            i1=2; i2=1
        !!    endif
                    i1=1; i2=2
         endif

!!         Syntax of gnuplot set arrow:
!!         set arrow {<tag>} {from <position>} {to|rto <position>}
!!         { {arrowstyle | as <arrow_style>}
!!            | { {nohead | head | backhead | heads}
!!                {size <length>,<angle>{,<backangle>}}
!!                {filled | empty | nofilled}
!!                {front | back}
!!                { {linestyle | ls <line_style>}
!!                  | {linetype | lt <line_type>}
!!                    {linewidth | lw <line_width} 
!!                } 
!!              } 
!!         }

         if(abs(edge_value(ll))<=skip_edges_threshould) then
            write(dir_arrows_unit,'(a,i0,a,i0)') '# small edges_value for directed edge ',edges(ll,i1),'->',edges(ll,i2)
            !!write(dir_arrows_unit,'(a)',advance='no') '# '
         else
            write(dir_arrows_unit,'(a,i0,a,i0)') '# directed edge ',edges(ll,i1),'->',edges(ll,i2)
            !!write(dir_arrows_unit,'(a)') '# edges_value for directed edge'
         endif
            write(dir_arrows_unit,'(a    )',advance='no') 'set arrow from '                 

            if( nmax >= 3 )  then
            n1=edges(ll,i1); n2=edges(ll,i2)
            x1=nodes(n1,1); x2=nodes(n2,1)
            y1=nodes(n1,2); y2=nodes(n2,2)
            z1=nodes(n1,3); z2=nodes(n2,3)

            write(dir_arrows_unit,'(e12.5)',advance='no') x1+vue*(x2-x1)
            write(dir_arrows_unit,'(a    )',advance='no') ', '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') y1+vue*(y2-y1)
            write(dir_arrows_unit,'(a    )',advance='no') ', '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') z1+vue*(z2-z1)
            write(dir_arrows_unit,'(a    )',advance='no') ' to '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') x1+mue*(x2-x1)
            write(dir_arrows_unit,'(a    )',advance='no') ', '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') y1+mue*(y2-y1)
            write(dir_arrows_unit,'(a    )',advance='no') ', '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') z1+mue*(z2-z1)
            else
            n1=edges(ll,i1); n2=edges(ll,i2)
            x1=nodes(n1,1); x2=nodes(n2,1)
            y1=nodes(n1,2); y2=nodes(n2,2)

            write(dir_arrows_unit,'(e12.5)',advance='no') x1+vue*(x2-x1)
            write(dir_arrows_unit,'(a    )',advance='no') ', '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') y1+vue*(y2-y1)
            write(dir_arrows_unit,'(a    )',advance='no') ' to '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') x1+mue*(x2-x1)
            write(dir_arrows_unit,'(a    )',advance='no') ', '                 
            write(dir_arrows_unit,'(e12.5)',advance='no') y1+mue*(y2-y1)
            endif

         if(present(line_type)) then
         if(abs(edge_value(ll))<=skip_edges_threshould) then
            write(dir_arrows_unit,'(2a)',advance='no') ' ',trim(arrow_dir)
            !!write(dir_arrows_unit,'(a,e10.2)',advance='no') ' size ',arrow_head_size
            write(dir_arrows_unit,'(a,i0,a,e10.2)',advance='yes') &
                                  & ' lt ',line_type(edges(ll,1)),' lw ',linewidth_scale*edge_value(ll)
            has_directed_small=.true.
         else
            !!write(dir_arrows_unit,'(2a,a,i0,a,e10.2)') ' ',trim(arrow_dir),' lt ',line_type(edges(ll,1)),' lw ',linewidth_scale*edge_value(ll)
            write(dir_arrows_unit,'(2a)',advance='no') ' ',trim(arrow_dir)
            !!write(dir_arrows_unit,'(a,e10.2)',advance='no') ' size ',arrow_head_size
            write(dir_arrows_unit,'(a,i0,a,e10.2)',advance='yes') &
                           &  ' lt ',line_type(edges(ll,1)),' lw ',linewidth_scale*edge_value(ll)
            has_directed=.true.
         endif
         else
         if(abs(edge_value(ll))<=skip_edges_threshould) then
            !!write(dir_arrows_unit,'(2a,a,i0,a,e10.2)') ' ',trim(arrow_dir),' lt ',lt_directed_small,' lw ',linewidth_scale*edge_value(ll)
            write(dir_arrows_unit,'(2a)',advance='no') ' ',trim(arrow_dir)
            !!write(dir_arrows_unit,'(a,e10.2)',advance='no') ' size ',arrow_head_size
          write(dir_arrows_unit,'(a,i0,a,e10.2)',advance='yes') ' lt ',lt_directed_small,' lw ',linewidth_scale*edge_value(ll)
            has_directed_small=.true.
         else
            !!write(dir_arrows_unit,'(2a,a,i0,a,e10.2)') ' ',trim(arrow_dir),' lt ',lt_directed,' lw ',linewidth_scale*edge_value(ll)
            write(dir_arrows_unit,'(2a)',advance='no') ' ',trim(arrow_dir)
            !!write(dir_arrows_unit,'(a,e10.2)',advance='no') ' size ',arrow_head_size
            write(dir_arrows_unit,'(a,i0,a,e10.2)',advance='yes')&
                                 &  ' lt ',lt_directed,' lw ',linewidth_scale*edge_value(ll)
            has_directed=.true.
         endif
         endif
       enddo

         if(present(line_type)) then
!! if(has_directed) then
!! write(dir_arrows_unit,'(2(a,e10.2),a,i0)') ' set label "directed"         at ',xlabel,',',ylabel-2*diff_label,' left tc lt ',line_type(edges(ll,1))
!! endif
!! if(has_directed_small) then
!! write(dir_arrows_unit,'(2(a,e10.2),a,i0)') ' set label "directed_small"   at ',xlabel,',',ylabel-3*diff_label,' left tc lt ',line_type(edges(ll,1))
!! endif
!write(*,*) 'min_lt=',min_lt,' max_lt=',max_lt
do ii=min_lt,max_lt
!!write(*,*) ii,ubound(all_lt),lbound(all_lt)
 write(dir_arrows_unit,'(a,i0,2(a,e10.2),a,i0)') &
                       & ' set label "',all_lt(ii),'"   at ',xlabel,',',ylabel-(ii-min_lt)*diff_label,' left tc lt ',all_lt(ii)
 enddo
         else
 if(has_directed) then
 write(dir_arrows_unit,'(2(a,e10.2),a,i0)')&
                      &  ' set label "directed"         at ',xlabel,',',ylabel-2*diff_label,' left tc lt ',lt_directed
 endif
 if(has_directed_small) then
 write(dir_arrows_unit,'(2(a,e10.2),a,i0)')&
                      &  ' set label "directed_small"   at ',xlabel,',',ylabel-3*diff_label,' left tc lt ',lt_directed_small
 endif
         endif

            close(dir_arrows_unit)

end subroutine plot_edges
!*************************************************************************************
subroutine determine_multiple_edges(filename,edges,group_start,group_index)
integer,dimension(:,:)             :: edges
integer,dimension(:),allocatable   :: key
integer,dimension(:),allocatable   :: sort_index
integer,dimension(:),allocatable,intent(out)   :: group_start
!!integer,dimension(:),allocatable,intent(out)   :: group_end
integer,dimension(:),allocatable,intent(out)   :: group_index
integer,dimension(:),allocatable   :: next
integer                            :: max_temp
integer                            :: max_num
integer                            :: row,ind
integer                            :: i1,i2   
integer                            :: kk,ll,mm   
integer                            :: max_next
integer                            :: group_size,max_group_size
character(len=*)                  :: filename



      max_temp=size(edges,1)
      if(max_temp==0 ) then
         max_group_size=0
         allocate(group_start(max_group_size+1))
         allocate(group_index(max_temp))
         goto 1111
      endif

      allocate(key(max_temp))

      max_num=maxval(edges(:,1))+1
    !!  call print('determine_multiple_edges: edges for '//trim(filename),edges)

      ! key independent on the direction of the edge
      do mm=1,max_temp
            row=edges(mm,1); ind=edges(mm,2)   
            i1=min(row,ind); i2=max(row,ind)
            key(mm)=i1*max_num+i2-1  ! is unique for both edges connecting the node row and the node ind
      enddo

      call sort_up(key,sort_index)
      if(max_temp>size(sort_index)) then
              write(*,*) 'plot_edges: sort_index too small'
              stop 'determine_multiple_edges: sort_index too small'
      endif

!!      call print('determine_multiple_edges: sorted edges for '//trim(filename),edges(sort_index(1:max_temp),:))


      allocate(next(max_temp+1))
              kk=1
              next(1)=1
      do ll=2,max_temp 
         if(key(sort_index(ll))/=key(sort_index(ll-1)) ) then
              kk=kk+1
              next(kk)=ll ! next(kk) points to the begin of the kk-th constant section
         endif
      enddo
       max_next=kk
       next(max_next+1)=max_temp+1 ! points to the first index after all actual indices
!!       call print('determine_multiple_edges: next for '//trim(filename),next(1:max_next+1))

       max_group_size=0
       do kk=1,max_next
          max_group_size=max(max_group_size,next(kk+1) - next(kk)) 
       enddo
     !!  write(*,'(2(a,i0))') 'determine_multiple_edges: max_next=',max_next
     !!  write(*,'(2(a,i0))') 'determine_multiple_edges: edges with up to max_group_size=',max_group_size,' repetitions'

       !!allocate(group_end(0:max_group_size))
       allocate(group_start(max_group_size+1))
       allocate(group_index(max_temp))
       mm=0
       !!group_end(0)=0
       do group_size=1,max_group_size
     !!  write(*,'(2(a,i0))') 'test group_size=',group_size
          group_start(group_size)=mm+1
          do kk=1,max_next
             if(next(kk+1) - next(kk) == group_size) then
                mm=mm+1
                group_index(mm)=sort_index(next(kk)) ! points to the first edge in each group
             endif
          enddo
        !!  group_end(group_size)=mm
     !!  write(*,'(2(a,i0))') 'group_end(',group_size,')=',group_end(group_size)
       enddo
          group_start(max_group_size+1)=mm+1


!!    !!  call print('determine_multiple_edges: group_end ',group_end)                                
!!      call print('determine_multiple_edges: group_start ',group_start)                                
!!      do group_size=1,max_group_size
!!       write(*,'(2(a,i0))') 'group for group_size=',group_size
!!     !! call print('determine_multiple_edges: group_index ',group_index(group_end(group_size-1)+1:group_end(group_size)))                                
!!      call print('determine_multiple_edges: group_index ',group_index(group_start(group_size):group_start(group_size+1)-1))                                
!!      enddo

   !!   write(*,'(a,i0)') 'determine_multiple_edges: size(group_start)=',size(group_start)


 1111 continue

end subroutine determine_multiple_edges
!***************************************************************************
subroutine plot_arrows(string,nodes,v_edges,selection,filename,file_extension,arrow_direction,append,gnu_script_unit)
type(valued_egde_type)             :: v_edges
character(len=*)                   :: string
real(kind=rk)                      :: nodes(:,:)
character(len=*)                   :: filename,file_extension
logical,optional                   :: append
logical                            :: loc_append
integer                            :: gnu_script_unit
!!integer,dimension(:),optional      :: region
integer,dimension(:),optional      :: selection
integer,optional                   :: arrow_direction
integer                            :: arrows_unit 
integer                            :: i1,i2
real(kind=rk)                      :: mue,vue       
integer                            :: n1,n2
real(kind=rk)                      :: x1,y1,z1,x2,y2,z2
integer                            :: ll,mm
integer                            :: min_lt,max_lt
integer                            :: nmax
integer                            :: actual_line_type
integer                            :: max_edges
character(len=16)                   :: arrow_dir
!!character(len=8)                   :: arrow_dir='backhead'
!!character(len=8)                   :: arrow_dir='head'
!!character(len=8)                   :: arrow_dir='heads'
!!character(len=8)                   :: arrow_dir='nohead'
!real(kind=rk)                      :: arrow_head_size=0.15

            loc_append=.false.
            if(present(append)) loc_append=append

            write(*,*) 'plot_arrows: will open file ',trim(filename),' at ',trim(string)
            if(loc_append) then
               open(newunit=arrows_unit,file=trim(filename)//'.'//trim(file_extension),access='append')
            else
               open(newunit=arrows_unit,file=trim(filename)//'.'//trim(file_extension))
               write(gnu_script_unit,'(5a)') &
                              & 'if( plot_',trim(file_extension),'==1) load "',trim(filename)//'.'//trim(file_extension),'"'
               call delete_file(trim(filename)//'.'//trim(file_extension))
            endif



               ! Folgendes falsch, wenn selection present
                 min_lt=1
                 min_lt=minval(v_edges%type)
                 max_lt=maxval(v_edges%type)

            if(present(arrow_direction)) then
                    select case(arrow_direction) 
                    case(-1)
                      arrow_dir='backhead as 1'
                    case(0)
                      arrow_dir='nohead'
                    case(+1)
                      arrow_dir='head as 1'
                    case default
                          write(*,*) 'plot_arrows: error; arrow_direction must have the values -1,0,+1'
                          stop       'plot_arrows: error; arrow_direction must have the values -1,0,+1'
                    end select
            endif

            nmax=size(nodes,2)


          if(present(selection)) then
             max_edges=size(selection)
          !!  write(*,'(a,i0,3a)') '# ',max_edges,' edges for "',trim(filename)//'.'//trim(file_extension),'"'
            !!    call print('plot_arrows: v_edges%index(selection,:) ',v_edges%index(selection,:))
             if(maxval(selection) > v_edges%number ) then
                write(*,'(a,i0)') 'plot_arrows: selection has an element larger than v_edges%number=',v_edges%number
                call print('plot_arrows: edges for '//trim(filename),v_edges%index)
                call print('plot_arrows: selection ',selection)
                stop 'plot_arrows: selection has an element larger than v_edges%number'
             elseif(minval(selection) <= 0  ) then
                write(*,'(a,i0)') 'plot_arrows: selection has an element < 1 '
                call print('plot_arrows: edges for '//trim(filename),v_edges%index)
                call print('plot_arrows: selection ',selection)
                stop 'plot_arrows: selection has an element < 1 '
             endif
          else
             max_edges=v_edges%number


          if(maxval(v_edges%index) > size(nodes,1) .or. minval(v_edges%index) < 1 ) then
                  write(*,'(a,i0)') '---------------------------------------------------------------------------------------- '
                  write(*,'(a,i0)') 'plot_arrows: error, v_edges%index does not fit into the addressing area of nodes '
                  write(*,'(3a)') 'for "',trim(filename),'"'
                  write(*,'(3a)') 'at "',trim(string),'"'
                call print('plot_arrows: edges for '//trim(filename),v_edges%index)
                call print('plot_arrows: nodes ',nodes)
                stop 'plot_arrows: v_edges%index does not fit into the addressing area of nodes '
          endif
          endif


            write(arrows_unit,'(a,i0,3a)') '# ',max_edges,' edges for "',trim(string),'"'



                 mue=0.97
                 vue=0.03

                 mue=0.97
                 vue=0.10  ! head

                 mue=0.99
                 vue=0.01  

       do mm=1,max_edges 

          if(present(selection)) then
            ll=selection(mm)
          else
            ll=mm
          endif

 !!        if(present(region)) then
 !!           if(region(v_edges%index(ll,1)) == 0  ) then
 !!                   i1=2; i2=1
 !!           elseif(region(v_edges%index(ll,2)) == 0  ) then
 !!                   i1=1; i2=2
 !!           elseif(region(v_edges%index(ll,1)) <  region(v_edges%index(ll,2)) ) then
 !!                   i1=1; i2=2
 !!           else
 !!                   i1=2; i2=1
 !!           endif
 !!        else
                    i1=1; i2=2
 !!        endif

!!         Syntax of gnuplot set arrow:
!!         set arrow {<tag>} {from <position>} {to|rto <position>}
!!         { {arrowstyle | as <arrow_style>}
!!            | { {nohead | head | backhead | heads}
!!                {size <length>,<angle>{,<backangle>}}
!!                {filled | empty | nofilled}
!!                {front | back}
!!                { {linestyle | ls <line_style>}
!!                  | {linetype | lt <line_type>}
!!                    {linewidth | lw <line_width} 
!!                } 
!!              } 
!!         }

            write(arrows_unit,'(a,i0,a,i0)') '# edge ',v_edges%index(ll,i1),'->',v_edges%index(ll,i2)

            write(arrows_unit,'(a    )',advance='no') 'set arrow from '                 

            n1=v_edges%index(ll,i1); n2=v_edges%index(ll,i2)
            if( nmax >= 3 )  then
               x1=nodes(n1,1); x2=nodes(n2,1)
               y1=nodes(n1,2); y2=nodes(n2,2)
               z1=nodes(n1,3); z2=nodes(n2,3)

               write(arrows_unit,'(e12.5,a,e12.5,a,e12.5)',advance='no') x1+vue*(x2-x1),', ',y1+vue*(y2-y1),', ',z1+vue*(z2-z1)
               write(arrows_unit,'(a    )',advance='no') ' to '                 
               write(arrows_unit,'(e12.5,a,e12.5,a,e12.5)',advance='no') x1+mue*(x2-x1),', ',y1+mue*(y2-y1),', ',z1+mue*(z2-z1)

            else
               x1=nodes(n1,1); x2=nodes(n2,1)
               y1=nodes(n1,2); y2=nodes(n2,2)

               write(arrows_unit,'(e12.5,a,e12.5,a,e12.5)',advance='no') x1+vue*(x2-x1),', ',y1+vue*(y2-y1)
               write(arrows_unit,'(a    )',advance='no') ' to '                 
               write(arrows_unit,'(e12.5,a,e12.5,a,e12.5)',advance='no') x1+mue*(x2-x1),', ',y1+mue*(y2-y1)
            endif

                 actual_line_type=v_edges%type(ll)

                 write(arrows_unit,'(2a)',advance='no')       ' ',     trim(arrow_dir)
               !!write(arrows_unit,'(a,e10.2)',advance='no')  ' size ',arrow_head_size
                 write(arrows_unit,'(a,i0)',advance='no')    ' lt ',  actual_line_type
                 write(arrows_unit,'(a,e10.2)',advance='yes') ' lw ',  linewidth_scale*v_edges%value(ll)
       enddo

!!  if(present(line_type)) then
!!
!!      do ii=min_lt,max_lt
!!      !!write(*,*) ii,ubound(all_lt),lbound(all_lt)
!!          write(arrows_unit,'(a,i0,2(a,e10.2),a,i0)') ' set label "',all_lt(ii),'"   at ',xlabel,',',ylabel-(ii-min_lt)*diff_label,' left tc lt ',all_lt(ii)
!!       enddo
!!  else

!?    if(has_selection) then
!?    write(arrows_unit,'(2(a,e10.2),a,i0)') ' set label "selection"         at ',xlabel,',',ylabel-2*diff_label,' left tc lt ',lt_selection
!?    endif

!!  endif

            close(arrows_unit)

end subroutine plot_arrows
!***************************************************************************
subroutine append_plot(filename,suppress)
character(len=*)  :: filename
logical,optional          :: suppress
logical                   :: loc_suppress
integer           :: unit
logical,save     :: first_time=.true.

   loc_suppress=.false.
   if(present(suppress)) loc_suppress=suppress

   if(first_time) then
   open(newunit=unit,file='all_plots.gnu_script')
      write(unit,'(a,a,a)') 'load "parameter.gnu_script"'
   call delete_file('all_plots.gnu_script')
   call delete_file(filename)
   else
   open(newunit=unit,file='all_plots.gnu_script',position='append')
   call delete_file(filename)
   endif
   first_time=.false.
   if(.not. loc_suppress) then
      write(unit,'(a,a,a)') 'load "',trim(filename),'"'
   else
           write(unit,'(a,a,a)') 'if(suppress_large_pictures==0) load "',trim(filename),'"'
   endif
 !!  write(0,'(a,a,a)') 'append plot file "',trim(filename),'"'
   close(unit)

end subroutine append_plot
!***************************************************************************
subroutine delete_file(filename)
! to be called directly after opening a temporary file
character(len=*)  :: filename
integer           :: unit
logical,save     :: first_time=.true.

   if(first_time) then
      open(newunit=unit,file='all_delete.sh')
      write(unit,'(a,a,a)') 'echo off '
      write(unit,'(a,a,a)') 'rm mm '
   else
      open(newunit=unit,file='all_delete.sh',position='append')
   endif
   first_time=.false.
   write(unit,'(a,a,a)') 'rm ',trim(filename)
   close(unit)

end subroutine delete_file
!***************************************************************************
SUBROUTINE write_plot_script_begin(filename,title,gnu_script_unit)
character(len=*)          :: filename
character(len=*),optional :: title
integer                   :: gnu_script_unit

     open(newunit=gnu_script_unit,file=trim(filename)//'.gnu_script')
     write(gnu_script_unit,'(4a)') 'load "','parameter.gnu_script','"'
     write(gnu_script_unit,'(4a)') 'unset label'
     write(gnu_script_unit,'(4a)') 'unset xlabel'
     write(gnu_script_unit,'(4a)') 'unset ylabel'
     write(gnu_script_unit,'(4a)') 'unset arrow'
     if(present(title)) then
             write(gnu_script_unit,'(4a)') 'if(plot_title ==1) set title "',trim(title),'"'
     else
             write(gnu_script_unit,'(4a)') 'if(plot_title ==1) set title test_case." ',trim(filename),'"'
     endif

end subroutine write_plot_script_begin
!***************************************************************************
subroutine write_plot_script_end(filename,dim,minmax,suppress,gnu_script_unit)
character(len=*)          :: filename
integer                   :: gnu_script_unit
logical,optional          :: suppress
integer,optional          :: dim
real(kind=rk),dimension(:,:),optional :: minmax
integer                   :: loc_dim

   if(present(dim)) then
        loc_dim=dim
   else
        loc_dim=2
   endif

     if(loc_dim==1 .or.loc_dim==2) then
        if(present(minmax)) then
           write(gnu_script_unit,'(3(a,e12.2))') 'set xrange[',minmax(1,1),':',minmax(1,2),']'
           write(gnu_script_unit,'(3(a,e12.2))') 'set yrange[',minmax(2,1),':',minmax(2,2),']'
        else
           write(gnu_script_unit,'(4a)') 'set xrange[xmin:xmax]'
           write(gnu_script_unit,'(4a)') 'set yrange[ymin:ymax]'
        endif
           write(gnu_script_unit,'(4a)') 'plot "',trim(filename),'.gnu_dat" title ""'
     else
        if(present(minmax)) then
           write(gnu_script_unit,'(3(a,e12.2))') 'set xrange[',minmax(1,1),':',minmax(1,2),']'
           write(gnu_script_unit,'(3(a,e12.2))') 'set yrange[',minmax(2,1),':',minmax(2,2),']'
           if(size(minmax,1) >=3) then
              write(gnu_script_unit,'(3(a,e12.2))') 'set zrange[',minmax(3,1),':',minmax(3,2),']'
           endif
        else
           write(gnu_script_unit,'(4a)') 'set xrange[xmin:xmax]'
           write(gnu_script_unit,'(4a)') 'set yrange[ymin:ymax]'
           write(gnu_script_unit,'(4a)') 'set zrange[zmin:zmax]'
        endif
        write(gnu_script_unit,'(4a)') 'splot "',trim(filename),'.gnu_dat" title ""'
     endif

     write(gnu_script_unit,'(4a)') 'pause -1 "',trim(filename),'"'

     call write_gnuplot_trailer(gnu_script_unit,filename)

     close(gnu_script_unit)

     call append_plot(trim(filename)//'.gnu_script',suppress)
end subroutine write_plot_script_end
!***************************************************************************
subroutine write_gnuplot_trailer(gnu_script_unit,filename)
integer                    :: gnu_script_unit
character(len=*)           :: filename

     write(gnu_script_unit,'(*(a))') 'if(plot_type==1) set output "',trim(filename),'.svg"'
     write(gnu_script_unit,'(*(a))') 'if(plot_type==1) set terminal svg size 750,600 dynamic noenhanced '

     write(gnu_script_unit,'(*(a))') 'if(plot_type==2) set output "',trim(filename),'.pdf"'
     write(gnu_script_unit,'(*(a))') 'if(plot_type==2) set terminal pdf '

     write(gnu_script_unit,'(*(a))') 'if(plot_type==3) set output "',trim(filename),'.png"'
     write(gnu_script_unit,'(*(a))') 'if(plot_type==3) set terminal png crop '

     write(gnu_script_unit,'(*(a))') 'if(plot_type!=0) replot'
     write(gnu_script_unit,'(*(a))') 'if(plot_type!=0) set terminal pop'

end subroutine write_gnuplot_trailer
!***************************************************************************
end module
