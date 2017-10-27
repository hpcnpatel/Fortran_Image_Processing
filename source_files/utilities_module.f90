module utilities_module
use constants_module
implicit none

character(len=*),parameter :: digits='0123456789'
character(len=*),parameter :: alpha= 'abcdefghijklmnopqrstuvwxyz' &
                                   //'ABCDEFGHIJKLMNOPQRSTUVWXYZ' &
                                   //'_'
character(len=*),parameter :: alpha_num=digits//alpha

integer               :: line_length=140

interface extract_command_parameter
   module procedure extract_command_parameter_string
   module procedure extract_command_parameter_integer
   module procedure extract_command_parameter_real
end interface

interface print
  module procedure print_array_1D
  module procedure print_array_2D
  module procedure print_array_1D_real
  module procedure print_dense_matrix
end interface print

contains
!*****************************************************************************************************
function string_of(integer)
! delivers an integer as string of correct length
integer  :: integer
!!character(len=:),allocatable   :: string_of
character(len=20)               :: string_of
character(len=64)              :: temp

write(temp,*) integer

string_of=trim(adjustl(temp))

end function string_of
!*****************************************************************************************************
subroutine test_all_keywords(command_line,syntax)
character(len=*)                  :: command_line
character(len=*)                  :: syntax
character(len=32),dimension(20)   :: keywords
character(len=32),dimension(20)   :: cmd_keywords
integer                           :: nn,nmax
integer                           :: cmd_nn,cmd_nmax

 call find_command_keywords(syntax,keywords,nmax) 
 call find_command_keywords(command_line,cmd_keywords,cmd_nmax) 

 command_loop: do cmd_nn=1,cmd_nmax
 do nn=1,nmax
   if(cmd_keywords(cmd_nn)==keywords(nn)) cycle command_loop
 enddo
       write(*,'(a,a,a)') 'test_all_keywords: "',trim(cmd_keywords(cmd_nn)),'" is not an allowed keyword in the command line '
       write(*,'(3a)') '"',trim(command_line),'"'
       write(*,'(a)') 'The syntax of the command line is'
       write(*,'(3a)') '"',trim(syntax),'"'
       write(*,'(a)') 'The allowed keywords are'
       do nn=1,nmax
       write(*,*) trim(keywords(nn))
       enddo
       stop 'cmd_keywords(nn): wrong switch in command line'

 enddo command_loop

end subroutine test_all_keywords
!*****************************************************************************************************
subroutine find_command_keywords(syntax,keywords,nmax)
character(len=*)               :: syntax
character(len=*),dimension(:)  :: keywords
integer                        :: nn,nmax
integer                        :: ind,ind_old,ind2,ind3
integer                        :: len_syntax

 len_syntax=len(syntax)


 nmax=0
 ind_old=1
 do  nn=1,size(keywords)                         ! loop must be left by 'exit'
    ind=index(syntax(ind_old:),' -')+ind_old-1   ! we search for ' -' in front of a keyword
    if(ind==ind_old-1) then                      ! in this case there is no other keyword
            !!write(*,*) 'exit 2'
            exit
    else

      ind2=scan(syntax(ind+2:ind+2),alpha)+ind+1

      if(ind2>ind+1) then                       ! the first character of the string is a letter
         ind3=verify(syntax(ind+2:),alpha_num)+ind+1
         if(ind3==ind+1) then                    ! the remaining string is a keyword
                 keywords(nn)=syntax(ind+1:len_syntax)
                 nmax=nn
                 !!write(*,*) 'exit 1'
                 exit                            ! this will be the last keyword
         else                                    ! the string up to ind3-1 is a keyword
                 keywords(nn)=syntax(ind+1:ind3-1)
                 nmax=nn
                 ind_old=ind3
         endif

      endif
    endif
 enddo
    if(nmax== size(keywords)) then
            write(*,*) 'find_command_keywords: size(keywords) "',size(keywords),'" is too small'
            stop       'find_command_keywords: size keywords is too small'
    endif
   !!    do nn=1,nmax
   !!    write(*,*) nn,trim(keywords(nn))
   !!    enddo

end subroutine find_command_keywords
!*****************************************************************************************************
function extract_command_parameter_string(command_line,switch,stop_on_error,value,syntax) result(result_value)
integer           :: result_value
logical,optional  :: stop_on_error
logical           :: loc_stop_on_error
character(len=*)  :: command_line,syntax
character(len=*)  :: switch
character(len=*)  :: value
integer           :: ind
integer           :: ind_end
integer           :: lensw
integer,save      :: nmax=0
integer           :: nn,nn_switch,nn_next
character(len=32),dimension(20),save   :: keywords

  loc_stop_on_error=.false.
if(present(stop_on_error)) then
  loc_stop_on_error=stop_on_error
endif


if(nmax==0) call find_command_keywords(syntax,keywords,nmax) ! called only the first time

 call scan_for_keyword(switch,keywords,nmax,nn_switch,1,ind_end)
       !!  write(*,'(a,i0,a,a,a)') ' found keyword(',nn_switch,') "',trim(keywords(nn_switch)),'"'
if(nn_switch==0) then
       write(*,'(a,a,a)') 'extract_command_parameter: "',trim(switch),'" is not an allowed keyword '
       write(*,'(a)') trim(command_line)
       write(*,'(a)') 'syntax of command line is'
       write(*,'(a)') trim(syntax)
       write(*,'(a)') ' allowed keywords are'
       do nn=1,nmax
       write(*,*) trim(keywords(nn))
       enddo
       stop 'extract_command_parameter: wrong switch in command line'
endif



 ind=index(command_line,trim(keywords(nn_switch)))
 lensw=len_trim(keywords(nn_switch))

 ! we scan now the begin of all potential keywords to find the end of the string after the switch
 call scan_for_keyword(command_line,keywords,nmax,nn_next,ind+lensw,ind_end)

 if(ind==0) then
    if(loc_stop_on_error) then
       print*,'extract_command_parameter: did not found "',trim(switch),'" in command line'
       print*,trim(command_line)
       write(*,'(a)') 'syntax of command line is'
       write(*,'(a)') trim(syntax)
       stop 'extract_command_parameter: parameter not found in command line'
    endif
    result_value=1
 else
    result_value=0
    value=adjustl(command_line(ind+lensw:ind_end-1))
    write(*,'(a,a,a,a,a)') 'extract_command_parameter: found "',trim(value),'" for "',trim(switch),'"'
 endif


end function extract_command_parameter_string
!*****************************************************************************************************
subroutine scan_for_keyword(string,keywords,nmax,nn_switch,start,end)
character(len=*)              :: string
character(len=*),dimension(:) :: keywords
integer                       :: start
integer                       :: temp_end
integer                       :: end
integer                       :: nmax
integer                       :: nn,nn_switch

 nn_switch=0
 end=huge(end)
 !! We look for the first keyword to appear
 do nn=1,nmax
          !!  write(*,*) 'scan_for_keyword: testing for ',trim(keywords(nn))
    temp_end=index(string(start:),trim(keywords(nn)))+start-1
    if(temp_end>start-1) then
            if(temp_end < end) then
               end=temp_end
             !!  write(*,*) 'scan_for_keyword: found ',trim(keywords(nn))
               nn_switch=nn
            endif
    endif
 enddo
 if(nn_switch==0 ) then
         end=len(string)
 !!else
       !!  write(*,'(a,i0,a,a,a)') ' found next keyword(',nn_switch,') "',trim(keywords(nn_switch)),'"'
 endif
 end subroutine scan_for_keyword
!*****************************************************************************************************
function extract_command_parameter_real(command_line,switch,stop_on_error,value,syntax) result(result_value)
integer           :: result_value
logical,optional  :: stop_on_error
logical           :: loc_stop_on_error
character(len=*)  :: command_line,syntax
character(len=*)  :: switch
real(kind=rk)     :: value
integer           :: ind,ind_end
integer           :: lensw
character(len=32),dimension(20),save   :: keywords
integer,save      :: nmax=0
integer           :: nn
integer           :: nn_switch,nn_next
integer           :: status

  loc_stop_on_error=.false.
if(present(stop_on_error)) then
  loc_stop_on_error=stop_on_error
endif


if(nmax==0) call find_command_keywords(syntax,keywords,nmax) ! called only the first time

 call scan_for_keyword(switch,keywords,nmax,nn_switch,1,ind_end)

  ! first test if switch is allowed
if(nn_switch==0) then
       write(*,'(a,a,a)') 'extract_command_parameter_real: "',trim(switch),'" is not an allowed keyword '
       write(*,'(a)') trim(command_line)
       write(*,'(a)') 'syntax of command line is'
       write(*,'(a)') trim(syntax)
       write(*,'(a)') ' allowed keywords are'
       do nn=1,nmax
       write(*,*) trim(keywords(nn))
       enddo
       stop 'extract_command_parameter_real: wrong switch in command line'
endif

 ind=index(command_line,trim(keywords(nn_switch)))
 lensw=len_trim(keywords(nn_switch))

 ! we scan now the begin of all potential keywords to find the end of the string after the switch
 call scan_for_keyword(command_line,keywords,nmax,nn_next,ind+lensw,ind_end)

 if(ind==0) then
    if(loc_stop_on_error) then
       write(*,'(a,a,a)') 'extract_command_parameter_real: did not found "',trim(switch),'" in command line'
       write(*,'(a)') trim(command_line)
       write(*,'(a)') 'syntax of command line is'
       write(*,'(a)') trim(syntax)
       stop 'extract_command_parameter_real: parameter not found in command line'
    endif
    result_value=1
 else
    result_value=0
    read(command_line(ind+lensw:ind_end),*,iostat=status) value
    if(status/=0) then
       if(loc_stop_on_error) then
            write(*,'(3a)') 'extract_command_parameter_real: error reading integer parameter for switch "',trim(switch),'"'
            write(*,'(3a)') 'found no integer in "',command_line(ind+lensw:ind_end),'"'
            stop 'extract_command_parameter_real: error reading integer parameter for switch'
       endif
    endif
    write(*,'(a,g0.3,a,a,a)') 'extract_command_parameter_real: found "',value,'" for switch "',trim(switch),'"'
 endif


end function extract_command_parameter_real
!*****************************************************************************************************
function extract_command_parameter_integer(command_line,switch,stop_on_error,value,syntax) result(result_value)
integer           :: result_value
logical,optional  :: stop_on_error
logical           :: loc_stop_on_error
character(len=*)  :: command_line,syntax
character(len=*)  :: switch
integer           :: value
integer           :: ind,ind_end
integer           :: lensw
character(len=32),dimension(20),save   :: keywords
integer,save      :: nmax=0
integer           :: nn
integer           :: nn_switch,nn_next
integer           :: status

  loc_stop_on_error=.false.
if(present(stop_on_error)) then
  loc_stop_on_error=stop_on_error
endif


if(nmax==0) call find_command_keywords(syntax,keywords,nmax) ! called only the first time

 call scan_for_keyword(switch,keywords,nmax,nn_switch,1,ind_end)

  ! first test if switch is allowed
if(nn_switch==0) then
       write(*,'(a,a,a)') 'extract_command_parameter: "',trim(switch),'" is not an allowed keyword '
       write(*,'(a)') trim(command_line)
       write(*,'(a)') 'syntax of command line is'
       write(*,'(a)') trim(syntax)
       write(*,'(a)') ' allowed keywords are'
       do nn=1,nmax
       write(*,*) trim(keywords(nn))
       enddo
       stop 'extract_command_parameter: wrong switch in command line'
endif

 ind=index(command_line,trim(keywords(nn_switch)))
 lensw=len_trim(keywords(nn_switch))

 ! we scan now the begin of all potential keywords to find the end of the string after the switch
 call scan_for_keyword(command_line,keywords,nmax,nn_next,ind+lensw,ind_end)

 if(ind==0) then
    if(loc_stop_on_error) then
       write(*,'(a,a,a)') 'extract_command_parameter: did not found "',trim(switch),'" in command line'
       write(*,'(a)') trim(command_line)
       write(*,'(a)') 'syntax of command line is'
       write(*,'(a)') trim(syntax)
       stop 'extract_command_parameter: parameter not found in command line'
    endif
    result_value=1
 else
    result_value=0
         !!   write(*,'(3a)') 'try to read integer from  "',command_line(ind+lensw:ind_end),'"'
    read(command_line(ind+lensw:ind_end),*,iostat=status) value
    if(status/=0) then
       if(loc_stop_on_error) then
            write(*,'(3a)') 'extract_command_parameter_integer: error reading integer parameter for switch "',trim(switch),'"'
            write(*,'(3a)') 'found no integer in "',command_line(ind+lensw:ind_end),'"'
            stop 'extract_command_parameter_integer: error reading integer parameter for switch'
       endif
    endif
    write(*,'(a,i0,a,a,a)') 'extract_command_parameter: found integer "',value,'" for switch "',trim(switch),'"'
 endif


end function extract_command_parameter_integer
!*****************************************************************************************************
function extract_from(string,open,close,start_end,value) result(result_value)
character(len=*)             :: string
character(len=*),optional    :: open
character(len=*)             :: close
!!character(len=:),allocatable :: value
character(len=*)           :: value
integer                      :: result_value
integer                      :: start_end
integer                      :: len_open
integer                      :: ind,ind_end

result_value=0
len_open=0
ind=start_end

   if(present(open)) then
      len_open=len(open)
      ind=index(string(start_end:),open) + start_end-1
      if(ind==start_end-1) then
              goto 1000
      endif
   endif

   ind_end=index(string(ind+len_open:),close) + ind+len_open-1
   if(ind_end/=ind) then
           value=string(ind+len_open:ind_end-1)
           print*,'value=',trim(value)
           result_value=ind_end-ind-len_open
   else
           goto 1000
   endif

start_end=ind_end+len(close)

1000 continue

end function extract_from
!*********************************************************************************************************
subroutine print_array_1D(string,array,lower_bound)
character(len=*)      :: string
integer,dimension(:)  :: array
integer,optional      :: lower_bound
integer               :: loc_lower_bound
integer               :: nblock
integer               :: mm,nn
integer               :: size_array

character(len=7)      :: fmt_int
character(len=7)      :: fmt_chr
integer               :: field_length
integer               :: max_array,max_field_number

size_array=size(array)
write(*,'(2a)') trim(string)
if(size_array==0) then
   write(*,'(a)') ' has no elements'
   goto 2222
endif
write(*,'(a,i0,a)') '',size_array,' elements'

loc_lower_bound=1
if(present(lower_bound) ) loc_lower_bound=lower_bound

max_array=maxval(array)
max_field_number=max(loc_lower_bound-1+size_array,max_array)

!!write(*,'(a,i0)') 'max_array=',max_array
!!write(*,'(a,i0)') 'size_array=',size_array

field_length=int(log10(real(max_field_number)),4) + 1 + 1
if(minval(array) < 0 ) field_length=field_length+1

!!write(*,'(a,i0)') 'field_length=',field_length

nblock=line_length/field_length

fmt_int='(i'//trim(string_of(field_length))//')'
fmt_chr='(a'//trim(string_of(field_length))//')'

!!print*,' fmt_int=',fmt_int
!!print*,' fmt_chr=',fmt_chr

do nn=1,size_array,nblock
     write(*,*)
     write(*,fmt_chr,advance='no') 'no:'
  do mm=nn,min(size_array,nn+nblock-1)
     write(*,fmt_int,advance='no') mm+loc_lower_bound-1
  enddo
     write(*,'(a)',advance='yes') ''

     write(*,fmt_chr,advance='no') '   '
  do mm=nn,min(size_array,nn+nblock-1)
     write(*,fmt_chr,advance='no') '---------------------------'
  enddo
     write(*,'(a)',advance='yes') ''

     write(*,fmt_chr,advance='no') '   '
  do mm=nn,min(size_array,nn+nblock-1)
     write(*,fmt_int,advance='no') array(mm)
  enddo
     write(*,'(a)',advance='yes') ''
enddo

2222 continue
     write(*,'(2a)') 'end: ',trim(string)
     write(*,'(a)',advance='yes') ''

end subroutine print_array_1D
!*********************************************************************************************************
subroutine print_array_1D_real(string,array,threshould)
character(len=*)      :: string
real(kind=rk),dimension(:)  :: array
integer               :: nblock
integer               :: mm,nn
integer               :: size_array
character(len=7)      :: fmt_int
character(len=10)      :: fmt_fp
character(len=7)      :: fmt_chr
integer               :: field_length
real(kind=rk),optional          :: threshould
real(kind=rk)                   :: loc_threshould


!loc_threshould=print_threshould
if(present(threshould)) loc_threshould=threshould

size_array=size(array)
write(*,'(2a)') trim(string)
if(size_array==0) then
   write(*,'(a)') ' has no elements'
   goto 2222
endif
write(*,'(a,i0,a)') '',size_array,' elements'
write(*,'(a,a,2(e10.3,a))') trim(string),': values in [',minval(array),',',maxval(array),']'



field_length=10
if(minval(array) < 0 ) field_length=field_length+1

nblock=line_length/field_length

fmt_int='(i'//trim(string_of(field_length))//')'
fmt_fp='(e'//trim(string_of(field_length))//'.3)'
fmt_chr='(a'//trim(string_of(field_length))//')'

!!print*,' fmt_int=',fmt_int
!!print*,' fmt_fp=',fmt_fp
!!print*,' fmt_chr=',fmt_chr

do nn=1,size_array,nblock
     write(*,*)
     write(*,'(a4)',advance='no') 'no:'
  do mm=nn,min(size_array,nn+nblock-1)
     write(*,fmt_int,advance='no') mm
  enddo
     write(*,'(a)',advance='yes') ''
     write(*,'(a4)',advance='no') '   '
  do mm=nn,min(size_array,nn+nblock-1)
     if(abs(array(mm)) <= loc_threshould) then
        write(*,fmt_chr,advance='no') '-'
     else
        write(*,fmt_fp,advance='no') array(mm)
     endif
  enddo
     write(*,'(a)',advance='yes') ''
enddo

2222 continue
     write(*,'(2a)') 'end: ',trim(string)
     write(*,'(a)',advance='yes') ''

end subroutine print_array_1D_real
!*********************************************************************************************************
subroutine print_array_2D(string,array,mask)
character(len=*)      :: string
integer,dimension(:,:):: array
logical,dimension(:),optional  :: mask
integer               :: nblock
integer               :: mm,nn,ll
integer               :: size_array_1,size_array_2

character(len=7)      :: fmt_int
character(len=7)      :: fmt_chr
integer               :: field_length
integer               :: max_array,max_field_number

size_array_1=size(array,1)
size_array_2=size(array,2)
write(*,'(2a)') trim(string)
write(*,'(2(a,i0),a)') 'array with ',size_array_1,'x',size_array_2,' elements'

if(size_array_1<=0 .or. size_array_2 <=0) goto 1000

max_array=maxval(array)
max_field_number=max(size_array_2,max_array)

field_length=int(log10(real(max_field_number))) + 1 + 1
if(minval(array) < 0 ) field_length=field_length+1

nblock=line_length/field_length

fmt_int='(i'//trim(string_of(field_length))//')'
fmt_chr='(a'//trim(string_of(field_length))//')'

!!print*,' field_length=',field_length
!!print*,' fmt_int=',fmt_int
!!print*,' fmt_chr=',fmt_chr
!!print*,' nblock=',nblock

if(.not. present(mask)) then
do nn=1,size_array_2,nblock
     write(*,*)
     write(*,fmt_chr,advance='no') 'no:'
  do mm=nn,min(size_array_2,nn+nblock-1)
     write(*,fmt_int,advance='no') mm
  enddo
     write(*,'(a)',advance='yes') ''

     write(*,fmt_chr,advance='no') '   '
  do mm=nn,min(size_array_2,nn+nblock-1)
     write(*,fmt_chr,advance='no') '---------------------------'
  enddo
     write(*,'(a)',advance='yes') ''

  do ll=1,size_array_1
     !!write(*,fmt_chr,advance='no') '   '
     write(*,fmt_int,advance='no') ll
  do mm=nn,min(size_array_2,nn+nblock-1)
     write(*,fmt_int,advance='no') array(ll,mm)
  enddo
     write(*,'(a)',advance='yes') ''
  enddo
     write(*,'(a)',advance='yes') ''
enddo

else
do nn=1,size_array_1,nblock
     write(*,fmt_chr,advance='no') 'no:'
  do mm=nn,min(size_array_2,nn+nblock-1)
  if(mask(mm)) write(*,fmt_int,advance='no') mm
  enddo
     write(*,'(a)',advance='yes') ''
  do ll=1,size_array_1
     !!write(*,fmt_chr,advance='no') '   '
     write(*,fmt_int,advance='no') ll
  do mm=nn,min(size_array_1,nn+nblock-1)
  if(mask(mm))    write(*,fmt_int,advance='no') array(ll,mm)
  enddo
     write(*,'(a)',advance='yes') ''
  enddo
     write(*,'(a)',advance='yes') ''
enddo
endif

1000 continue
     write(*,'(2a)') 'end: ',trim(string)
     write(*,'(a)',advance='yes') ''

end subroutine print_array_2D
!*********************************************************************************************************
subroutine print_dense_matrix(string,matrix,index,threshould)
character(len=*)                :: string
real(kind=rk),dimension(:,:)    :: matrix
real(kind=rk),optional          :: threshould
real(kind=rk)                   :: loc_threshould
integer,dimension(:),optional   :: index
integer                         :: ll,nn,mm
integer                         :: max_1,max_2
integer                         :: ll_diff=15


!loc_threshould=print_threshould
if(present(threshould)) loc_threshould=threshould

!!write(*,'(2(a,i0))') 'print_dense_matrix: size(matrix,1)=',size(matrix,1),' size(matrix,2)=',size(matrix,2)
max_1=size(matrix,1)
max_2=size(matrix,2)

!!write(*,'(10e9.2)') matrix
write(*,*)
write(*,'(a)',advance='no')  'begin:'
write(*,*) trim(string)

if(maxval(abs(matrix)) > loc_threshould) then
do ll=1,max_2,ll_diff
   if(present(index)) then
     write(*,'(5x,a)',advance='no') 
     write(*,'(5x,a)',advance='no') '|'
     do mm=ll,min(ll+ll_diff-1,max_2)
        write(*,'(i6,4x)',advance='no') index(mm)
     enddo
     write(*,'(a)',advance='yes') ''
     write(*,'(5x,a)',advance='no') 
   endif

   write(*,'(5x,a)',advance='no') '|'
   do mm=ll,min(ll+ll_diff-1,max_2)
      write(*,'(i6,4x)',advance='no') mm
   enddo
   write(*,'(a)',advance='yes') ''

   write(*,'(a)',advance='yes') repeat('_',10*min(ll_diff,max_2)+5+2)
do nn=1,max_1
   if(present(index)) write(*,'(i5)',advance='no') index(nn)
   write(*,'(i5,a)',advance='no') nn,'|'
  do mm=ll,min(ll+ll_diff-1,max_2)
 !!!     write(*,'(2(a,i0))') 'nn=',nn,'  mm=',mm 
      if(abs(matrix(nn,mm)) < loc_threshould) then
      write(*,'(a10)',advance='no')    '    -    '
      else
      write(*,'(e10.3)',advance='no') matrix(nn,mm)
      endif
   enddo
   write(*,'(a)',advance='yes') ''
enddo
   write(*,'(a)',advance='yes') repeat('_',10*min(ll_diff,max_2)+5+2)
enddo
else
   write(*,*) 'matrix has no elements exceeding the threshould of ',loc_threshould
endif

write(*,'(a)',advance='no')  'end:'
write(*,*) trim(string)

end subroutine print_dense_matrix
!*********************************************************************************************************
function versioning_of_filenames(filename) result(versioned_filename)
character(len=*)                                      :: filename
character(len=len_trim(filename)+6)                   :: versioned_filename
integer,parameter                                     :: max_filenames=1000
character(len=128),dimension(max_filenames),save      :: list_of_filenames
integer,dimension(max_filenames),save                 :: version_of_filenames
integer,save                                          :: no_list=0
integer                                               :: mm

   if(no_list==0) version_of_filenames=2
   do mm=1,no_list
      if(filename==list_of_filenames(mm)) then
          write(versioned_filename,'(a,a,i0)') trim(filename),'_v',version_of_filenames(mm)
          version_of_filenames(mm)=version_of_filenames(mm)+1
          if(version_of_filenames(mm)>9999) then
              write(*,*) 'versioning_of_filenames: error, version_of_filenames>9999'
              write(*,'(3a)') 'for file "',trim(filename),'"'
              stop 'versioning_of_filenames:  error, version_of_filenames>9999'  
          endif
          goto 777
      endif
   enddo
   no_list=no_list+1
   if(no_list>max_filenames) then
     write(*,*) 'versioning_of_filenames: list list_of_filenames is exhausted'
           stop 'versioning_of_filenames: list list_of_filenames is exhausted'
   endif
   list_of_filenames(no_list)=filename
   versioned_filename=filename
   777 continue
end function versioning_of_filenames
!*****************************************************************************************************
end module utilities_module
