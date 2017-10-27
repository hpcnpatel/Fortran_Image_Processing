!==============================================================================
!> \file vtk_writter.f90
!> Module for vtk output in Fortran
!> The file holds a module with routines for vtk output
!>
!> \author Ralf Schneider, changed by Nisarg Patel
!> \date 19.08.2014

Module vtkwritter

  use constants_module
  use variable_module

  Implicit None

!  Private give_new_unit
!  Private file_err

!  interface write_vtk

!     Module Procedure write_vtk_structured_points_int4
!     Module Procedure write_vtk_structured_points_int8
!     Module Procedure write_vtk_structured_points_real8
!     Module Procedure write_vtk_structured_points_real8_tensor
!     Module Procedure write_vtk_unstructured_grid
!     Module Procedure write_vtk_polydata_int4_1D
!     Module Procedure write_vtk_polydata_int4_3D

!  End interface write_vtk

contains
 !============================================================================
subroutine write_vtk_polydata(file_name,x_vec,y_vec,z_vec)

INTEGER :: unit
CHARACTER(c_len),intent(IN)::file_name
REAL(KIND=rk8), dimension(:),allocatable ::x_vec,y_vec,z_vec

open(newunit=unit,file=trim(file_name),action='write',status='replace') 
call write_vtk_head_ascii(unit)
call write_vtk_polydata_line(unit,x_vec,y_vec,z_vec)
 CLOSE(unit)
end subroutine
!============================================================================
  !Subroutine which writes a vtk polydata 
  subroutine write_vtk_polydata_line(unit,x_vec,y_vec,z_vec)

INTEGER :: unit
!    INTEGER(Kind=rk4)              , INTENT(IN) :: un_out
!    REAL(Kind=rk4)   , Dimension(3), INTENT(IN) :: spacing,origin
!    INTEGER(kind=rk4), Dimension(3), INTENT(IN) :: extend

!    CHARACTER(c_len)               :: tmp_line, tmp_origin,tmp_real
    CHARACTER(c_len)               :: tmp_pointdata,tmp_line
!    INTEGER(Kind=rk4)              :: pointdata
REAL(KIND=rk8), dimension(:),allocatable ::x_vec,y_vec,z_vec
INTEGER ::dim,i

dim=size(x_vec)

write(unit,'(A)') 'DATASET POLYDATA'!,achar(10)
write(unit,'(A,A,I0,A,A)') 'POINTS',' ',dim+1,' ','float'

!this loop writes coordinates from 1:n, n+1 is again the first value to complete the polyline.
Do i=1,dim
write(tmp_pointdata,'(I0,A,I0,A,I0)')nint(x_vec(i)),' ',nint(y_vec(i)),' ',nint(z_vec(i))!,' ',achar(10)
write(unit,'(A)')trim(tmp_pointdata)
if(i==dim)then
    write(tmp_pointdata,'(I0,A,I0,A,I0)')nint(x_vec(1)),' ',nint(y_vec(1)),' ',nint(z_vec(1))!,' ',achar(10)
    write(unit,'(A)')trim(tmp_pointdata)
endif

enddo

write(unit,'(A,A,A,A,I0)') 'LINES',' ','1',' ',dim+2

tmp_line=''
tmp_pointdata=''

Do i=1,dim
if(i==1)then
    write(tmp_line,'(I0,A,I0)')dim+1,' ',i-1
elseif(i==dim)then
    write(tmp_line,'(A,I0,A,I0)')' ',i-1,' ',i
else
    write(tmp_line,'(A,I0)')' ',i-1
endif

!write(unit,'(A)')trim(tmp_pointdata)
tmp_pointdata=trim(tmp_pointdata)//trim(tmp_line)
!if(i==dim)write(tmp_pointdata,'(A,A,I0,A,A)')'2',' ',i,' ','0'
enddo

write(unit,'(A)')trim(tmp_pointdata)


End subroutine write_vtk_polydata_line
  !============================================================================
 !============================================================================
subroutine WRITE_VTK(ct,data,file_name)
   
 TYPE(tCt),INTENT(INOUT)                                   ::ct
 Real(kind=rk8),dimension(:,:,:),allocatable,INTENT(INOUT) ::data
 CHARACTER(c_len),intent(IN)::file_name
 REAL(Kind=rk4)   , Dimension(3) :: spacing,origin
 INTEGER(kind=rk4), Dimension(3) :: extend

 open(999,file=trim(file_name),action='write',access='stream'&
                                 &,status='replace',convert='big_endian')

spacing=(/ct%dx,ct%dy,ct%dz/)
origin =(/ct%offset_x ,ct%offset_y ,ct%offset_z/)
extend =(/ct%pixel_x ,ct%pixel_y ,ct%slices/)

call write_vtk_head_binary(999)
call write_vtk_structured_points_head(999,extend,spacing,origin)
write(999)int(data,2),achar(10)

 Close(999)
 end subroutine write_vtk
 !============================================================================
  !Subroutine which writes a head for a vtk structured_points dataset
  subroutine write_vtk_structured_points_head(un_out,  extend, &
                                              spacing, origin   )

    INTEGER(Kind=rk4)              , INTENT(IN) :: un_out
    REAL(Kind=rk4)   , Dimension(3), INTENT(IN) :: spacing,origin
    INTEGER(kind=rk4), Dimension(3), INTENT(IN) :: extend

    CHARACTER(c_len)               :: tmp_line, tmp_origin,tmp_real
    CHARACTER(c_len)               :: tmp_pointdata
    INTEGER(Kind=rk4)              :: pointdata

    write(un_out) 'DATASET STRUCTURED_POINTS',achar(10)

!# Dimension
!tmp_real=""
write(tmp_real,'(A,I0,A1,I0,A1,I0,A)')'DIMENSIONS ',extend(1),'',extend(2),'',extend(3),achar(10)
write(un_out) trim(tmp_real)

!# spacing
!tmp_line=""
Write(tmp_line,'(A,F6.4,A1,F6.4,A1,F6.4,A)')'SPACING ',spacing(1),'',spacing(2),'',spacing(3),achar(10)
write(un_out)trim(tmp_line)

!# origin
!tmp_origin=""
write(tmp_origin,'(A,F10.3,A1,F10.3,A1,F10.3,A)') 'ORIGIN ',origin(1),'',origin(2),'',origin(3),achar(10)
write(un_out)trim(tmp_origin)

!# pointdata
pointdata = extend(1) * extend(2) * extend(3)
!tmp_pointdata=""
write(tmp_pointdata,'(A,I0,A)') 'POINT_DATA ',pointdata,achar(10)
write(un_out)trim(tmp_pointdata)

!# Temp addition! 13/06/2014
!# Scalars
write(un_out) 'SCALARS scalars short',achar(10)
!write(un_out) 'SCALARS scalars double',achar(10)

!# Lockup_table
write(un_out) 'LOOKUP_TABLE default',achar(10)

End subroutine write_vtk_structured_points_head
  !============================================================================
  !> Subroutine which writes a head for a vtk binary file
  subroutine write_vtk_head_binary(un_out)

    integer(kind=ik4)                         , intent(in) :: un_out

    write(un_out) '# vtk DataFile Version 3.0',achar(10)
    write(un_out) 'vtk output',achar(10)
    write(un_out) 'BINARY',achar(10)

  End subroutine write_vtk_head_binary

!============================================================================
  !> Subroutine which writes a head for a vtk ASCII file
  subroutine write_vtk_head_ascii(nunit)

    integer                         , intent(in) :: nunit

    write(nunit,'(A)') '# vtk DataFile Version 3.0'!,achar(10)
    write(nunit,'(A)') 'vtk output'!,achar(10)
    write(nunit,'(A)') 'ASCII'!,achar(10)


  End subroutine write_vtk_head_ascii

 !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
!subroutine write_vtk_unstructured_grid(nodes, elems, no_nodes, no_elems,filename)
!
!Real(Kind=ik8)    , Dimension(:,:), intent(in) :: nodes
!Integer(Kind=ik8) , Dimension(:,:), intent(in) :: elems

!Integer(Kind=ik4)                 , intent(in) :: no_nodes, no_elems
!Character(Len=*)                 , Intent(in)  :: filename
!character(c_len)                               :: tmp_line
!integer(kind=ik4)                              :: un_out
!integer(kind=ik4)                              :: ii
!
!    un_out = give_new_unit()
!
!    call Open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!         action='write', status='replace')

!    write(un_out) '# vtk DataFile Version 3.0',achar(10)
!    write(un_out) 'vtk output',achar(10)
!    write(un_out) 'BINARY',achar(10)
!    write(un_out) 'DATASET UNSTRUCTURED_GRID',achar(10)

!    tmp_line=""
!    write(tmp_line,"(A,1X,I0,1X,A,A)")'POINTS',no_nodes,'double',achar(10)
!    write(un_out)trim(tmp_line)

!    Do ii = 1, no_nodes
!       Write(un_out)nodes(:,ii)
!    End Do

!    tmp_line=""
!    write(tmp_line,"(A,1X,I0,1X,I0,A)")'CELLS',no_elems,no_elems*4,achar(10)
!    write(un_out)trim(tmp_line)

!    Do ii = 1, no_elems
!       Write(un_out)3_4,Int(elems(:,ii)-1,4)
!    End Do

!    tmp_line=""
!    write(tmp_line,"(A,1X,I0,A)")'CELL_TYPES',no_elems,achar(10)
!    write(un_out)trim(tmp_line)

!    Do ii = 1, no_elems
!       Write(un_out)5_4
!    End Do

!    call close_big_endian_stream(un_out)

!  end subroutine write_vtk_unstructured_grid

!============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_unstructured_grid_nodelist(nodes, no_nodes, &
!       filename)
!
!    Real(Kind=ik)    , Dimension(:,:), intent(in) :: nodes
!
!    Integer(Kind=ik)                 , intent(in) :: no_nodes
!    Character(Len=*)                 , Intent(in) :: filename
!
!    character(len=mcl)                :: tmp_line
!    integer(kind=4)                   :: un_out
!    integer(kind=ik)                  :: ii
!
!    un_out = give_new_unit()
!
!    call Open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!         action='write', status='replace')
!
!    write(un_out) '# vtk DataFile Version 3.0',achar(10)
!    write(un_out) 'vtk output',achar(10)
!    write(un_out) 'BINARY',achar(10)
!    write(un_out) 'DATASET UNSTRUCTURED_GRID',achar(10)
!
!    tmp_line=""
!    write(tmp_line,"(A,1X,I0,1X,A,A)")'POINTS',no_nodes,'double',achar(10)
!    write(un_out)trim(tmp_line)
!
!    Do ii = 1, no_nodes
!       Write(un_out)nodes(:,ii)
!    End Do
!
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_unstructured_grid_nodelist
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_structured_points_int4(matrix, filename, extend, &
!                                              spacing, origin, desc, &
!                                              head, init)
!
!    integer(kind=4) , Dimension(:,:,:), intent(in) :: matrix
!    character(len=*)                  , intent(in) :: filename
!    Real(kind=ik)   , Dimension(3)    , intent(in) :: spacing,origin
!    integer(kind=ik), Dimension(3)    , intent(in) :: extend
!    character(len=*), optional        , intent(in) :: desc
!    Logical         , optional        , intent(in) :: head, init
!
!    !character(len=mcl)               :: tmp_line
!    !character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head, loc_init
!    !integer(kind=ik)                 :: pointdata
!
!    if (present(init)) then
!       loc_init = init
!    Else
!       loc_init = .FALSE.
!    End if
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!    Else
!       loc_desc = "scalars_i4"
!    End if
!
!    un_out = give_new_unit()
!    if (loc_init) then
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='replace')
!       call write_vtk_head(un_out)
!    Else
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='old', position='append')
!    End if
!
!
!    If (loc_head) then
!       !# -------
!       !# Header
!       call write_vtk_structured_points_head(un_out,  extend, &
!                                             spacing, origin   )
!    End If
!    !# -------
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
!
!    !# -------
!    !# Lockup_table
!    write(un_out) 'LOOKUP_TABLE default',achar(10)
!
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_structured_points_int4
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_structured_points_int8(matrix, filename, extend, &
!                                              spacing, origin, desc, &
!                                              head, init)
!
!    integer(kind=8) , Dimension(:,:,:), intent(in) :: matrix
!    character(len=*)                  , intent(in) :: filename
!    Real(kind=ik)   , Dimension(3)    , intent(in) :: spacing,origin
!    integer(kind=ik), Dimension(3)    , intent(in) :: extend
!    character(len=*), optional        , intent(in) :: desc
!    Logical         , optional        , intent(in) :: head, init
!
!    !character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
!    !character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head, loc_init
!!   integer(kind=ik)                 :: pointdata
!
!    if (present(init)) then
!       loc_init = init
!    Else
!       loc_init = .FALSE.
!    End if
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!    Else
!       loc_desc = "scalars_i8"
!    End if
!
!    un_out = give_new_unit()
!    if (loc_init) then
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='replace')
!       call write_vtk_head_binary(un_out)
!    Else
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='old', position='append')
!    End if
!
!    If (loc_head) then
!       !# -------
!       !# Header
!       call write_vtk_structured_points_head(un_out,  extend, &
!                                             spacing, origin   )
!    End If
!
!    !# -------
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'SCALARS '//trim(loc_desc)//' long',achar(10)
!    !# -------
!    !# Lockup_table
!    write(un_out) 'LOOKUP_TABLE default',achar(10)
!
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_structured_points_int8
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_structured_points_real8(matrix,  filename, extend, &
!                                               spacing, origin, desc, &
!                                               head, init)
!
!    Real(kind=8)    , Dimension(:,:), intent(in) :: matrix
!    character(len=*)                  , intent(in) :: filename
!    Real(kind=ik)   , Dimension(3)    , intent(in) :: spacing,origin
!    integer(kind=ik), Dimension(3)    , intent(in) :: extend
!    character(len=*), optional        , intent(in) :: desc
!    Logical         , optional        , intent(in) :: head, init
!
!    !character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
!    !character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head, loc_init
!
!!   integer(kind=ik)                 :: pointdata
!
!    if (present(init)) then
!       loc_init = init
!    Else
!       loc_init = .FALSE.
!    End if
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!
!    Else
!       loc_desc = "scalars_r8"
!
!    End if
!
!    un_out = give_new_unit()
!
!    if (loc_init) then
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='replace')
!       call write_vtk_head_binary(un_out)
!    Else
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='old', position='append')
!    End if
!    If (loc_head) then
!       !# -------
!       !# Header
!       call write_vtk_structured_points_head(un_out,  extend, &
!                                             spacing, origin   )
!    End If
!
!    !# -------
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'SCALARS '//trim(loc_desc)//' double',achar(10)
!    !# -------
!    !# Lockup_table
!    write(un_out) 'LOOKUP_TABLE default',achar(10)
!
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_structured_points_real8
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_structured_points_real8_3D1D(matrix,  filename, extend, &
!                                               spacing, origin, desc, &
!                                               head, init)
!
!    Real(kind=8)    , Dimension(:)    , intent(in) :: matrix
!    character(len=*)                  , intent(in) :: filename
!    Real(kind=ik)   , Dimension(3)    , intent(in) :: spacing,origin
!    integer(kind=ik), Dimension(3)    , intent(in) :: extend
!    character(len=*), optional        , intent(in) :: desc
!    Logical         , optional        , intent(in) :: head, init
!
!!   character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
!!   character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head, loc_init
!
!!   integer(kind=ik)                 :: pointdata
!
!    if (present(init)) then
!       loc_init = init
!    Else
!       loc_init = .FALSE.
!    End if
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!    Else
!       loc_desc = "scalars_r8"
!    End if
!
!    un_out = give_new_unit()
!    if (loc_init) then
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='replace')
!       call write_vtk_head_binary(un_out)
!    Else
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='old', position='append')
!    End if
!
!    If (loc_head) then
!       !# -------
!       !# Header
!       call write_vtk_structured_points_head(un_out,  extend, &
!                                             spacing, origin   )
!    End If
!
!    !# -------
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'SCALARS '//trim(loc_desc)//' double',achar(10)
!    !# -------
!    !# Lockup_table
!    write(un_out) 'LOOKUP_TABLE default',achar(10)
!
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_structured_points_real8_3D1D
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_structured_points_real8_tensor (&
!       matrix,  filename, extend, spacing, origin, desc, head, init)
!
!    Real(kind=8)    , Dimension(:,:,:,:,:), intent(in) :: matrix
!    character(len=*)                      , intent(in) :: filename
!    Real(kind=ik)   , Dimension(3)        , intent(in) :: spacing,origin
!    integer(kind=ik), Dimension(5)        , intent(in) :: extend
!    character(len=*), optional            , intent(in) :: desc
!    Logical         , optional            , intent(in) :: head, init
!
!!   character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
!!   character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head, loc_init
!
!!   integer(kind=ik)                 :: pointdata
!
!    if (present(init)) then
!       loc_init = init
!    Else
!       loc_init = .FALSE.
!    End if
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!    Else
!       loc_desc = "tensors_r8"
!    End if
!
!    un_out = give_new_unit()
!    if (loc_init) then
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='replace')
!       call write_vtk_head_binary(un_out)
!    Else
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='old', position='append')
!    End if
!
!    If (loc_head) then
!       !# -------
!       !# Header
!       call write_vtk_structured_points_head(un_out,  extend(3:5), &
!                                             spacing, origin    )
!    End If
!
!    !# -------
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'TENSORS '//trim(loc_desc)//' double',achar(10)
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_structured_points_real8_tensor
!
  !============================================================================
 !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_polydata_int4_1D (matrix, filename, desc, head)
!
!    Integer(kind=4) , Dimension(:)            , intent(in) :: matrix
!    character(len=*)                          , intent(in) :: filename
!    character(len=*), optional                , intent(in) :: desc
!    Logical         , optional                , intent(in) :: head
!
!    character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head!, loc_init
!
!!   integer(kind=ik)                 :: pointdata
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!    Else
!       loc_desc = "scalars_i8"
!    End if
!
!    un_out = give_new_unit()
!    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!                                   action='write', status='old', &
!                                   position='append')
!    If (loc_head) then
!       !# -------
!       !# Header
!       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix)
!       write(un_out)trim(tmp_pointdata)
!
!    End If
!
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
!    !# -------
!    !# Lockup_table
!    write(un_out) 'LOOKUP_TABLE default',achar(10)
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_polydata_int4_1D
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_polydata_int4_3D (matrix, filename, desc, head)
!
!    Integer(kind=4) , Dimension(:,:,:)        , intent(in) :: matrix
!    character(len=*)                          , intent(in) :: filename
!    character(len=*), optional                , intent(in) :: desc
!    Logical         , optional                , intent(in) :: head
!
!    character(len=mcl)               :: tmp_pointdata
!    character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_head!, loc_init
!
!!   integer(kind=ik)                 :: pointdata
!
!    if (present(head)) then
!       loc_head = head
!    Else
!       loc_head = .FALSE.
!    End if
!
!    if (present(desc)) then
!       loc_desc = desc(1:min(12,len_trim(desc)))
!    Else
!       loc_desc = "scalars_i8"
!    End if
!
!    un_out = give_new_unit()
!    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!                                   action='write', status='old', &
!                                   position='append')
!    If (loc_head) then
!       !# -------
!       !# Header
!       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix)
!       write(un_out)trim(tmp_pointdata)
!
!    End If
!
!    !# Scalars
!    write(un_out)char(10)
!    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
!    !# -------
!    !# Lockup_table
!    write(un_out) 'LOOKUP_TABLE default',achar(10)
!
!    write(un_out) matrix
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_polydata_int4_3D
!
!  !============================================================================
!  !> Subroutine which writes a vtk structured_points dataset
!  subroutine write_vtk_polydata_grid(grids,  filename, orig, init)
!
!    Real(kind=ik)   , Dimension(:,:)          , intent(in) :: grids
!    character(len=*)                          , intent(in) :: filename
!    Real(kind=ik)   , Dimension(3), optional  , intent(in) :: orig
!    Logical         , optional                , intent(in) :: init
!
!    character(len=mcl)               :: tmp_char
!!   character(len=12)                :: loc_desc
!    integer(kind=4)                  :: un_out
!    logical                          :: loc_init
!
!    integer(kind=ik)                 :: no_points,ii
!
!    if (present(init)) then
!       loc_init = init
!    Else
!       loc_init = .FALSE.
!    End if
!
!    no_points = size(grids(1,:))
!
!    un_out = give_new_unit()
!    if (loc_init) then
!
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='replace')
!       call write_vtk_head_binary(un_out)
!
!    Else
!
!       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
!            action='write', status='old', position='append')
!
!    End if
!
!    write(un_out) 'DATASET POLYDATA',achar(10)
!
!    write(tmp_char,'(A,I0,A)')"POINTS ",no_points," double"
!    write(un_out)trim(tmp_char),achar(10)
!
!    if (present(orig)) then
!       Do ii = 1, no_points
!          write(un_out)grids(:,ii)+orig
!       End Do
!    Else
!       write(un_out)grids
!    End if
!
!    write(un_out)achar(10)
!
!    call close_big_endian_stream(un_out)
!
!  end subroutine write_vtk_polydata_grid
!
  !============================================================================
  !> Function which returns new free unit
  function give_new_unit() result(new_unit)

    Integer(kind=4) :: new_unit

    Integer(kind=4) :: ii

    Logical :: unit_is_open

    Do ii = 3000, huge(new_unit)-1

       inquire(unit=ii, opened=unit_is_open)

       if( .not.unit_is_open ) then
          new_unit = ii
          Exit
       end if

    End Do

    if ( unit_is_open ) then

       WRITE(*,*)
       WRITE(*,*)'Something bad and unexpected happened during search ',&
            'for free unit'
       WRITE(*,*)'Could not find a new unit between 3000 and huge(Int(kind=4))'
       WRITE(*,*)' '
       WRITE(*,*)'PROGRAM STOPPED'
       STOP
    END IF

  End function give_new_unit

  !============================================================================
!  !> Subroutine for I/O error handling while operating on files
  SUBROUTINE file_err(in_file,io_stat)

    INTEGER             :: io_stat
    CHARACTER (LEN=*)   :: in_file

    IF (io_stat /= 0) Then
       WRITE(*,*)
       WRITE(*,"(80('='))")
       WRITE(*,"('EE ',A,T77,' EE')")   'Operation on file :'
       WRITE(*,"('EE ',A          )")   in_file
       WRITE(*,"('EE ',A,T77,' EE')")   'faild !!'
       WRITE(*,"('EE ',A,I0,T77,' EE')")'With I/O Status ',io_stat
       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
       STOP
    End IF

  END SUBROUTINE file_err

  !============================================================================
!> Subroutine for opening files with big endian encoding
!>
Subroutine open_as_big_endian_stream(unit,file,action,status,position)

    integer(kind=ik4) ,intent(in)           :: unit
    character(len=*),intent(in)           :: file,action,status
    character(len=*),intent(in), optional :: position

!   integer(kind=4)             :: ier
    character(c_len)            :: loc_pos

    If (present(position)) then
       loc_pos = position
    Else
       loc_pos='rewind'
    End If

    !**************************************************************************
    !** GFortran, Intel implementation
    Open(unit=unit, file=trim(file), action=trim(action), &
         status=trim(status), &
         access="stream", convert="big_endian", position=trim(loc_pos))

!!!$    !**************************************************************************
!!!$    !** CRAY-Fortran implementation
!!!$    call asnunit (unit,"-N swap_endian",ier)
!!!$    IF (ier /= 0) Then
!!!$       WRITE(*,*)
!!!$       WRITE(*,"(80('='))")
!!!$       WRITE(*,"('EE ',A,T77,' EE')")  'Asign operation on file :'
!!!$       WRITE(*,"('EE ',A          )")  file
!!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")   'Connected to unit :',unit
!!!$       WRITE(*,"('EE ',A,T77,' EE')")   'faild !!'
!!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")'With error flag ',ier
!!!$       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
!!!$       STOP
!!!$    End IF
!!!$
!!!$    Open(unit=unit, file=trim(file), action=trim(action), &
!!!$         status=trim(status), access="stream")

  End Subroutine open_as_big_endian_stream

  !============================================================================
!  !> Subroutine for closing files with big endian encoding
!  !>
  Subroutine close_big_endian_stream(unit)

    integer(kind=4) ,intent(in) :: unit
!   integer(kind=4)             :: ier

   !**************************************************************************
    !** GFortran, Intel implementation
    CLOSE(unit=unit)

!!!$
!!!$    !**************************************************************************
!!!$    !** CRAY-Fortran implementation
!!!$    call asnunit (unit,"-R",ier)
!!!$    IF (ier /= 0) Then
!!!$       WRITE(*,*)
!!!$       WRITE(*,"(80('='))")
!!!$       WRITE(*,"('EE ',A,T77,' EE')")  'Asign release on unit :',unit
!!!$       WRITE(*,"('EE ',A,T77,' EE')")  'faild !!'
!!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")'With error flag ',ier
!!!$       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
!!!$       STOP
!!!$    End IF
!!!$
!!!$    CLOSE(unit=unit)
!
  End Subroutine close_big_endian_stream
  !> @}
!  !# End of memeber group "vtkio private auxilliary routines" #################

End Module vtkwritter
