Module filters

USE vtkwritter
USE ghost
USE constants_module
USE variable_module

IMPLICIT NONE

CONTAINS
!===============================================================================
!===============================================================================
SUBROUTINE contrast(ct,ct_data,ct_new,file_name)

    TYPE(tCt),INTENT(INOUT)                                   ::ct
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::ct_data
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)   ::ct_new
    CHARACTER(c_len),INTENT(IN)                               ::file_name
    CHARACTER(c_len)                                          ::fname


        allocate(ct_new(ct%pixel_x,ct%pixel_y,ct%slices))

           ct_new= ct_data*2.0+50
        
        write(fname,'(A)')'contrast.vtk'
       CALL write_vtk(ct,ct_new,fname)


END SUBROUTINE contrast
!===============================================================================
!===============================================================================
SUBROUTINE sobel_op(ct,ct_data)

    Real(kind=8),dimension(3,3)                               ::s_patch_x
    Real(kind=8),dimension(3,3)                               ::s_patch_y
    TYPE (tCt),INTENT(INOUT)                                  ::ct
    Real(kind=8),dimension(:,:,:),allocatable,INTENT(INOUT)   ::ct_data
    Real(kind=8),dimension(:,:,:),allocatable                 ::temp
    REAL(KIND=8)                                              ::Gx,Gy
    INTEGER                                                   ::ctr
    INTEGER                                                   ::ii,jj
    INTEGER                                                   ::unit=100
    CHARACTER(len=256)                                        ::file_name
    REAL(Kind=rk4),Dimension(3)                               :: spacing,origin
    INTEGER(kind=rk4),Dimension(3)                            :: extend

WRITE(*,*)'========================================================================'
WRITE(*,*)'Dimensions in sobel .....'
WRITE(*,*)'========================================================================'

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'

!CALL ghost_region_2d(ct,ct_data,temp,3)

ctr=1! (len-1)/2
s_patch_x=RESHAPE((/-1,-2,-1,0,0,0,1,2,1/),(/3,3/))
s_patch_y=RESHAPE((/1,0,-1,2,0,-2,1,0,-1/),(/3,3/))

         Do jj=1,ct%pixel_y

                Do ii=1,ct%pixel_x

                Gx=sum(s_patch_x(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))
                Gy=sum(s_patch_y(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))

                ct_data(ii,jj,1)=sqrt(Gx**2+Gy**2)
                END DO

        END DO

write(file_name,'(A)')'sobel.vtk'

open(unit,file=trim(file_name),action='write',access='stream'&
                                 &,status='replace',convert='big_endian')

        spacing=(/ct%dx,ct%dy,ct%dz/)
        origin =(/ct%offset_x ,ct%offset_y ,ct%offset_z/)
        extend =(/ct%pixel_x ,ct%pixel_y ,ct%slices/)

!        call write_vtk_head(unit)
!        call write_vtk_structured_points_head(unit&
!                   &,extend,spacing,origin)

        WRITE(unit)int(ct_data,2)
Close(unit)

END SUBROUTINE
!=============================================================================
!===============================================================================
SUBROUTINE scharr_op(ct,ct_data)

    Real(kind=8),dimension(3,3)                               ::s_patch_x
    Real(kind=8),dimension(3,3)                               ::s_patch_y
    TYPE (tCt),INTENT(INOUT)                                  ::ct
    Real(kind=8),dimension(:,:,:),allocatable,INTENT(INOUT)   ::ct_data
    Real(kind=8),dimension(:,:,:),allocatable                 ::temp
    REAL(KIND=8)                                              ::Gx,Gy
    INTEGER                                                   ::ctr
    INTEGER                                                   ::ii,jj
    INTEGER                                                   ::unit=100
    CHARACTER(len=256)                                        ::file_name
    REAL(Kind=rk4),Dimension(3)                               :: spacing,origin
    INTEGER(kind=rk4),Dimension(3)                            :: extend

WRITE(*,*)'========================================================================'
WRITE(*,*)'Dimensions in scharr_op .....'
WRITE(*,*)'========================================================================'

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'

!CALL ghost_region_2d(ct,ct_data,temp,3)

ctr=1! (len-1)/2
s_patch_x=RESHAPE((/3,0,-3,10,0,-10,3,0,-3/),(/3,3/))
s_patch_y=RESHAPE((/3,10,3,0,0,0,-3,-10,-3/),(/3,3/))

         Do jj=1,ct%pixel_y

                Do ii=1,ct%pixel_x

                Gx=sum(s_patch_x(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))
                Gy=sum(s_patch_y(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))

                ct_data(ii,jj,1)=sqrt(Gx**2+Gy**2)
                END DO

        END DO

write(file_name,'(A)')'scharr.vtk'

open(unit,file=trim(file_name),action='write',access='stream'&
                                 &,status='replace',convert='big_endian')

        spacing=(/ct%dx,ct%dy,ct%dz/)
        origin =(/ct%offset_x ,ct%offset_y ,ct%offset_z/)
        extend =(/ct%pixel_x ,ct%pixel_y ,ct%slices/)

!        call write_vtk_head(unit)
!        call write_vtk_structured_points_head(unit&
!                   &,extend,spacing,origin)

        WRITE(unit)int(ct_data,2)
Close(unit)

END SUBROUTINE
!=============================================================================
!===============================================================================
SUBROUTINE prewitt_op(ct,ct_data)

    Real(kind=8),dimension(3,3)                               ::s_patch_x
    Real(kind=8),dimension(3,3)                               ::s_patch_y
    TYPE (tCt),INTENT(INOUT)                                  ::ct
    Real(kind=8),dimension(:,:,:),allocatable,INTENT(INOUT)   ::ct_data
    Real(kind=8),dimension(:,:,:),allocatable                 ::temp
    REAL(KIND=8)                                              ::Gx,Gy
    INTEGER                                                   ::ctr
    INTEGER                                                   ::ii,jj
    INTEGER                                                   ::unit=100
    CHARACTER(len=256)                                        ::file_name
    REAL(Kind=rk4),Dimension(3)                               :: spacing,origin
    INTEGER(kind=rk4),Dimension(3)                            :: extend

WRITE(*,*)'========================================================================'
WRITE(*,*)'Dimensions in prewitt_op.....'
WRITE(*,*)'========================================================================'

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'

!CALL ghost_region_2d(ct,ct_data,temp,3)

ctr=1! (len-1)/2

s_patch_x=RESHAPE((/-1,-1,-1,0,0,0,1,1,1/),(/3,3/))
s_patch_y=RESHAPE((/1,0,-1,1,0,-1,1,0,-1/),(/3,3/))

         Do jj=1,ct%pixel_y

                Do ii=1,ct%pixel_x

                Gx=sum(s_patch_x(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))
                Gy=sum(s_patch_y(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))

                ct_data(ii,jj,1)=sqrt(Gx**2+Gy**2)
                END DO

        END DO

write(file_name,'(A)')'prewitt.vtk'

open(unit,file=trim(file_name),action='write',access='stream'&
                                 &,status='replace',convert='big_endian')

        spacing=(/ct%dx,ct%dy,ct%dz/)
        origin =(/ct%offset_x ,ct%offset_y ,ct%offset_z/)
        extend =(/ct%pixel_x ,ct%pixel_y ,ct%slices/)

!        call write_vtk_head(unit)
!        call write_vtk_structured_points_head(unit&
!                   &,extend,spacing,origin)

        WRITE(unit)int(ct_data,2)
Close(unit)

END SUBROUTINE
!=============================================================================
!===============================================================================
SUBROUTINE laplacian_op(ct,ct_data)

    TYPE (tCt),INTENT(INOUT)                                  ::ct
    Real(kind=8),dimension(:,:,:),allocatable,INTENT(INOUT)   ::ct_data
    Real(kind=8),dimension(:,:,:),allocatable                 ::temp
    REAL(KIND=8),dimension(3,3)                               ::patch
    INTEGER                                                   ::ctr
    INTEGER                                                   ::ii,jj
    INTEGER                                                   ::unit=100
    CHARACTER(len=256)                                        ::file_name
    REAL(Kind=rk4),Dimension(3)                               :: spacing,origin
    INTEGER(kind=rk4),Dimension(3)                            :: extend

WRITE(*,*)'========================================================================'
WRITE(*,*)'Dimensions in laplacian op .....'
WRITE(*,*)'========================================================================'

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'

!CALL ghost_region_2d(ct,ct_data,temp,3)

ctr=1! (len-1)/2

!patch=RESHAPE((/-0.25,-1.0,-0.25,-1.0,4.0,-1.0,-0.25,-1.0,-0.25/),(/3,3/))
patch=RESHAPE((/0.25,1.0,0.25,1.0,-4.0,1.0,0.25,1.0,0.25/),(/3,3/))
!patch=RESHAPE((/0.0,1.0,0.0,1.0,-4.0,1.0,0.0,1.0,0.0/),(/3,3/))


         Do jj=1,ct%pixel_y

                Do ii=1,ct%pixel_x

                ct_data(ii,jj,1)=sum(patch(:,:)*temp(ii-ctr:ii+ctr,jj-ctr:jj+ctr,1))
                ct_data(ii,jj,1)=sqrt(ct_data(ii,jj,1))

                END DO

        END DO

write(file_name,'(A)')'laplacian.vtk'

open(unit,file=trim(file_name),action='write',access='stream'&
                                 &,status='replace',convert='big_endian')

        spacing=(/ct%dx,ct%dy,ct%dz/)
        origin =(/ct%offset_x ,ct%offset_y ,ct%offset_z/)
        extend =(/ct%pixel_x ,ct%pixel_y ,ct%slices/)

!        call write_vtk_head(unit)
!        call write_vtk_structured_points_head(unit&
!                   &,extend,spacing,origin)

        WRITE(unit)int(ct_data,2)
Close(unit)

END SUBROUTINE
!=============================================================================
END MODULE
