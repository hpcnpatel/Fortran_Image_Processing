Module interpolate_ctdata

    USE CONSTANTS_MODULE
    USE VARIABLE_MODULE
    USE vtkwritter

IMPLICIT NONE

CONTAINS

!ct_data is the original data here.
!Z1 and Z2 are two slices in Z direction.

! Z1  Z1   Z1
! |   |
! |   +  = Z112
! |   |
! + = Z12  Z12
! |   |
! |   +  = Z122
! |   |
! Z2  Z2   Z2

!===================================================================
!===================================================================
SUBROUTINE in_data_3d_Z_zero(ct,ct_data,ct_interpolate,file_name)

    TYPE(tCt),INTENT(INOUT)                                   ::ct
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::ct_data
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)   ::ct_interpolate
    CHARACTER(c_len),INTENT(IN)                               ::file_name
    CHARACTER(c_len)                                          ::fname
    INTEGER                                                   ::i,j,k

        allocate(ct_interpolate(ct%pixel_x,ct%pixel_y,ct%slices))

!do k=1,ct%slices
!do j=1,ct%pixel_y
!do i=1,ct%pixel_x
!
!   IF(ct_data(i,j,k) > 1250.0_8)then
!        ct_data(i,j,k)=3000.0_8
!   ENDIF
!   IF(ct_data(i,j,k) < 950.0_8)then
!        ct_data(i,j,k)=0.0_8
!   ENDIF
!
!end do
!end do
!end do
     ct_interpolate= ct_data
        
        write(fname,'(A)')'gray_data.vtk'
       CALL write_vtk(ct,ct_interpolate,fname)


END SUBROUTINE in_data_3d_Z_zero
!=================================================================================
!=================================================================================
SUBROUTINE in_data_3d_Z_one(ct,ct_data,ct_interpolate,file_name)

    TYPE(tCt),INTENT(INOUT)                                   ::ct
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::ct_data
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)   ::ct_interpolate
    INTEGER                                                   :: kk
    CHARACTER(c_len),INTENT(IN)                               ::file_name
    CHARACTER(c_len)                                          ::fname

!===========================================================================

        write(fname,'(A)')'gray_data.vtk'
        CALL write_vtk(ct,ct_data,fname)

!=================================================================================

write(*,*)'We increase number of slices in Z direction by a factor of 2&
          &by linear Interpolation'

        allocate(ct_interpolate(ct%pixel_x,ct%pixel_y,ct%slices*2))

        ct_interpolate=0.0

    Do kk=1,ct%slices

        IF(kk == ct%slices)then

            ct_interpolate(:,:,kk*2)=ct_data(:,:,kk)
            ct_interpolate(:,:,kk*2-1)=ct_data(:,:,kk)

        ELSE

            ct_interpolate(:,:,kk*2-1)=ct_data(:,:,kk)
            ct_interpolate(:,:,kk*2)=(ct_data(:,:,kk)+ct_data(:,:,kk+1))/2.0

        END IF        
    End Do

ct%dz = ct%dz*0.5
ct%slices = ct%slices*2 

ct_data=ct_interpolate

write(*,*)'ct_mod bounds in IN_OUT module after interpolating:'
write(*,*)'lower bound',lbound(ct_interpolate)
write(*,*)'upper bound',ubound(ct_interpolate)
write(*,*)'========================================================================'

!below two variables changes because of insertion of new slices via
!interpolation..

END SUBROUTINE in_data_3d_Z_one
!=================================================================================
!Call this subroutine if two slices in z direction is required
SUBROUTINE in_data_3d_Z_two(ct,ct_data,ct_interpolate,file_name)

    TYPE(tCt),INTENT(INOUT)                                   ::ct
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::ct_data
    REAL(KIND=rk8),ALLOCATABLE,DIMENSION(:,:,:),INTENT(OUT)   ::ct_interpolate
    INTEGER                                                   :: kk
    CHARACTER(c_len),INTENT(IN)                               ::file_name
    CHARACTER(c_len)                                          ::fname

!===========================================================================

        write(fname,'(A)')'gray_data.vtk'
        CALL write_vtk(ct,ct_data,fname)

!=================================================================================
write(*,*)'We increase number of slices in Z direction by a factor of 4'
!Linear interpolation in Z direction thus by reducing the distance by 1/4
!between two slices in Z direction.

        allocate(ct_interpolate(ct%pixel_x,ct%pixel_y,ct%slices*4))

        ct_interpolate=0.0

    Do kk=1,ct%slices

        IF (kk == ct%slices) then
            ct_interpolate(:,:,kk*4)=ct_data(:,:,kk)
            ct_interpolate(:,:,kk*4-1)=ct_data(:,:,kk)
            ct_interpolate(:,:,kk*4-2)=ct_data(:,:,kk)
            ct_interpolate(:,:,kk*4-3)=ct_data(:,:,kk)
        ELSE
            ct_interpolate(:,:,kk*4-3)=ct_data(:,:,kk)
            ct_interpolate(:,:,kk*4-1)=(ct_data(:,:,kk)+ct_data(:,:,kk+1))/2.0
            ct_interpolate(:,:,kk*4-2)=(((ct_data(:,:,kk)+ct_data(:,:,kk+1))/2.0)&
            &+ct_data(:,:,kk))/2
            ct_interpolate(:,:,kk*4)=(((ct_data(:,:,kk)+ct_data(:,:,kk+1))/2.0)&
            &+ct_data(:,:,kk+1))/2
        END IF
    End Do

ct%dz = ct%dz*0.25
ct%slices = ct%slices*4
ct_data=ct_interpolate

        !write(fname,'(A)')'gray_data.vtk'
        !CALL write_vtk(ct,ct_data,fname)

write(*,*)'ct_mod bounds in IN_OUT module after interpolating:'
write(*,*)'lower bound',lbound(ct_interpolate)
write(*,*)'upper bound',ubound(ct_interpolate)
write(*,*)'========================================================================'

END SUBROUTINE in_data_3d_Z_two
!=================================================================================
END MODULE interpolate_ctdata
