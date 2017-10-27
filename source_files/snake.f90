module snake

    USE sort_module
USE CONSTANTS_MODULE
USE VARIABLE_MODULE
USE vtkwritter
use OMP_LIB
use ghost
    USE MOD_grad_of_image
!USE vtkio
!use mod_type
!use vandermonde_matrix
!use QR_module

IMPLICIT NONE

CONTAINS

!=============================================================================
SUBROUTINE snake_2d(ct,arr,len,logi,np)

    TYPE (tCt)                                                ::gd_coff
    TYPE (tCt),INTENT(INOUT)                                  ::ct
    Real(kind=rk8),dimension(:,:,:),allocatable,INTENT(INOUT) ::arr
    LOGICAL                                                   ::logi
    CHARACTER(c_len)                                          ::file_name,fname
    REAL(KIND=rk8), dimension(:),allocatable ::x_vec,y_vec,z_vec
    INTEGER,INTENT(IN)                                        ::len
    INTEGER                                                   ::c_points
    INTEGER                                                   ::deg,aa,ctr,a
    INTEGER                                                   ::t1,t2,np
    INTEGER ::ii,jj,xx,yy,tt,iter
    REAL                                                      ::cput1,cput2
    Real(kind=rk8),dimension(:,:,:),allocatable                ::temp_2d,ct_ghost
!        Real(kind=rk8),dimension(:,:,:),allocatable                ::g_curve
!        Real(kind=rk8),dimension(:,:,:),allocatable                ::g_curve_1
!        Real(kind=rk8),dimension(:,:,:),allocatable                ::snake_curve
!        Real (kind=rk8), dimension(:,:),allocatable ::A_mat,A1_mat,A2_mat
!        Real (kind=rk8), dimension(:,:),allocatable                ::B_mat
!        Real (kind=rk8), dimension(:,:),allocatable                ::I_mat
!        Real (kind=rk8), dimension(:),allocatable                  ::patch_1d
        Real (kind=rk8), dimension(:),allocatable                  ::temp
    INTEGER,ALLOCATABLE,DIMENSION(:)                          ::SORT_INDEX
!        REAL (kind=rk8), dimension(:),allocatable                  ::x_mat,y_mat
!        REAL (kind=rk8), dimension(:),allocatable ::x_mat_res,y_mat_res
        !Real (kind=rk8), dimension(:,:)  ,allocatable             ::patch_temp
!        CHARACTER(len=256)                                         ::name_file
        !CHARACTER(len=256)                                        ::fmt
!        real(kind=rk8)        ::alp,beta,gamma
!        INTEGER                                                    ::count
!        REAL(KIND=R_K),parameter                         ::pi=3.14159265
write(*,*)'========================================================================'
write(*,*)'Dimensions in snake.....'
write(*,*)'========================================================================'

        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in X-dir: ',ct%pixel_x,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Y-dir: ',ct%pixel_y,' *****'
        WRITE(*,'(A,I0,A)')'***** Num. of Pixel in Z-dir: ',ct%slices,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in X-dir: ',ct%dx,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Y-dir: ',ct%dy,' *****'
        WRITE(*,'(A,F10.7,A)')'***** Spacing in Z-dir: ',ct%dz,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset X-dir starts: ',ct%offset_x,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Y-dir starts: ',ct%offset_y,' *****'
        WRITE(*,'(A,F15.7,A)')'***** Offset Z-dir starts: ',ct%offset_z,' *****'

write(*,*)'========================================================================'
temp_2d=arr
    CALL grad(CT,temp_2d,.false.,np=1)
    CALL ghost_region_2d(ct,temp_2d,ct_ghost,len)


!========================================
! make a dummy curve
!ct_temp=1000.0
!do j=0,98
!do i=0,98
!ct_temp(150+i:350-i,150+j:350-j,1)=1000+i*5+j*5
!enddo
!enddo
!========================================
!                do deg=0,359,2
!                        aa=aa+1
!i=nint(250.0+(1.0+aa)*cos(deg*pi/180.0))
!ii=nint(250.0-(1.0+aa)*cos(deg*pi/180.0))
!j=nint(250.0+(1.0+aa)*sin(deg*pi/180.0))
!jj=nint(250.0-(1.0+aa)*sin(deg*pi/180.0))
!ct_temp(i:ii,jj:j,1)=4000-i*30
!                end do
!========================================

c_points=20
                allocate(x_vec(c_points))
                allocate(y_vec(c_points))
                allocate(z_vec(c_points))
                z_vec=0.0
                !========================================


                aa=0
                do deg=0,359,18
                    aa=aa+1
                    x_vec(aa)=250.0+50.0*cos(deg*pi/180.0)
                    y_vec(aa)=350.0+45.0*sin(deg*pi/180.0)
                    !!snake_curve(nint(x_vec(aa)),nint(y_vec(aa)),1)=3000.0
                    !arr(nint(x_vec(aa)),nint(y_vec(aa)),1)=4000.0
                end do

allocate(temp(len**2))
allocate(SORT_INDEX(len**2))

ctr=(len-1)/2

!**********************************************
do iter=1,5
    write(*,*)'iteration starts-------',iter
do tt=1,c_points

    ii=nint(x_vec(tt))
    jj=nint(y_vec(tt))

       aa=0
        do yy=-ctr,ctr
        do xx=-ctr,ctr
            aa=aa+1
            temp(aa)=ct_ghost(ii+xx,jj+yy,1)
        enddo
        enddo
        call SORTRX_real8(temp(1:aa),SORT_INDEX)
        a=0
        do yy=-ctr,ctr
        do xx=-ctr,ctr
            a=a+1
            if(SORT_INDEX(aa)==a)then
            ii=ii+xx 
            jj=jj+yy 
            endif
        enddo
        enddo
        x_vec(tt)=ii
        y_vec(tt)=jj
!        if(iter==5)ct_ghost(ii,jj,1)=500 
!            write(*,*)ii,',',jj,',',500
enddo
enddo
!**********************************************
arr(:,:,1)=ct_ghost(1:ct%pixel_x,1:ct%pixel_y,1)
!IF (logi .eqv. .true.)then
!       write(fname,'(A)')'snake.vtk'
!       CALL write_vtk(ct,arr,fname)
!ENDIF
IF (logi .eqv. .true.)then
       write(fname,'(A)')'snake.vtk'
       CALL write_vtk_polydata(fname,x_vec,y_vec,z_vec)
ENDIF

END subroutine
!=============================================================================
end module snake
