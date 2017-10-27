!===========================================================
        CALL vandermonde_3d(V_arr,len)

        dim_mm=size(V_arr,2)
        dim_nn=size(V_arr,1)
        temp_V=V_arr

        IF(allocated(cc))deallocate(cc)
        IF(allocated(dd))deallocate(dd)
        allocate(cc(1:dim_mm))
        allocate(dd(1:dim_mm))

       CALL QR_dcmp_non_square(temp_V,dim_nn,dim_mm,cc,dd,sing)

        IF(allocated(temp_1d))deallocate(temp_1d)
        allocate(temp_1d(dim_nn))
                temp_1d=0.0
                temp_1d(1:dim_mm)=V_arr((dim_nn+1)/2,1:dim_mm)

        IF(allocated(patch_1d))deallocate(patch_1d)
        allocate(patch_1d(1:dim_nn))

        CALL QR_solv_non_square(temp_V,dim_nn,dim_mm,cc,dd,temp_1d,patch_1d)
!===========================================================

        ctr=(len-1)/2

DO ll=1,4
IF(ll == 1)GV=min
IF(ll == 2)GV=CV3
IF(ll == 3)GV=CV2
IF(ll == 4)GV=CV1

n
call cpu_time(cput1)
call system_clock(t1)

        num_threads=np
        call omp_set_num_threads(num_threads)
        gl_length = ct%slices
        local_length = gl_length/num_threads

!$OMP PARALLEL private(thread_num,kk_start,kk_end,ii,jj,kk,i,j,k,add,a)

        thread_num=omp_get_thread_num()
        kk_start=thread_num*local_length+1
        kk_end=kk_start+local_length-1

        if(thread_num .eq. num_threads-1) kk_end=gl_length

        write(*,'(5(A,I0))')"thread_num:",thread_num, &
                            "; num_threads:",omp_get_num_threads(),&
                            "; kk_start:",kk_start,&
                            "; kk_end:", kk_end,&
                            "; last should be:",ct%slices
DO kk=kk_start,kk_end

    if(mod(kk,10) .eq. 0 .and. thread_num .eq. 0) then
    write(*,'(3(A,I0))') "thread_num:",thread_num, &
    "; kk:",kk, " from ",kk_end
    end if

!        a=0
!        do k=-ctr,ctr
!        do j=-ctr+abs(k),ctr-abs(k)
!        do i=-ctr+abs(j+k),ctr-abs(j+k)
!             a=a+1
!           patch_3d(i,j,k)=patch_1d(a)
!        end do
!        end do
!        end do


!Do kk=1,ct%slices
Do jj=1,ct%pixel_y
Do ii=1,ct%pixel_x

    IF(ct_grad(ii,jj,kk) .gt. GV)then

add=0.0
a=0
        do k=-ctr,ctr
        do j=-ctr+abs(k),ctr-abs(k)
        do i=-ctr+abs(j+k),ctr-abs(j+k)
             a=a+1
add=add+patch_1d(a)*ct_ghost(ii+i,jj+j,kk+k)
        end do
        end do
        end do
ct_data(ii,jj,kk)=add+((1-ct_grad(ii,jj,kk)/max)*(add-ct_data(ii,jj,kk)))

    ENDIF

END DO
END DO
END DO


