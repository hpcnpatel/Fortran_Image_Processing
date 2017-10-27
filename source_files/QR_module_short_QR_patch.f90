module QR_module_patch
use constants_module
use utilities_module
! http://cc.oulu.fi/~tf/tiedostot/pub/nrf/
!  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#&lt;0(9p#3.
!
! Gesamtaufwand  m*n^2 - n^3/3
!
implicit none
!!      integer,parameter :: rk=8

type QR_matrix_type
   real(kind=rk),allocatable,dimension(:) :: cc
   real(kind=rk),allocatable,dimension(:) :: dd
   integer,allocatable,dimension(:)       :: index
   logical                                :: sing
   integer                                :: max_ind
end type QR_matrix_type

      contains
!************************************************************************************
      subroutine QR_dcmp(a,n,np,cc,dd,sing)

! logical sing
! Constructs the QR decomposition of a(1:n,1:n), with physical dimension np. The upper
! triangular matrix R is returned in the upper triangle of a, except for the diagonal elements
! of R which are returned in dd(1:n). The orthogonal matrix Q is represented as a product of
! n-1 Householder matrices Q_1 . . .Q_{n-1}, where Q_j = 1 - u_j x u_j/c_j. The i-th component
! of u_j is zero for i = 1, . . . , j - 1 while the nonzero components are returned in a(i,j) for
! i = j, . . . ,n. sing returns as true if singularity is encountered during the decomposition,
! but the decomposition is still completed in this case.

      integer           :: n,np
      real(kind=rk)     :: a(np,np),cc(n),dd(n)
      logical           :: sing
      integer           :: i,j,k
      real(kind=rk)     :: scale,sigma,sum,tau



      sing=.false.
do_main:      do  k=1,n-1

          scale=0._rk
        do_column_scale: do i=k,n
          scale=max(scale,abs(a(i,k)))
         ! va(:)=a(:,k)   spielt in allen Schleifen eine Rolle, unterliegt aber Reduktionen
        enddo do_column_scale

        if(scale.eq.0._rk) then
          sing=.true.
          cc(k)=0._rk
          dd(k)=0._rk
        else

               ! do 12 und do 13 koennen zusammengefasst werden
          sum=0._rk
          do_scaling_and_norm: do i=k,n
            a(i,k)=a(i,k)/scale
            sum=sum+a(i,k)**2
          enddo do_scaling_and_norm

          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          cc(k)=sigma*a(k,k)        ! nur das Inverse wird spaeter verwendet
          dd(k)=-scale*sigma        ! nur das Inverse wird spaeter verwendet

             ! zeitintensiver Teil
          do_rotation:  do j=k+1,n

            sum=0._rk
            do_dot_product: do i=k,n
              sum=sum+a(i,k)*a(i,j)    ! a(i,k) nur im unteren Dreieck mit Diagonale, a(i,j) nur im echten oberen Dreieck
            enddo do_dot_product

            tau=sum/cc(k)

            do_projection: do i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
            enddo do_projection

          enddo do_rotation

        endif
     enddo do_main

      dd(n)=a(n,n)

      if(dd(n).eq.0._rk) sing=.true.

      end subroutine QR_dcmp
!************************************************************************************
      subroutine QR_solv(a,n,np,cc,dd,b)

! Solves the set of n linear equations Ax = b, where a is a matrix with physical dimension np.
! a, c, and d are input as the output of the routine QR_dcmp and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer     ::  n,np
      real(kind=rk)     :: a(np,np),b(n),cc(n),dd(n)
      integer     ::  i,j
      real(kind=rk)     :: sum,tau

!  form Q^T b
      do 13 j=1,n-1

        sum=0._rk
        do 11 i=j,n
          sum=sum+a(i,j)*b(i)
11      continue

        tau=sum/cc(j)

        do 12 i=j,n
          b(i)=b(i)-tau*a(i,j)
12      continue

13    continue

! solve R^(-1) b
      call QR_R_solv(a,n,np,dd,b)


      end subroutine QR_solv
!************************************************************************************
      subroutine QR_updt(r,qt,n,np,u,v)

! Given the QR decomposition of some n x n matrix, calculates the QR decomposition of
! the matrix Q (R + u x v). The matrices r and qt have physical dimension np. Note that
! QT is input and returned in qt.

      integer     ::  n,np
      real(kind=rk)     :: r(np,np),qt(np,np),u(np),v(np)
!U    USES Jacobi_rotate
      integer     ::  i,j,k

      do 11 k=n,1,-1
        if(u(k).ne.0._rk)goto 1
11    continue

      k=1
1     do 12 i=k-1,1,-1
        call Jacobi_rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0._rk)then
          u(i)=abs(u(i+1))
        else if(abs(u(i)).gt.abs(u(i+1)))then
          u(i)=abs(u(i))*sqrt(1._rk+(u(i+1)/u(i))**2)
        else
          u(i)=abs(u(i+1))*sqrt(1._rk+(u(i)/u(i+1))**2)
        endif
12    continue

      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue

      do 14 i=1,k-1
        call Jacobi_rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue


      end subroutine QR_updt
!************************************************************************************
      subroutine Jacobi_rotate(r,qt,n,np,i,a,b)

! Given n x n matrices r and qt of physical dimension np, carry out a Jacobi rotation on rows i
! and i+1 of each matrix. a and b are the parameters of the rotation: cos theta = a/sqrt(a^2 + b^2),
! sin theta= b/sqrt(a^2 + b^2).

      integer     ::  n,np,i
      real(kind=rk)     :: a,b,r(np,np),qt(np,np)
      integer     ::  j
      real(kind=rk)     :: c,fact,s,w,y

      if(a.eq.0._rk)then
        c=0._rk
        s=sign(1._rk,b)
      else if(abs(a).gt.abs(b))then
        fact=b/a
        c=sign(1._rk/sqrt(1._rk+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1._rk/sqrt(1._rk+fact**2),b)
        c=fact*s
      endif

      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i  ,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue

      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i  ,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue


      end subroutine Jacobi_rotate
!************************************************************************************
      subroutine QR_R_solv(a,n,np,dd,b)

! real a(np,np),b(n),dd(n)
! Solves the set of n linear equations R x = b, where R is an upper triangular matrix stored
! in a and d. a and d are input as the output of the routine QR_dcmp and are not modified.
! b(1:n) is input as the right-hand side vector, and is overwritten with the solution vector
! on output.

! R =
! d(1) a(1,2) a(1,3) a(1,4) ....
!      d(2)   a(2,3) a(2,4) ....
!             d(3)   a(3,4) ....
!                    d(4)   ....
!                           ....

      integer     ::  n,np
      real(kind=rk)     :: a(np,np),b(n),dd(n)
      integer     ::  i,j
      real(kind=rk)     :: sum

      b(n)=b(n)/dd(n)

      do 12 i=n-1,1,-1

        sum=0._rk
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue

        b(i)=(b(i)-sum)/dd(i)       ! versagt, wenn dd(i)=0
12    continue


      end subroutine QR_R_solv
!************************************************************************************
      subroutine QR_dcmp_non_square(a,nn,mm,cc,dd,sing)

! logical sing
! Constructs the QR decomposition of a(1:nn,1:mm), with physical dimension np. The upper
! triangular matrix R is returned in the upper triangle of a, except for the diagonal elements
! of R which are returned in dd(1:nn). The orthogonal matrix Q is represented as a product of
! nn-1 Householder matrices Q_1 . . .Q_{nn-1}, where Q_j = 1 - u_j x u_j/c_j. The i-th component
! of u_j is zero for i = 1, . . . , j - 1 while the nonzero components are returned in a(i,j) for
! i = j, . . . ,nn. sing returns as true if singularity is encountered during the decomposition,
! but the decomposition is still completed in this case.

! remark that Q is constructed by a sequence of Householder transformations
! Every transformation remains sparse.

      integer,intent(in)                           :: nn,mm
      real(kind=rk),intent(inout)                  :: a(:,:)
      real(kind=rk) ,intent(out)                   :: cc(:)
      real(kind=rk),intent(out)                    :: dd(:)
      logical                                      :: sing
      integer                                      :: i,j,k
      real(kind=rk)                                :: scale,sigma,sum,tau

      if(mm > nn) then
         print*,'QR_dcmp_non_square: error'
         print*,'mm > nn not allowed'
         print*,'mm=',mm
         print*,'nn=',nn
         stop
      endif

      sing=.false.
! we will have n-1 Householder vectors in every case: nn < mm , nn=mm, nn > mm

! possible, that intermediate columns disappear
!original:      do 17 k=1,mm-1
     do_main: do k=1,mm   ! we also allow for the last column

! Householder transforms will be formed
        scale=0._rk
       do_column_scale: do i=k,nn
          scale=max(scale,abs(a(i,k)))
       enddo do_column_scale

        if(scale.eq.0._rk)then
          sing=.true.
          cc(k)=0._rk
          dd(k)=0._rk
        else

! these are the Householder vectors; they are only scaled
          sum=0._rk
do_scaling_and_norm: do i=k,nn
          a(i,k)=a(i,k)/scale          !ation
          sum=sum+a(i,k)**2
          enddo do_scaling_and_norm
          sigma=sign(sqrt(sum),a(k,k))   != ||a|| a_1/|a_1|


! this is the sole part of the Householder vector which is changed
       a(k,k)=a(k,k)+sigma            != ( a_1 + ||a|| a_1/|a_1| )
       cc(k)=sigma*a(k,k)              != ||a|| * ( |a_1| + ||a|| ) ; koennte auch das Inverse abspeichern
       dd(k)=-scale*sigma              != - scale* ||a|| a_1/|a_1| ; verschwindet, wenn a=0
! end of forming Householder transforms
! cc() wird nur zu Division verwandt!
! dd() wird nur zu Division verwandt!

          ! part where R is formed
          do_rotation: do j=k+1,mm                 ! will not be executed for k=mm
          ! Hier auch eine Auswahl von Spalten, auf die die Transformation angewandt werden muss
          ! calculation of a^T a(:,j)
            sum=0._rk
            do_dot_product: do i=k,nn
              sum=sum+a(i,k)*a(i,j)
            enddo do_dot_product

            tau=sum/cc(k)

!calculation of a(:,j) = a(:,j) - tau*a
            do_projection: do i=k,nn
              a(i,j)=a(i,j)-tau*a(i,k)
            enddo do_projection
          enddo do_rotation
        endif
      enddo do_main

      end subroutine QR_dcmp_non_square
!************************************************************************************
      subroutine QR_QTxb_non_square(a,nn,mm,cc,b)

! calculates Q^T b
! a, c, and d are input as the output of the routine QR_dcmp_non_square and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer           :: nn,mm
      real(kind=rk)     :: a(:,:),b(:),cc(:)
      integer           :: i,j
      real(kind=rk)     :: sum,tau

      if(size(b) < nn) then
         write(*,*) 'QR_QTxb_non_square: error; size(b) < nn'
         write(*,'(a,i0)') 'size(b)=',size(b)
         write(*,'(a,i0)') 'nn     =',nn
         stop 'QR_QTxb_non_square: error; size(b) < nn'
      endif

!  form Q^T b
      do 13 j=1,mm

        sum=0._rk
        do 11 i=j,nn
          sum=sum+a(i,j)*b(i)
11      continue

!print*,'cc(j),dd(j)=',cc(j),dd(j)
! Diese Unterscheidung ist neu aufgenommen
        if(abs(cc(j)) == 0._rk) then
           tau=0._rk
        else
           tau=sum/cc(j)
        endif

        do 12 i=j,nn
          b(i)=b(i)-tau*a(i,j)
12      continue

13    continue

      end subroutine QR_QTxb_non_square
!************************************************************************************
      subroutine QR_Qxb_non_square(a,nn,mm,cc,b)

! calculates Q b
! a, c, and d are input as the output of the routine QR_dcmp_non_square and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer           :: nn,mm
      real(kind=rk)     :: a(:,:),b(:),cc(:)
      integer           :: i,j
      real(kind=rk)     :: sum,tau

      ! to ensure that b has the correct values
      ! b has to be large enough!
      do i=mm+1,nn
         b(i)=0._rk
      enddo

      if(size(b) < nn) then
         write(*,*) 'QR_Qxb_non_square: error; size(b) < nn'
         write(*,'(a,i0)') 'size(b)=',size(b)
         write(*,'(a,i0)') 'nn     =',nn
         stop 'QR_Qxb_non_square: error; size(b) < nn'
      endif

!  form Q b
      do 13 j=mm,1,-1  ! in contrast to QR_QTxb_non_square
!!do 13 j=1,mm ! as in QR_QTxb_non_square

        sum=0._rk
        do 11 i=j,nn
          sum=sum+a(i,j)*b(i)
11      continue

!print*,'cc(j),dd(j)=',cc(j),dd(j)
! Diese Unterscheidung ist neu aufgenommen
        if(abs(cc(j)) == 0._rk) then
           tau=0._rk
        else
           tau=sum/cc(j)
        endif

        do 12 i=j,nn
          b(i)=b(i)-tau*a(i,j)
12      continue

13    continue

      end subroutine QR_Qxb_non_square
!************************************************************************************
      subroutine QR_R_solv_non_square(a,mm,dd,b)

! real a(np,np),b(n),dd(n)
! Solves the set of n linear equations R x = b, where R is an upper triangular matrix stored
! in a and d. a and d are input as the output of the routine QR_dcmp and are not modified.
! b(1:n) is input as the right-hand side vector, and is overwritten with the solution vector
! on output.

! R =
! d(1) a(1,2) a(1,3) a(1,4) ....
!      d(2)   a(2,3) a(2,4) ....
!             d(3)   a(3,4) ....
!                    d(4)   ....
!                           ....

      integer           ::  mm
      real(kind=rk)     :: a(:,:),b(:),dd(:)
      integer           ::  i,j
      real(kind=rk)     :: sum

      if(size(b) < mm) then
         write(*,*) 'QR_R_solv_non_square: error; size(b) < mm'
         write(*,'(a,i0)') 'size(b)=',size(b)
         write(*,'(a,i0)') 'mm     =',mm
         stop 'QR_R_solv_non_square: error; size(b) < mm'
      endif

      b(mm)=b(mm)/dd(mm)

      do i=mm-1,1,-1

        sum=0._rk
        do j=i+1,mm
          sum=sum+a(i,j)*b(j)
        enddo

        b(i)=(b(i)-sum)/dd(i)       ! versagt, wenn dd(i)=0
        enddo


      end subroutine QR_R_solv_non_square
!************************************************************************************
      subroutine QR_RT_solv_non_square(a,mm,dd,b)

! real a(np,np),b(n),dd(n)
! Solves the set of n linear equations R^T x = b, where R is an upper triangular matrix stored
! in a and d. a and d are input as the output of the routine QR_dcmp and are not modified.
! b(1:n) is input as the right-hand side vector, and is overwritten with the solution vector
! on output.

! R =
! d(1) a(1,2) a(1,3) a(1,4) ....
!      d(2)   a(2,3) a(2,4) ....
!             d(3)   a(3,4) ....
!                    d(4)   ....
!                           ....

      integer           ::  mm
      real(kind=rk)     :: a(:,:),b(:),dd(:)
      integer           ::  i,j
      real(kind=rk)     :: sum

      if(size(b) < mm) then
         write(*,*) 'QR_RT_solv_non_square: error; size(b) < mm'
         write(*,'(a,i0)') 'size(b)=',size(b)
         write(*,'(a,i0)') 'mm     =',mm
         stop 'QR_RT_solv_non_square: error; size(b) < mm'
      endif

      b(1)=b(1)/dd(1) ! we start with 1 instead of m in QR_R_solv_non_square

      do i=2,mm ! forward instead of backward in QR_R_solv_non_square

        sum=0._rk
        do  j=1,i-1 ! instead of fo j=i+1,mm in QR_R_solv_non_square
          sum=sum+a(j,i)*b(j)  ! instead of sum=sum+a(i,j)*b(j) in QR_R_solv_non_square
        enddo

        b(i)=(b(i)-sum)/dd(i)       ! versagt, wenn dd(i)=0
        enddo


      end subroutine QR_RT_solv_non_square
!************************************************************************************
      subroutine QR_solv_non_square_FOR_PATCH(a,n,mm,cc,dd,b,result)
!Edited by NISARG on 11/08/2015
!To solve the equations to get the QR PATCH
! Solves the set of n linear equations Ax = b, where a is a matrix with physical dimension np.
! a, c, and d are input as the output of the routine QR_dcmp and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer           ::  n,mm
      real(kind=rk)     :: a(:,:),b(:),cc(:),dd(:)
      real(kind=rk)     :: result(:)
      real(kind=rk),allocatable :: temp(:)

!!  form Q^T b
!      temp=b
!      call QR_Qxb_non_square(a,n,mm,cc,temp)

!! solve R^(-1) b
!      call QR_RT_solv_non_square(a,mm,dd,temp)
!      result(1:mm)=temp(1:mm)

      temp=b

                call QR_RT_solv_non_square(a,mm,dd,temp)

                call QR_Qxb_non_square(a,n,mm,cc,temp)

                result(1:n)=temp(1:n)

      end subroutine QR_solv_non_square_FOR_PATCH
!************************************************************************************
      subroutine QR_solv_non_square(a,n,mm,cc,dd,b,result)

! Solves the set of n linear equations Ax = b, where a is a matrix with physical dimension np.
! a, c, and d are input as the output of the routine QR_dcmp and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer           ::  n,mm
      real(kind=rk)     :: a(:,:),b(:),cc(:),dd(:)
      real(kind=rk)     :: result(:)
      real(kind=rk),allocatable :: temp(:)
      !real(kind=rk)     :: sum,tau

!  form Q^T b
      temp=b
      call QR_QTxb_non_square(a,n,mm,cc,temp)

! solve R^(-1) b
      call QR_R_solv_non_square(a,mm,dd,temp)
      result(1:mm)=temp(1:mm)


      end subroutine QR_solv_non_square
!************************************************************************************
      subroutine QR_QTxb_with_selection(a,nn,QR,b)
      !Q^T b
! a, c, and d are input as the output of the routine QR_dcmp and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer              :: nn
      real(kind=rk)        :: a(:,:)
      real(kind=rk)        :: b(:)
      integer              :: i,j,col
      real(kind=rk)        :: sum,tau
      type(QR_matrix_type) :: QR

      ! Q = Q_index(1) Q_index(2) Q_index(3) Q_index(4) .... Q_index(max_ind)
      ! Q_i = I - a 1/cc a^T   symmetrisch
      ! Q^T = Q_index(max_ind) ... Q_index(1)
 !!print*,'QR_QTxb_with_selection: QR%max_ind=',QR%max_ind
!  form Q^T b
do j=1,QR%max_ind

        col=QR%index(j)

        sum=0._rk
        do i=j,nn
          sum=sum+a(i,col)*b(i)
        enddo

!!print*,'QR_QTxb_with_selection: QR%cc(j)=',QR%cc(j)

! Diese Unterscheidung ist neu aufgenommen
        if(abs(QR%cc(j)) == 0._rk) then
           tau=0._rk
        else
           tau=sum/QR%cc(j)
        endif

        do i=j,nn  ! oder QR%max_ind?
          b(i)=b(i)-tau*a(i,col)
        enddo

enddo



      end subroutine QR_QTxb_with_selection
!************************************************************************************
      subroutine QR_Qxb_with_selection(a,nn,QR,b)

! calculates Q^T b
! a, c, and d are input as the output of the routine QR_dcmp and are not modified. b(1:n)
! is input as the right-hand side vector, and is overwritten with the solution vector on output.

      integer              :: nn
      real(kind=rk)        :: a(:,:)
      real(kind=rk)        :: b(:)
      integer              :: i,j,col
      real(kind=rk)        :: sum,tau
      type(QR_matrix_type) :: QR

      ! Q = Q_index(1) Q_index(2) Q_index(3) Q_index(4) .... Q_index(max_ind)
      ! Q_i = I - a 1/cc a^T   symmetrisch
      ! Q^T = Q_index(max_ind) ... Q_index(1)
 !!print*,'QR_Qxb_with_selection: QR%max_ind=',QR%max_ind
!  form Q b
do j=QR%max_ind,1,-1   ! sole difference

        col=QR%index(j)

        sum=0._rk
        do i=j,nn
          sum=sum+a(i,col)*b(i)
        enddo

!!print*,'QR_Qxb_with_selection: QR%cc(j)=',QR%cc(j)

! Diese Unterscheidung ist neu aufgenommen
        if(abs(QR%cc(j)) == 0._rk) then
           tau=0._rk
        else
           tau=sum/QR%cc(j)
        endif

        do i=j,nn  ! oder QR%max_ind?
          b(i)=b(i)-tau*a(i,col)
        enddo

enddo



      end subroutine QR_Qxb_with_selection
!************************************************************************************
      subroutine QR_test
      real(kind=rk),allocatable,dimension(:)   :: b,cc,dd
      real(kind=rk),allocatable,dimension(:)   :: vec
      real(kind=rk),allocatable,dimension(:)   :: ww
      real(kind=rk),allocatable,dimension(:,:) :: a,a_mod,a_orig
      real(kind=rk),allocatable,dimension(:,:) :: result
      logical                                  :: sing
      integer                                  :: nmax,mmax,np
      integer                                  :: nn,mm

       nmax=5;mmax=3
       np=nmax
      allocate(b(max(nmax,mmax)))
      allocate(cc(max(nmax,mmax)))
      allocate(dd(max(nmax,mmax)))
      allocate(a(nmax,mmax))
      allocate(a_orig(nmax,mmax))
      allocate(a_mod(nmax,mmax))
      allocate(result(mmax,mmax))

      do nn=1,nmax
         do mm=1,mmax
            a(nn,mm)=1._rk+ 1._rk/(nn+mm)
          !  a(nn,mm)=nn - mm
          !  a(nn,mm)=0._rk
          !  if(nn <= mm) a(nn,mm)=1._rk
         enddo
      enddo
       !  a(2,2)=0._rk
       !  a(4,4)=0._rk
         do mm=1,mmax
            b(mm)=1._rk
         enddo

       a_orig=a

      a_mod=a

      ! eigentliche Q*R=a_mod  Zerlegung
      ! Input                   A  nmax X mmax Matrix
      call QR_dcmp_non_square(a_mod,nmax,mmax,cc,dd,sing)
      ! Output                 Q R

      do mm=1,nmax
         vec=[(0._rk,nn=1,nmax)]
         vec(mm)=1._rk
         call print('vorher vec_'//trim(string_of(mm)),vec)
         call QR_Qxb_non_square(a_mod,nmax,mmax,cc,vec)
         call QR_QTxb_non_square(a_mod,nmax,mmax,cc,vec)
         call print('vec_'//trim(string_of(mm)),vec)
      enddo

         allocate(ww(max(mmax,nmax)))
      do mm=1,nmax
         vec=[(0._rk,nn=1,nmax)]
         vec(mm)=1._rk
         call QR_Qxb_non_square(a_mod,nmax,mmax,cc,vec)
         ww=matmul(transpose(a_orig),vec)
         call print('vorher ww_'//trim(string_of(mm)),ww)
         call QR_RT_solv_non_square(a_mod,mmax,dd,ww)
         call print('base unit ww_'//trim(string_of(mm)),ww)
      enddo
    !  call QR_dcmp(a_mod,nmax,np,cc,dd,sing)

     !!    call QR_RT_solv_non_square(a_mod,mmax,dd,xx)
     !!    call QR_Qxb_non_square(a_mod,nmax,mmax,cc,xx)
     !!    yy=xx

      write(*,*)
      write(*,*) 'triangular matrix R including diagonal dd'
               write(*,fmt='(a)',advance='yes') repeat('_',mmax*11)
      do nn=1,nmax
         do mm=1,mmax
            if(mm < nn ) then
               write(*,fmt='(a11  )',advance='no') ' . '
            elseif(mm == nn ) then
               write(*,fmt='(e11.3)',advance='no') dd(nn)
            else
               write(*,fmt='(e11.3)',advance='no') a_mod(nn,mm)
            endif
         enddo
               write(*,fmt='(a)',advance='yes')
               if(nn==mmax) then
               write(*,fmt='(a)',advance='yes') repeat('-',mmax*11)
               endif
      enddo
               write(*,fmt='(a)',advance='yes') repeat('_',mmax*11)

       !  call QR_QTxb_non_square(a_mod,nmax,mmax,cc,dd,b)
 !!     do mm=1,mmax
 !!        call QR_QTxb_non_square(a_mod,nmax,mmax,cc,dd,a(:,mm))
 !!     enddo

      ! Test: Berechnung von R^{-1} Q^T A oder äquivalent ((A^T A)^{-1} A^T) A
      do mm=1,mmax
         !                       Q R                     b      (A^T A)^{-1} A^T b = R_1^{-1} Q^T b
         call QR_solv_non_square(a_mod,nmax,mmax,cc,dd,a(:,mm),result(:,mm))
      enddo

      write(*,*)
      write(*,*) 'original matrix'
               write(*,fmt='(a)',advance='yes') repeat('_',mmax*11)
      do nn=1,nmax
         do mm=1,mmax
              if(abs(a_orig(nn,mm)) < 1.e-7) then
                  write(*,fmt='(a11  )',advance='no') ' - '
              else
                  write(*,fmt='(e11.3)',advance='no') a_orig(nn,mm)
              endif
         enddo
               write(*,fmt='(a)',advance='yes')
      enddo
               write(*,fmt='(a)',advance='yes') repeat('_',mmax*11)

      write(*,*)
      write(*,*) 'after calculation: the following matrix must be the identity matrix '
               write(*,fmt='(a)',advance='yes') repeat('_',mmax*11)
      do nn=1,mmax
         do mm=1,mmax
              if(abs(result(nn,mm)) < 1.e-7) then
                  write(*,fmt='(a11  )',advance='no') ' - '
              else
                  write(*,fmt='(e11.3)',advance='no') result(nn,mm)
              endif
         enddo
               write(*,fmt='(a)',advance='yes')
      enddo
               write(*,fmt='(a)',advance='yes') repeat('_',mmax*11)

      end subroutine QR_test
!************************************************************************************
      end module
!************************************************************************************
!    program test_QR_module
!     use QR_module
!     call QR_test
!    end program test_QR_module
