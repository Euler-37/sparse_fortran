program check_coo
   use sparse_fortran
   implicit none
   block
      type(sparse_coo)::a
      integer::i,j
      write(*,*)"Init general"
      call a%sparse_init(10,10,4)
      do i=1,9
         call a%sparse_append(i,i+1,real(i,8))
      end do
      do i=1,size(a)
         write(*,*)a%row_(i),a%col_(i),a%val_(i)
      end do
   end block

   block
      type(sparse_coo)::a
      integer::ix(10),jx(10)
      real(8)::val(10)
      integer::i
      write(*,*)"Init Vector"
      ix=[1,2,1,3,4,5,1,6,8,9]
      jx=[9,10,7,2,4,5,3,8,2,2]
      val=[1,2,3,4,5,6,7,8,9,10]
      call a%sparse_init(ix,jx,val,10,10)
      do i=1,size(a)
         write(*,*)a%row_(i),a%col_(i),a%val_(i)
      end do
   end block

   block
      type(sparse_coo)::a
      real(8)::x(10,10)
      integer::i
      write(*,*)"Init Matrix"
      x=0.d0
      do i=1,9
         x(i,i+1)=i
      end do
      call a%sparse_init(x)
      do i=1,size(a)
         write(*,*)a%row_(i),a%col_(i),a%val_(i)
      end do
   end block

   block
      type(sparse_coo)::a
      real(8)::x(10,10),y(10)
      integer::i
      write(*,*)"Matrix Product"
      x=0.d0
      do i=1,9
         x(i,i+1)=i
      end do
      y=[real::(i,i=1,10)]
      call a%sparse_init(x)
      write(*,"(*(F7.4,1x))")matmul(x,y)
      y=a%sparse_product(y)
      write(*,"(*(F7.4,1x))")y
   end block

   block
      type(sparse_coo)::a
      integer::ix(10),jx(10)
      integer::i
      real(8)::val(10)
      write(*,*)"Init None Square Vec"
      ix=[1,2,1,3,4,5,1,6,8,9]
      jx=[9,10,7,2,4,5,3,8,2,2]
      val=[1,2,3,4,5,6,7,8,9,10]
      call a%sparse_init(ix,jx,val,20,10)
      do i=1,size(a)
         write(*,*)a%row_(i),a%col_(i),a%val_(i)
      end do
   end block

   block
      type(sparse_coo)::a
      real(8)::x(20,10),y(10)
      integer::i
      write(*,*)"Init None Square Array"
      x=0.d0
      do i=1,9
         x(2*i,i+1)=i
      end do
      y=[real::(i,i=1,10)]
      call a%sparse_init(x)
      write(*,"(*(F7.4,1x))")matmul(x,y)
      write(*,"(*(F7.4,1x))")a%sparse_product(y)
   end block

   block
      type(sparse_coo)::a
      integer::ix(10),jx(10)
      integer::i
      real(8)::val(10)
      write(*,*)"Init Vec Order"
      ix =[1,2,1,3,4,5,1,6,8,14]
      jx =[9,10,7,2,4,5,3,8,2,8]
      val=[1,2,3,4,5,6,7,8,9,10]
      call a%sparse_init(ix,jx,val,dimx=20,dimy=10,order="F")
      do i=1,size(a)
         write(*,*)a%row_(i),a%col_(i),a%val_(i)
      end do
   end block

   block
      type(sparse_coo)::a
      real(8)::x(20,10),y(20)
      integer::i
      write(*,*)"Init Square Order"
      x=0.d0
      do i=1,9
         x(2*i,i+1)=i
      end do
      y=[real::(i,i=1,20)]
      call a%sparse_init(x,order="F")
      do i=1,size(a)
         write(*,*)a%row_(i),a%col_(i),a%val_(i)
      end do
      write(*,"(*(F9.4,1x))")matmul(y,x)
      write(*,"(*(F9.4,1x))")a%sparse_product(y)
   end block
end program check_coo
