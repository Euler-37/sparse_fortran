program main
   use sparse_fortran
   implicit none
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
end program main
