module sparse_csr_fortran
   use utils_mod
   implicit none
   private
   public::sparse_csr,size
   type::sparse_csr
      integer,allocatable::idx_(:)
      integer,allocatable::offset_(:)
      real(8),allocatable::val_(:)
      integer::dimx_
      integer::dimy_
      integer::num_
      integer::capacity_
      character::order_
   contains
      generic::sparse_init=> sparse_init_none, sparse_init_vec, sparse_init_array
      procedure,pass::sparse_init_none
      procedure,pass::sparse_init_vec
      procedure,pass::sparse_init_array
      procedure,pass::sparse_search
      procedure,pass::sparse_set_capacity
      procedure,pass::sparse_append
      procedure,pass::sparse_val
      procedure,pass::sparse_product
   end type sparse_csr
   real(8),parameter,private::eps=1.0d-15
   integer,parameter,private::max_init=1000
   interface size
      module procedure sparse_size
   end interface size
contains
   function sparse_size(this)result(n)
      type(sparse_csr),intent(in)::this
      integer::n
      n=this%num_
   end function sparse_size

   subroutine sparse_init_none(this,dimx,dimy,capacity,order)
      class(sparse_csr),intent(inout)::this
      integer,intent(in)             ::dimx
      integer,optional,intent(in)    ::dimy
      integer,optional,intent(in)    ::capacity
      character,optional,intent(in)  ::order
      this%dimx_=dimx
      this%dimy_=optval(dimy,dimx)
      this%num_=0
      this%capacity_=optval(capacity,max_init)
      this%order_=optval(order,"C")
      allocate(this%idx_(this%capacity_))
      allocate(this%val_(this%capacity_))
      if(this%order_=="C")then
         allocate(this%offset_(0:this%dimx_),source=0)
      else
         allocate(this%offset_(0:this%dimy_),source=0)
      end if
   end subroutine sparse_init_none

   subroutine sparse_init_array(this,a,order)
      class(sparse_csr),intent(inout)::this
      real(8),intent(in)             ::a(:,:)
      character,optional,intent(in)  ::order
      integer::i,j,k
      this%dimx_ =size(a,dim=1)
      this%dimy_ =size(a,dim=2)
      this%num_  =count(abs(a)>eps)
      this%order_=optval(order,"C")
      this%capacity_=this%num_
      allocate(this%idx_(this%capacity_))
      allocate(this%val_(this%capacity_))
      if(this%order_=="C")then
         allocate(this%offset_(0:this%dimx_),source=0)
      else
         allocate(this%offset_(0:this%dimy_),source=0)
      end if
      if(this%order_=="C")then
         k=0
         this%offset_(0)=0
         do i=1,this%dimx_
            do j=1,this%dimy_
               if(abs(a(i,j))>eps)then
                  k=k+1
                  this%idx_(k)=j
                  this%val_(k)=a(i,j)
               end if
            end do
            this%offset_(i)=k
         end do
      else
         k=0
         this%offset_(0)=0
         do j=1,this%dimy_
            do i=1,this%dimx_
               if(abs(a(i,j))>eps)then
                  k=k+1
                  this%idx_(k)=i
                  this%val_(k)=a(i,j)
               end if
            end do
            this%offset_(j)=k
         end do
      end if
   end subroutine sparse_init_array

   subroutine sparse_init_vec(this,row,col,val,dimx,dimy,order)
      class(sparse_csr),intent(inout)::this
      integer,intent(in)             ::row(:)
      integer,intent(in)             ::col(:)
      real(8),intent(in)             ::val(:)
      integer,intent(in)             ::dimx
      integer,optional,intent(in)    ::dimy
      character,optional,intent(in)  ::order
      integer::i
      call this%sparse_init_none(dimx,dimy,size(val),order)
      do i=1,this%capacity_
         call this%sparse_append(row(i),col(i),val(i))
      end do
   end subroutine sparse_init_vec

   subroutine sparse_append(this,row,col,val)
      class(sparse_csr),intent(inout)::this
      integer,intent(in)::row
      integer,intent(in)::col
      real(8),intent(in)::val
      integer::idx
      logical::judge
      integer::capacity
      integer::i,j
      if(row>this%dimx_)then
         write(*,"(A,1x,g0,1x,A,1x,g0)")"Sparse Append Error: Index",row, &
            &"of dimension 1 of array 'sparse_csr' above upper bound of",this%dimx_
         error stop
      end if
      if(col>this%dimy_)then
         write(*,"(A,1x,g0,1x,A,1x,g0)")"Sparse Append Error: Index",col, &
            &"of dimension 2 of array 'sparse_csr' above upper bound of",this%dimy_
         error stop
      end if
      judge=abs(val)>eps
      idx=this%sparse_search(row,col)
      if(this%order_=="C")then
         i=row
         j=col
      else
         i=col
         j=row
      end if
      if(idx<0.and.judge)then
         ! append
         idx=abs(idx)
         if(this%num_+1>this%capacity_)then
            capacity=this%capacity_+this%capacity_/3
            call this%sparse_set_capacity(capacity)
            this%capacity_=capacity
         end if
         this%idx_(idx+1:this%num_+1)=this%idx_(idx:this%num_)
         this%idx_(idx)=j
         this%val_(idx+1:this%num_+1)=this%val_(idx:this%num_)
         this%val_(idx)=val
         this%offset_(i:)=this%offset_(i:)+1
         this%num_=this%num_+1
      elseif(idx>0.and.judge)then
         this%val_(idx)=val
      elseif(idx>0.and.(.not.judge))then
         this%idx_(idx:this%num_-1)=this%idx_(idx+1:this%num_)
         this%idx_(this%num_)=0
         this%val_(idx:this%num_-1)=this%val_(idx+1:this%num_)
         this%val_(this%num_)=0
         this%offset_(i:)=this%offset_(i:)-1
         this%num_=this%num_-1
      end if
   end subroutine sparse_append

   function sparse_product(this,a)result(b)
      class(sparse_csr),intent(inout)::this
      real(8),intent(in) ::a(:)
      real(8),allocatable::b(:)
      integer::i,j
      if(this%order_=="C")then
         if(size(a,dim=1)==this%dimy_)then
            allocate(b(this%dimx_))
            b=0.d0
            do i=1,this%dimx_
               do j=this%offset_(i-1)+1,this%offset_(i)
                  b(i)=b(i)+this%val_(j)*a(this%idx_(j))
               end do
            end do
         else
            error stop "Different shape for arguments 'vector_a' and 'sparse_dim' "
         end if
      else
         if(size(a,dim=1)==this%dimx_)then
            allocate(b(this%dimy_))
            b=0.d0
            do i=1,this%dimy_
               do j=this%offset_(i-1)+1,this%offset_(i)
                  b(i)=b(i)+this%val_(j)*a(this%idx_(j))
               end do
            end do
         else
            error stop "Different shape for arguments 'vector_a' and 'sparse_dim' "
         end if
      end if
   end function sparse_product

   function sparse_val(this,i,j)result(res)
      class(sparse_csr),intent(inout)::this
      integer,intent(in)::i
      integer,intent(in)::j
      integer::idx
      real(8)::res
      idx=this%sparse_search(i,j)
      if(idx<0)then
         res=0.0d0
      else
         res=this%val_(idx)
      end if
   end function sparse_val

   subroutine sparse_set_capacity(this,capacity)
      class(sparse_csr),intent(inout)::this
      integer,intent(in)::capacity
      integer,allocatable::tmp_idx(:)
      real(8),allocatable::tmp_val(:)
      allocate(tmp_idx(capacity))
      allocate(tmp_val(capacity))
      tmp_idx(1:this%capacity_)=this%idx_
      tmp_val(1:this%capacity_)=this%val_
      deallocate(this%idx_)
      deallocate(this%val_)
      call move_alloc(from=tmp_idx,to=this%idx_)
      call move_alloc(from=tmp_val,to=this%val_)
   end subroutine sparse_set_capacity

   function sparse_search(this,row,col)result(idx)
      class(sparse_csr),intent(inout)::this
      integer,intent(in)::row
      integer,intent(in)::col
      integer::idx,i,j
      integer::idx_min,idx_max,ipos,ipos_new
      if(this%num_==0)then
         idx=-1
         return
      end if
      if(this%order_=="C")then
         i=row
         j=col
      else
         i=col
         j=row
      end if
      if(this%offset_(i)-this%offset_(i-1)==0)then
         idx=-(this%offset_(i)+1)
         return
      else
         idx_min=this%offset_(i-1)+1
         idx_max=this%offset_(i)
      end if
      if(this%idx_(idx_min)==j)then
         idx=idx_min
      elseif(this%idx_(idx_max)==j)then
         idx=idx_max
      elseif(this%idx_(idx_max)<j)then
         idx=-(idx_max+1)
      elseif(this%idx_(idx_min)>j)then
         idx=-(idx_min-1)
      else
         do
            idx=(idx_min+idx_max)/2
            if(idx==idx_min)then
               idx=-(idx+1)
               exit
            end if
            ipos_new=this%idx_(idx)
            if(j==ipos_new)then
               exit
            elseif(j>ipos_new)then
               idx_min=idx
            elseif(j<ipos_new)then
               idx_max=idx
            end if
         end do
      end if
   end function sparse_search
end module sparse_csr_fortran
