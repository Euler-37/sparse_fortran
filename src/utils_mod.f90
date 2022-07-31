module utils_mod
   implicit none
   interface optval
      module procedure optval_int
      module procedure optval_char
   end interface optval
contains
   function optval_int(i,default)result(res)
      integer,optional,intent(in)::i
      integer,intent(in)::default
      integer::res
      if(present(i))then
         res=i
      else
         res=default
      end if
   end function optval_int

   function optval_char(s,default)result(res)
      character,optional,intent(in)::s
      character,intent(in)::default
      character::res
      if(present(s))then
         res=s
      else
         res=default
      end if
   end function optval_char
end module utils_mod
