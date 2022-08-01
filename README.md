# sparse_fortran
## Init sparse matrix
- general method
``` fortran
type(sparse_coo)::a
call a%sparse_init(dimx=10,dimy=10,capacity=4,order="C")
call a%sparse_append(i,j,val)
```

sparse matrix will automatically increase
- coo vectors

``` fortran
type(sparse_coo)::a
integer::ix(9),jx(9)
real(8)::val(9)
call a%sparse_init(ix,jx,val,dimx=10,dimy=10)
```
- array
``` fortran
type(sparse_coo)::a
real(8)::x(20,10)
call a%sparse_init(x,order="F")

```
## Matrix Product

If calculate `matmul(A,x)`, Using `order="C"`.If calculate `matmul(x,A)`, Using `order="F"`.

``` fortran
real(8)::x(10)
type(sparse_csr)::a
!...
write(*,"(*(F9.4,1x))")a%sparse_product(x)
```
