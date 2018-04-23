module sistemas
    contains
        subroutine grad_conj(A,b,x)
            implicit none
            real(8), intent(in) :: A(:,:), b(:)
            real(8), intent(inout) :: x(:)
            real(8), allocatable :: r(:), v(:), r1(:), v1(:), t(:,:)
            real(8) :: tao, s
            integer :: n,i !iteraciones
            n=size(x)
            n=size(b)
            n=size(A,1)
            allocate(r(n))
            allocate(r1(n))
            allocate(v(n))
            allocate(v1(n))
            allocate(t(n,1))
            r1=matmul(A,x)-b
            v1=-r1
            do i=1,20
                t= reshape(v1,[n,1])
                tao = (norma(r1)**2)/norma(matmul(matmul(A,v1),t))
                x = x + tao*v1 
                r1 = r 
                r1 = r + tao*matmul(A,v1)
                if (norma(r)<norma(b)) EXIT
                s=(norma(r1)**2)/(norma(r)**2)
                v1=v
                v1=-r1 + s*v 
            end do
        end subroutine grad_conj
        function norma (v)
                    implicit none
                    integer:: n,i
                    real(8):: v(:), norma
                n=size(v)
                    norma=v(1)**2
                    do i=2,n
                        norma=v(i)**2+norma
                    end do
                    norma=sqrt(norma)
        end function norma

end module
