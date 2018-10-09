module reorder_mod


contains

    subroutine reorder( n, nn, A )
        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: nn
        integer, intent(inout) :: A( n )

        integer :: i, j
        
        call quicksort( A )
        
        i = 1
        do while ( A( i ) == 0 )
            i = i + 1
        end do

        nn = 1
        do j = i, n - 1

            if ( A( j ) == A( j + 1 ) ) cycle

            A( nn ) = A( j )
            nn = nn + 1
            
        end do
        if ( A( nn ) /= A( j ) ) A( nn ) = A( j )
        
    end subroutine reorder

    
    recursive subroutine quicksort( A )
        implicit none
        integer, intent(inout) :: A(:)
        integer :: p

        if ( size(A) > 1 ) then
            p = partition( A )
            call quicksort( A(:p-1) )
            call quicksort( A(p:) )
        end if
        
        return
    end subroutine quicksort


    function partition( A ) result( p )
        implicit none
        integer, intent(inout) :: A(:)
        integer :: p

        integer :: pivot, i, j
        
        i = 1
        j = size(A)
        pivot = A(1)

        do

            do while ( A(j) > pivot )
                j = j - 1
            end do

            do while ( A(i) < pivot )
                i = i + 1
            end do

            if ( i < j ) then
                call swap( A(i), A(j) )
            else if ( i == j ) then
                p = i + 1
                exit
            else
                p = i
                exit
            end if
            j = j - 1
            i = i + 1
            
        end do
        
        return
    end function partition

    
    subroutine swap( a, b )
        implicit none
        integer, intent(inout) :: a, b
        integer :: c

        c = a
        a = b
        b = c
        return
    end subroutine swap
    
end module reorder_mod
