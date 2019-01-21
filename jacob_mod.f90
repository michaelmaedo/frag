module jacob_mod

    use precision_mod
    use parameters_mod
    
    implicit none
    
contains

    subroutine check_jacob( ndim, numel, coord, eltype, connec )
        integer, intent(in) :: ndim
        integer, intent(in), allocatable :: eltype(:)
        integer, intent(inout) :: numel
        integer, intent(inout), allocatable :: connec(:,:)
        real(DOUBLE), allocatable, intent(in) :: coord(:,:)

        integer :: iel, iel_new, inode
        integer :: nnel, old_numel, nnodside
        integer :: zero_area_iel
        integer :: lnod(20)
        real(double) :: vol

        zero_area_iel = 0
        old_numel = numel
        do iel = 1, old_numel
            
            iel_new = iel - zero_area_iel
            lnod(:) = 0
            vol = ZERO

            if ( ndim == 2 ) then
            
                select case( eltype(iel) )
                case ( 1 )
                    nnel = 3
                case ( 5 )
                    nnel = 4
                case ( 12 )
                    nnel = 3
                case ( 8 )
                    cycle
                end select
                    
                lnod(1:nnel) = connec(1:nnel,iel) 
                lnod(nnel+1) = connec(1,iel)

                vol = area_2d( nnel, coord, connec, lnod )                
                
            else if ( ndim == 3 ) then

                select case( eltype(iel) )
                case ( 1 )
                    nnel = 4
                    nnodside = 3
                case ( 26 )
                    nnel = 6
                    nnodside = 3
                case ( 3 )
                    nnel = 8
                    nnodside = 4
                case ( 8 )
                    cycle
                end select
                
                lnod(1:nnodside) = connec(1:nnodside,iel) 
                lnod(nnodside+1) = connec(1,iel)

                ! vol = area_2d( nnodside, coord, connec, lnod )                
                ! vol = vol_3d( nnel, coord, connec, lnod )                
                
            end if
                
            ! write(*,*) vol
            connec(1:nnel,iel_new) = lnod(1:nnel)
            if ( vol < -1.0e-30 ) then
                do inode = 1, nnel
                    connec(inode,iel_new) = lnod(nnel-(inode-1))
                end do
            else if ( vol < 1.0e-20 ) then
                zero_area_iel = zero_area_iel + 1
                numel = numel - 1
            end if
            
        end do
        
    end subroutine check_jacob

    
    function area_2d( nnel, coord, connec, lnod ) result ( area )
        integer, intent(in) :: nnel
        integer, intent(in) :: lnod(20)
        integer, intent(in), allocatable :: connec(:,:)
        real(double), intent(in), allocatable :: coord(:,:)

        integer :: inode
        integer :: ii, jj
        real(DOUBLE) :: x1, y1, x2, y2, area

        
        area = 0.0D+00
        do inode = 1, nnel

            ii = lnod(inode)
            jj = lnod(inode+1)
            
            x1 = coord(1,ii)
            y1 = coord(2,ii)

            x2 = coord(1,jj)
            y2 = coord(2,jj)

            area = area + (x1*y2 - y1*x2)

        end do
        
    end function area_2d
    
end module jacob_mod
