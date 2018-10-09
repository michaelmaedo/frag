module mesh_frag_mod

contains

    subroutine find_neighbors( nnode, nelem, mnnel, nelfrag,  &
        nnodeface, nfacel, connec, elset, elfrag, iordface, face )

        implicit none      

        integer, intent(in) :: nnode, nelem, mnnel, nelfrag, nnodeface, nfacel

        integer, intent(in) :: connec( mnnel, nelem ), &
                               elset( nnode ), &
                               elfrag( nnode ), &
                               iordface( nnodeface, nfacel )

        integer, intent(inout) :: face( nfacel, nelfrag )

        integer :: ii, iel, iface, ino_face, el_candidate, ino, node, nnel
        integer :: lnodface( nnodeface )
        
        nnel = 3
        do ii = 1, nelfrag
            iel = elfrag( ii )

            do iface = 1, nfacel

                do ino_face = 1, nnodeface
                    lnodface( ino_face ) = &
                        connec( iordface( ino_face, iface ), iel )
                end do

                do el_candidate = 1, nelem
                    if ( iel /= el_candidate ) then

                        node = 0
                        do ino = 1, nnel
                            do ino_face = 1, nnodeface
                                if ( connec( ino, el_candidate ) ==  &
                                    lnodface( ino_face ) ) then                             
                                    node = node + 1
                                end if
                            end do
                        end do

                        if ( node == nnodeface ) then

                            face( iface, ii ) = el_candidate                            
                            if ( elset(iel) /= elset(el_candidate) ) &
                                face( iface, ii ) = -el_candidate
                            exit

                        end if

                    end if
                end do
            end do
        end do

    end subroutine find_neighbors

    
    subroutine reduce_size_of_elements( &
        nelem, ndim, mnnel, nnodeface, nfacel, nelfrag, new_nnode, &
        thick, connec, iordface, elfrag, face, new_coord )
        implicit none

        real(kind(1.0D+00)), intent(in) :: thick
        integer,intent(in) ::  nelem, ndim, mnnel, nnodeface, &
                               nfacel, nelfrag, new_nnode 

        integer, intent(in) :: iordface( nnodeface, nfacel ), &
                               elfrag( nelfrag ), &
                               face( nfacel, nelfrag )
        integer, intent(inout) :: connec( mnnel, nelem )

        real(kind(1.0D+00)), intent(inout) :: new_coord( 3, new_nnode )

        real(kind(1.0D+00)) :: h, r, theta1, theta2
        integer :: ii, iel, iface, ino_face, node

        real(kind(1.0D+00)) :: x(mnnel), y(mnnel)
        real(kind(1.0D+00)) :: u(ndim), v(ndim), w(ndim)
        integer :: lnodface( nnodeface )

        w(:) = 0.0D+00
        do ii = 1, nelfrag
            iel = elfrag( ii )

            do iface = 1, nfacel
                h = thick/2
                if ( face(iface,ii) < 0 ) h = thick

                do ino_face = 1, nnodeface
                    lnodface( ino_face ) = &
                        connec( iordface( ino_face, iface ), iel )
                end do

                if (iface == 3 ) then
                    node = connec( iordface( 2, 1 ), iel )
                else
                    node = connec( iordface( 2, iface+1 ), iel )
                end if

                x(1) = new_coord( 1, lnodface(1) )
                y(1) = new_coord( 2, lnodface(1) )

                x(2) = new_coord( 1, lnodface(2) )
                y(2) = new_coord( 2, lnodface(2) )

                x(3) = new_coord( 1, node )
                y(3) = new_coord( 2, node )

                u(1) = x(2) - x(1)
                u(2) = y(2) - y(1)
                u(:) = u(:)/norm(u)

                v(1) = x(3) - x(2)
                v(2) = y(3) - y(2)
                v(:) = v(:)/norm(v)

                w(1) = x(3) - x(1)
                w(2) = y(3) - y(1)
                w(:) = w(:)/norm(w)

                theta1 = acos ( dot(v,u) )                
                theta2 = acos ( dot(w,u) )

                r = h/sin(theta2)
                new_coord( 1, lnodface(1) ) = x(1) + w(1)*r
                new_coord( 2, lnodface(1) ) = y(1) + w(2)*r

                r = h/sin(theta1)
                new_coord( 1, lnodface(2) ) = x(2) + v(1)*r
                new_coord( 2, lnodface(2) ) = y(2) + v(2)*r
            end do
        end do

    end subroutine reduce_size_of_elements


    subroutine insert_interface_elements(                          &
        nelem, ndim, mnnel, nmat, nelfrag, nfacel, nnodeface,      &
        new_nnode, new_nmat, new_nelem, new_coord, connec, elfrag, &
        iordface, face, new_connec, new_elset )

        integer, intent(in)    :: nelem, ndim, mnnel, nmat, nelfrag, &
                                  nfacel, nnodeface, new_nnode
        integer, intent(inout) :: new_nmat, new_nelem

        real(kind(1.0D+00)), intent(in) :: new_coord( 3, new_nnode )
        integer            , intent(in) :: connec( mnnel, nelem ), &
                                           iordface( nnodeface, nfacel ), &
                                           elfrag( nelfrag ), &
                                           face( nfacel, nelfrag )
        
        integer, intent(inout) :: new_connec( mnnel, new_nelem ), &
                                  new_elset( new_nelem )

        real(kind(1.0D+00)) ::  x1, y1, dist, dist_candidate
        integer :: ii, iel, iface, el, ino_face, ino

        real(kind(1.0D+00)) :: x( mnnel ), y( mnnel )
        integer :: lnodface( nnodeface ), neighbor_nodes( nnodeface ) 

        new_nelem = nelem
        nnel = 3
        do ii = 1, nelfrag
            iel = elfrag( ii )
            
            do iface = 1, nfacel
                el = face( iface, ii )
                
                if ( abs(el) > iel .or. el < 0 ) then

                    do ino_face = 1, nnodeface
                        lnodface( ino_face ) = &
                            connec( iordface( ino_face, iface ), iel )
                        
                        x1 = new_coord( 1, lnodface( ino_face ) ) 
                        y1 = new_coord( 2, lnodface( ino_face ) )

                        dist = 1.0D+9;
                        do ino = 1, nnel
                            x(ino) = new_coord ( 1, connec( ino, abs(el) ) )
                            y(ino) = new_coord ( 2, connec( ino, abs(el) ) )
                            
                            !distance between 2 points
                            dist_candidate = sqrt( (x1 - x(ino))**2 + (y1 - y(ino))**2 )

                            if ( dist > dist_candidate ) then
                                dist = dist_candidate
                                neighbor_nodes(ino_face) = connec( ino, abs(el) )
                            end if
                        end do

                    end do

                    new_nelem = new_nelem + 1
                    new_connec( 1 : nnodeface, new_nelem ) = lnodface( 1 : nnodeface )
!                    new_connec( nnodeface + 1, new_nelem ) = neighbor_nodes( 1 )
                    new_connec( nnodeface + 1, new_nelem ) = neighbor_nodes( 2 )

                    new_elset( new_nelem ) = nmat + 1
                    if ( el < 0 ) new_elset( new_nelem ) = nmat + 2

                    new_nelem = new_nelem + 1
                    new_connec( 1 : nnodeface, new_nelem ) = neighbor_nodes( 1 : nnodeface )
!                    new_connec( nnodeface + 1, new_nelem ) = lnodface( 2 )
                    new_connec( nnodeface + 1, new_nelem ) = lnodface( 1 )

                    new_elset( new_nelem ) = nmat + 1
                    if ( el < 0 ) new_elset( new_nelem ) = nmat + 2

                end if

            end do
        end do

        return
    end subroutine insert_interface_elements


    function dot(v1, v2) result (res)
        implicit none
        real(kind(1.0d+00)) :: v1(:), v2(:)
        real(kind(1.0d+00)) :: res
        integer :: idim, ndim

        ndim = size(v1)

        res = 0.0d+00
        do idim = 1, ndim
            res = res + v1(idim)*v2(idim)
        end do

    end function dot


    function norm(vec) result (res)
        implicit none
        real(kind(1.0d+00)) :: vec(:)
        real(kind(1.0d+00)) :: res
        integer :: idim, ndim

        ndim = size(vec)

        res = 0.0d+00
        do idim = 1, ndim
            res = res + ( vec(idim) * vec(idim) )
        end do
        res = sqrt( res )

    end function norm

end module mesh_frag_mod
