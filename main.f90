program main
    use mesh_frag_mod
    use reorder_mod
    use jacob_mod
    
    implicit none
    
    integer :: nnode, nelem, ndime, nnel, nmat, &
               new_nnode, new_nelem, new_nmat, &
               ino, iel, jno, new_ino, nbarnods, &
               ii, jj, kk
               
    integer, allocatable :: connec(:,:), &
                            elset(:), &
                            eltype(:), &
                            bcond(:,:), &
                            new_connec(:,:), &
                            new_elset(:), &
                            new_eltype(:), &
                            new_bcond(:,:), &
                            lnods(:), &
                            new_nods_position(:)

    real(kind(1.0D+00)) :: thick
    integer :: node, fragmat, nelfrag, nfacel, nnodeface, no_bar
    character(len = 500) str
    
    real(kind(1.0D+00)), allocatable :: coord(:,:), &
                                        new_coord(:,:)
    integer, allocatable :: elfrag(:), &
                            bar_bc(:), &
                            bar_nodes(:), &
                            iordface(:,:), &
                            lnodface(:), &
                            face(:,:)
    
    integer, parameter :: inpf = 10, outf = 20, grif = 30
                            
    call greeting
    
    open(unit = inpf, file = 'inp.dat', status = 'unknown')
    open(unit = outf, file = 'para.vtu', status = 'unknown')
    open(unit = grif, file = 'out.dat', status = 'unknown')

    write(*,'(A)') " Read data from the input file"
    write(*,'(A)') " -----------------------------"
    
    write(*,*) "Reading control parameters "
    read(inpf,*) nnode, nelem, ndime
    read(inpf,*) nmat, fragmat, thick
    read(inpf,*) no_bar
    
    allocate(coord(3,nnode))
    allocate(elset(nelem))
    allocate(connec(4,nelem))
    allocate(bcond(2,nnode))
    allocate(eltype(nelem))
    allocate(lnods(10))
    
    
!=====================================================================
!   Read from input file
!=====================================================================

    write(*,*) "Reading nodal coordinates "
    coord(:,:) = 0.0D+00
    do ino = 1, nnode
        read(inpf,*) ii, (coord(jj,ino), jj = 1, ndime), &
            (bcond(kk,ino), kk = 1, 2)
    end do
    
    write(*,*) "Reading connectivity list "
    nnel = 3
    do iel = 1, nelem
        read(inpf,'(6i10)') ii, elset(iel), eltype(iel), &
            (connec(jj,iel), jj = 1, nnel)
    end do
    close(inpf)
    write(*,'(A)') ""

!=====================================================================
!   Find bar nodes
!=====================================================================
    if (ndime == 2) then
        nnodeface = 2
        nfacel = 3 ! for triangles
    else if (ndime == 3) then
        nnodeface = 3
        nfacel = 4 ! for tetrahedron
    else
        write(*,*) 'Incorrect parameter: ndime'
        write(*,*) '...'
        stop
    end if
    
    allocate( elfrag(nelem) ); elfrag(:) = 0
    allocate( bar_bc(nnode) ); bar_bc(:) = 0

    nelfrag = 0
    do iel = 1, nelem
        if (elset(iel) == fragmat) then
            nelfrag = nelfrag + 1
            elfrag(nelfrag) = iel
        else
            do ino = 1, nnel
                bar_bc( connec(ino,iel) ) = -1
            end do
        end if
    end do

!=====================================================================
!   Find element's neighbor
!=====================================================================
    allocate( face( nfacel, nelfrag ) ); face(:,:) = 0
    allocate( lnodface(nnodeface) ); lnodface(:) = 0
    allocate( iordface(nnodeface,nfacel) ); iordface(:,:) = 0

    iordface(1,1) = 1;  iordface(2,1) = 2
    iordface(1,2) = 2;  iordface(2,2) = 3
    iordface(1,3) = 3;  iordface(2,3) = 1
    
    write(*,'(A)') " Fragment Finite Element Mesh"
    write(*,'(A)') " ----------------------------"

    write(*,*) "Finding the neighbors of the elements "
    call find_neighbors( nnode, nelem, 4, nelfrag,  &
        nnodeface, nfacel, connec, elset, elfrag, iordface, face )
    
!=====================================================================
!   Reduce size of the elements
!=====================================================================
    new_nnode = nnode + nnel*nelfrag
    allocate(new_coord(3,new_nnode)); new_coord(:,:) = 0.0D+00
    new_coord(1:3,1:nnode) = coord(1:3,1:nnode)

    allocate(bar_nodes(new_nnode)); bar_nodes(:) = 0
    allocate(new_bcond(2,new_nnode)); new_bcond(:,:) = 0
    new_bcond(1:2,1:nnode) = bcond(1:2,1:nnode)

    node = nnode
    do ii = 1, nelfrag
        iel = elfrag( ii )

        lnods(1:nnel) = connec(1:nnel,iel)
            
        do ino = 1, nnel

            node = node + 1
                
            connec(ino,iel) = node
                
            new_coord( 1, node ) = coord( 1, lnods(ino) )
            new_coord( 2, node ) = coord( 2, lnods(ino) )

            bar_nodes( node ) = lnods( ino )
            new_bcond( 1, node ) = bcond( 1, lnods(ino) )
        end do
        
    end do
    new_nnode = node
    
    write(*,*) "Reducing the size of the elements"
    call reduce_size_of_elements( &
        nelem, ndime, 4, nnodeface, nfacel, nelfrag, new_nnode, &
        thick, connec, iordface, elfrag, face, new_coord )
    
!=====================================================================
!   Insert elements with high aspect ratio
!=====================================================================
    
    new_nelem = nelem + 6*nelfrag + (new_nnode - nnode + 1)
    allocate( new_connec( 4, new_nelem ) );   new_connec(:,:) = 0
    new_connec(1:4,1:nelem) = connec(1:4,1:nelem)
    
    allocate(new_elset(new_nelem)); new_elset(:) = 0
    new_elset(1:nelem) = elset(1:nelem)

    allocate(new_eltype(new_nelem)); new_eltype(:) = 1
    new_eltype(1:nelem) = eltype(1:nelem)
    
    write(*,*) "Inserting finite elements with high aspect ratio"
    new_nmat = nmat
    call insert_interface_elements(                                &
        nelem, ndime, 4, nmat, nelfrag, nfacel, nnodeface,         &
        new_nnode, new_nmat, new_nelem, new_coord, connec, elfrag, &
        iordface, face, new_connec, new_elset )
    

!=====================================================================
!   Insert Bar elements
!=====================================================================
    if (no_bar == 1) then

        write(*,*) "Inserting 1D-Bar elements"
        do ino = nnode+1, new_nnode
            new_nelem = new_nelem + 1
            new_connec(1,new_nelem) = bar_nodes( ino )
            new_connec(2,new_nelem) = ino
            new_connec(3,new_nelem) = bar_nodes( ino )
            new_eltype(new_nelem) = 8
        end do

    else

        allocate( new_nods_position(new_nnode) )
        new_nods_position(:) = 0

!       Remove old nodes that should be used by 1D-Bar elements
        call reorder( new_nnode, nbarnods, bar_nodes )
        
!       Remove from 'bar_nodes' nodes that are shared with neighbor material
!       (i.e., preserve nodes from neighbor material)        
        node = 1
        do jno = 1, nbarnods
            do ino = 1, nnode
                if ( bar_bc(ino) == -1 .and.  &
                    ino == bar_nodes( jno ) ) then
                    bar_nodes( jno ) = 0
                    exit                    
                end if
            end do
        end do

        jno = 0
        do ino = 1, nbarnods
            if ( bar_nodes( ino ) /= 0 ) then
                jno = jno + 1
                bar_nodes( jno ) = bar_nodes( ino )
            end if
        end do
        nbarnods = jno
!
!       Reorder the node positions 
        new_ino = 0
        node = 1
        do ino = 1, new_nnode

            if ( ino == bar_nodes( node ) ) then
                node = node + 1
                if ( node > nbarnods ) node = nbarnods
            else
                new_ino = new_ino + 1
                new_coord(:,new_ino) = new_coord(:,ino)
                new_bcond(:,new_ino) = new_bcond(:,ino)
            end if
            new_nods_position(ino) = new_ino
        end do
        new_nnode = new_nnode - jno
!
!       Reorder the position of the nodes in the elements
        do iel = 1, new_nelem
            do ino = 1, nnel
                new_connec(ino,iel) = &
                    new_nods_position( new_connec( ino, iel ) )
            end do
        end do
        
    end if

!=====================================================================
!   Mechanical Boundary Conditions of the Bar Elements
!=====================================================================
    if (no_bar == 1) then
        write(*,*) "Correcting boundary of the bar element"
        do ino = 1, nnode
            if (bar_bc(ino) /= -1) then
                bcond(1,ino) = 5
            end if
        end do
    end if
    write(*,'(A)') ""

!=====================================================================
!   Write to grid file
!=====================================================================
    call check_jacob( ndime, new_nelem, new_coord, &
        new_eltype, new_connec)
    
!=====================================================================
!   Write to grid file
!=====================================================================
    write(*,'(A)') " Write fragmented mesh to the output file"
    write(*,'(A)') " ----------------------------------------"

100 format("(i10,", i3,"(' ',e20.13),2i10)")
    write(str,100) ndime 
    write(grif,*) new_nnode, new_nelem
    do ino = 1, new_nnode
        write(grif,str) ino, &
            (new_coord(jj,ino), jj = 1, 2), &
            new_bcond(1,ino), new_bcond(2,ino)
    end do

200 format("(3i10,", i3,"i10)")
    write(str,200) nnel 
    do iel = 1, new_nelem
        write(grif,str) iel, new_elset(iel), 1, &
            (new_connec(ino,iel), ino = 1, nnel)
    end do
    
    
!=====================================================================
!   Write to output file
!=====================================================================
   
    write(outf,500) new_nnode, new_nelem
500 format('<?xml version="1.0"?>'/&
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'/&
        '<UnstructuredGrid>'/&
        '<Piece NumberOfPoints="', i10,'" NumberOfCells="', i10,'">'/&
        '<Points>'/&
        '<DataArray type="Float64"  NumberOfComponents="3" format="ascii">')
    
    write(*,*) "Writing nodal coordinates "    
    do ino = 1, new_nnode
        write(outf,'(3(" ",e20.13))') &
            (new_coord(jj,ino), jj = 1, 3)
    end do
    write(outf,510)
510 format('</DataArray>',/,'</Points>',/,'<Cells>')

    
    write(*,*) "Writing connectivity list"       
520 format('<DataArray type="Int64"       Name="connectivity" format="ascii">')
    write(outf,520)
    do iel = 1, new_nelem
        do ino = 1, nnel
            lnods(ino) = new_connec(ino,iel) - 1
        end do
        write(outf,"(3i10)") (lnods(jj), jj = 1, 3)
    end do
    write(outf,530)
530 format('</DataArray>')

540 format('<DataArray type="Int64"            Name="offsets" format="ascii">')
    write(outf,540)
    kk = 0
    do iel = 1, new_nelem
        kk = kk + nnel
        write(outf,*) kk
    end do
    write(outf,530)
    
550 format('<DataArray type="Int64"            Name="types" format="ascii">')
    write(outf,550)
    do iel = 1, new_nelem
        write(outf,'(i5)') 5
    end do
    write(outf,530)

    
    write(outf,560)
560 format('</Cells>',/&
        '</Piece>',/&
        '</UnstructuredGrid>',/&
        '</VTKFile>')
    write(*,'(A)') ""

    close(grif)
    close(outf)
    
    call free(no_bar, connec, elset, eltype, bcond, lnods, &
        new_connec, new_elset, new_eltype, new_bcond, &
        new_nods_position, coord, new_coord, elfrag, bar_bc, &
        bar_nodes, iordface, lnodface, face )
    
contains

    subroutine greeting()

        write(*,'(A)')  ""
        write(*,'(A)')  "     F R A G"
        write(*,'(A)')  "     version 0"
        write(*,'(A)')  ""
        write(*,'(A)')  "     Insert finite elements with high aspect ratio,"
        write(*,'(A)')  "     also  called  interface  elements, between the"
        write(*,'(A)')  "     standard (bulk) elements of the original mesh"
        write(*,'(A)')  ""
        write(*,'(A)')  ""
        
    end subroutine greeting


    subroutine free(no_bar, connec, elset, eltype, bcond, lnods, &
        new_connec, new_elset, new_eltype, new_bcond, &
        new_nods_position, coord, new_coord, elfrag, bar_bc, &
        bar_nodes, iordface, lnodface, face )

        integer :: no_bar
        
        integer, allocatable :: connec(:,:), &
            elset(:), &
            eltype(:), &
            bcond(:,:), &
            new_connec(:,:), &
            new_elset(:), &
            new_eltype(:), &
            new_bcond(:,:), &
            lnods(:), &
            new_nods_position(:), &
            elfrag(:), &
            bar_bc(:), &
            bar_nodes(:), &
            iordface(:,:), &
            lnodface(:), &
            face(:,:)

        real(kind(1.0D+00)), allocatable :: &
            coord(:,:), &
            new_coord(:,:)

        deallocate( connec, elset, eltype, bcond, lnods, &
        new_connec, new_elset, new_eltype, new_bcond, &
        coord, new_coord, elfrag, bar_bc, &
        bar_nodes, iordface, lnodface, face )
        
        if ( no_bar == 0 ) deallocate( new_nods_position ) 
        
    end subroutine free
    
end program main
