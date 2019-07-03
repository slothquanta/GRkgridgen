Module grid_utils
  use num_types
  use vector_matrix_utilities, only: matrix_inverse, minkowski_reduce_basis, norm, cross_product, determinant
  use numerical_utilities, only: equal
  use kpointgeneration, only: generateIrredKpointList
  use symmetry, only: get_lattice_pointGroup
  implicit none
  private
  public transform_supercell, grid_selection, grid_metrics, compare_grids

CONTAINS

  !!<summary>Returns the possible symmetry preserving offsets given
  !!the niggli number.</summary>
  !!<parameter name="lat_id" regular="true">The niggli lattice
  !!id.</parameter>
  !!<parameter name="lattice" regular="true">The grid generating
  !!vectors..</parameter>
  !!<parameter name="eps" regular="true">Relative tolerance for
  !!comparisons.</parameter> 
  !!<parameter name="offsets" regular="true">The symmetry preserving
  !!offsets for this case.</parameter>
  SUBROUTINE get_offsets(lat_id, lattice, eps, offsets)
    integer, intent(in) :: lat_id
    real(dp), intent(in) :: lattice(3,3), eps
    real(dp), allocatable, intent(out) :: offsets(:,:)

    real(dp) :: norms(3)

    norms = norm(lattice)
    
    if ((lat_id==1) .or. (lat_id==3) .or. (lat_id==5) .or. (lat_id==32) .or. &
         (lat_id==33) .or. (lat_id==34) .or. (lat_id==35)) then
       allocate(offsets(8,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
       offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
       offsets(4,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
       offsets(5,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
       offsets(6,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
       offsets(7,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
       offsets(8,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
    elseif ((lat_id==12) .or. (lat_id==22) .or. (lat_id==4) .or. (lat_id==2) .or. &
         (lat_id==9) .or. (lat_id==24) .or. (lat_id==6) .or. (lat_id==7) .or. &
         (lat_id==15) .or. (lat_id==18)) then
       if (equal(norms(1), norms(2), eps)) then
          allocate(offsets(2,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
       elseif (equal(norms(1), norms(3), eps)) then
          allocate(offsets(2,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
       elseif (equal(norms(2), norms(3), eps)) then
          allocate(offsets(2,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
       else
          allocate(offsets(4,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
          offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
          offsets(4,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
       end if
    elseif ((lat_id==11) .or. (lat_id==21)) then
       if (equal(norms(1), norms(2), eps)) then
          allocate(offsets(4,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
          offsets(3,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
          offsets(4,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
       elseif (equal(norms(1), norms(3), eps)) then
          allocate(offsets(4,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
          offsets(3,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
          offsets(4,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
       elseif (equal(norms(2), norms(3), eps)) then
          allocate(offsets(4,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
          offsets(3,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
          offsets(4,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
       else
          allocate(offsets(8,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
          offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
          offsets(4,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
          offsets(5,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
          offsets(6,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
          offsets(7,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
          offsets(8,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
       end if
    elseif ((lat_id==16) .or. (lat_id==26)) then
       allocate(offsets(2,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
    elseif ((lat_id==19) .or. (lat_id==8) .or. (lat_id==42)) then
       allocate(offsets(4,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
       offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
       offsets(4,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
    elseif ((lat_id==40) .or. (lat_id==38) .or. (lat_id==36) .or. (lat_id==23) .or. &
         (lat_id==13) .or. (lat_id==10) .or. (lat_id==14) .or. (lat_id==17) .or. &
         (lat_id==20) .or. (lat_id==25) .or. (lat_id==27) .or. (lat_id==28) .or. &
         (lat_id==29) .or. (lat_id==30) .or. (lat_id==37) .or. (lat_id==39) .or. &
         (lat_id==41) .or. (lat_id==43)) then
       if ((norms(1) < norms(2)) .and. (norms(1) < norms(3))) then
          allocate(offsets(8,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
          offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
          offsets(4,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
          offsets(5,:) = (/0.75_dp, 0.25_dp, 0.0_dp/)
          offsets(6,:) = (/0.75_dp, 0.25_dp, 0.5_dp/)
          offsets(7,:) = (/0.25_dp, 0.25_dp, 0.0_dp/)
          offsets(8,:) = (/0.25_dp, 0.25_dp, 0.5_dp/)
       elseif ((norms(2) < norms(1)) .and. (norms(2) < norms(3))) then
          allocate(offsets(8,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
          offsets(3,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
          offsets(4,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
          offsets(5,:) = (/0.25_dp, 0.75_dp, 0.0_dp/)
          offsets(6,:) = (/0.25_dp, 0.75_dp, 0.5_dp/)
          offsets(7,:) = (/0.25_dp, 0.25_dp, 0.0_dp/)
          offsets(8,:) = (/0.25_dp, 0.25_dp, 0.5_dp/)
       elseif ((norms(3) < norms(1)) .and. (norms(3) < norms(2))) then
          allocate(offsets(8,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
          offsets(3,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
          offsets(4,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
          offsets(5,:) = (/0.0_dp, 0.25_dp, 0.75_dp/)
          offsets(6,:) = (/0.5_dp, 0.25_dp, 0.75_dp/)
          offsets(7,:) = (/0.5_dp, 0.25_dp, 0.25_dp/)
          offsets(8,:) = (/0.5_dp, 0.25_dp, 0.25_dp/)
       else 
          allocate(offsets(7,3))
          offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          offsets(2,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
          offsets(3,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
          offsets(4,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
          offsets(5,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
          offsets(6,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
          offsets(7,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
       end if
    elseif ((lat_id==31) .or. (lat_id==44)) then
       allocate(offsets(2,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
    else
       write(*,*) "Failed to find offsets for this case. Please report this occurance. Attempting to use defaults."
       allocate(offsets(2,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)       
    end if
  end SUBROUTINE get_offsets

  !!<summary>Finds the symmetry preserving supercell in the users
  !!basis.</summary>
  !!<parameter name="spHNF" regular="true">One of the symmetry
  !!preserving HNFs in our basis</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="UB" regular="true">The supercell in the users
  !!basis.</parameter>
  SUBROUTINE transform_supercell(spHNF,No,Nu,Co,Cu,O,UB)
    integer, intent(in) :: spHNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(out) :: UB(3,3)
    
    integer :: F(3,3)
    real(dp) :: Oinv(3,3), Noinv(3,3), Cuinv(3,3), L(3,3)

    call matrix_inverse(O,Oinv)
    call matrix_inverse(No,Noinv)
    call matrix_inverse(real(Cu,dp),Cuinv)
    L = matmul(matmul(O,spHNF),Oinv)
    F = nint(matmul(Noinv,matmul(matmul(L,O),Co)))
    UB = matmul(matmul(Nu,F),Cuinv)

  end SUBROUTINE transform_supercell

  !!<summary>Finds the rmin and the number of irreducible k-points for
  !!lattice given an HNF.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="HNF" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="rmin" regular="true">The rmin for this grid.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducbile
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="grid" regular="true">The grid from the HNF.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE grid_metrics(lat_vecs, B_vecs, at, HNF, No, Nu, Co, Cu, O, grid, rmin, &
       n_irr, symm_flag, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), optional, intent(in) :: eps_
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: HNF(3,3), symm_flag
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(out) :: rmin, grid(3,3)
    integer, intent(out) :: n_irr

    real(dp) :: supercell(3,3), shift(3), norms(3), reduced_grid(3,3)
    real(dp) :: grid_inv(3,3), lat_trans(3,3)
    real(dp) :: eps
    real(dp)              :: R(3,3)
    real(dp), pointer     :: rdKlist(:,:)
    integer, pointer      :: weights(:)

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    shift = 0.0_dp

    call transform_supercell(HNF, No, Nu, Co, Cu, O, supercell)

    grid_inv = transpose(supercell)
    call matrix_inverse(grid_inv, grid)
    lat_trans = transpose(lat_vecs)
    call matrix_inverse(lat_trans, R)

    call minkowski_reduce_basis(grid, reduced_grid, eps)
    norms(1) = sqrt(dot_product(reduced_grid(:,1), reduced_grid(:,1)))
    norms(2) = sqrt(dot_product(reduced_grid(:,2), reduced_grid(:,2)))
    norms(3) = sqrt(dot_product(reduced_grid(:,3), reduced_grid(:,3)))
    rmin = min(norms(1), norms(2), norms(3))
    call generateIrredKpointList(lat_vecs, B_vecs, at, grid, R, shift, rdKlist, weights, &
         reps_=eps, symm_=symm_flag)
    n_irr = size(rdKlist,1)

  end SUBROUTINE grid_metrics

  !!<summary>Compares a grid to the metrics of another to see which is
  !!has the better rmin or number of irreducible k-points.</summary>
  !!<summary>Finds the rmin and the number of irreducible k-points for
  !!lattice given an HNF.</summary>
  !!<parameter name="HNF" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="rmin" regular="true">The rmin for these grids.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="pac_limit_" regular="true">Minimal packing
  !!fraction grid must meet.</parameter>
  !!<parameter name="grids" regular="true">A list of grids with the best
  !!r_min.</parameter>
  !!<parameter name="best_hnfs" regular="true">The HNFs for the grids with
  !!this r_min.</parameter>
  !!<parameter name="ngrids" regular="true">The number of grids
  !!currently stored in the list of grids.</parameter>
  SUBROUTINE compare_grids(HNF, No, Nu, Co, Cu, O, grids, rmin, &
       best_HNFs, ngrids, pac_limit_, eps_, U_, space_group_)
    integer, intent(in) :: HNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(inout) :: rmin(2)
    real(dp), allocatable, intent(inout) :: grids(:,:,:,:)
    integer, intent(inout) :: ngrids(2)
    integer, allocatable, intent(inout) :: best_HNFs(:,:,:,:)
    real(dp), optional, intent(in) :: eps_, pac_limit_, U_(3,3)
    integer, optional, intent(in) :: space_group_

    real(dp) :: supercell(3,3), shift(3), reduced_grid(3,3), norms(3)
    real(dp) :: eps
    integer, allocatable :: ralloc_HNFs(:,:,:,:)
    real(dp), allocatable :: ralloc_grids(:,:,:,:)
    real(dp) :: temp_rmin, temp_grid(3,3)
    real(dp) :: pi, pac_limit, pac_frac, U(3,3)
    integer :: space_group, pgsize

    real(dp) :: recip(3,3)
    real(dp), pointer :: pg(:,:,:)
    
    pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164_dp

    if (present(space_group_)) then
       space_group = space_group_
    else
       space_group = 1
    end if
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    if (present(pac_limit_)) then
       pac_limit = pac_limit_
    else
       pac_limit = 0.3_dp
    end if

    if (present(U_)) then
       U = U_
    else
       U = 0.0_dp
    end if

    call matrix_inverse(transpose(U), recip)
    
    shift = 0.0_dp
    if (.not. allocated(grids)) then
       allocate(grids(3,3,10,2))
       grids = 0.0_dp
       ngrids = 0
    end if
    if (.not. allocated(best_HNFs)) then
       allocate(best_HNFs(3,3,10,2))
       best_HNFs = 0
    end if

    call transform_supercell(HNF, No, Nu, Co, Cu, O, supercell)

    call matrix_inverse(transpose(supercell), temp_grid)
    call minkowski_reduce_basis(temp_grid, reduced_grid, eps)
    norms(1) = sqrt(dot_product(reduced_grid(:,1), reduced_grid(:,1)))
    norms(2) = sqrt(dot_product(reduced_grid(:,2), reduced_grid(:,2)))
    norms(3) = sqrt(dot_product(reduced_grid(:,3), reduced_grid(:,3)))
    temp_rmin = min(norms(1), norms(2), norms(3))
    pac_frac = ((4.0_dp*pi*((temp_rmin/2.0_dp)**3))/3.0_dp)/determinant(reduced_grid)

    if (pac_frac >= pac_limit) then
       if (space_group > 1) then
          call get_lattice_pointGroup(reduced_grid, pg)
          pgsize = size(pg,3)
       else
          pgsize = 1
       end if
       if (pgsize >= space_group) then
          if (temp_rmin > (rmin(1)+temp_rmin*eps) .and. &
               (.not. equal(temp_rmin, rmin(1), eps))) then
             rmin(2) = rmin(1)
             ngrids(2) = ngrids(1)
             grids(:,:,:,2) = 0.0_dp
             best_HNFs(:,:,:,2) = 0
             best_HNFs(1:3,1:3,1:ngrids(2),2) = best_HNFs(1:3,1:3,1:ngrids(2),1)
             grids(1:3,1:3,1:ngrids(2),2) = grids(1:3,1:3,1:ngrids(2),1)
       
             rmin(1) = temp_rmin
             ngrids(1) = 1
             grids(:,:,:,1) = 0.0_dp
             best_HNFs(:,:,:,1) = 0
             grids(:,:,ngrids(1),1) = reduced_grid
             best_HNFs(:,:,ngrids(1),1) = HNF

          else if ((temp_rmin > (rmin(2)+temp_rmin*eps)) .and. &
               (.not. equal(temp_rmin, rmin(1), eps)) &
               .and. (.not. equal(temp_rmin, rmin(2), eps))) then
             rmin(2) = temp_rmin
             ngrids(2) = 1
             grids(:,:,:,2) = 0.0_dp
             best_HNFs(:,:,:,2) = 0
             grids(:,:,ngrids(2),2) = reduced_grid
             best_HNFs(:,:,ngrids(2),2) = HNF
          else if (equal(temp_rmin, rmin(1), eps)) then
             if (size(grids,3) == ngrids(1)) then
                allocate(ralloc_HNFs(3,3,ngrids(1),2), ralloc_grids(3,3,ngrids(1),2))
                ralloc_HNFs = best_HNFs
                ralloc_grids = grids
                deallocate(best_HNFs, grids)
                allocate(best_HNFs(3,3,2*ngrids(1),2), grids(3,3,2*ngrids(1),2))
                best_HNFs = 0
                grids = 0.0_dp
                best_HNFs(1:3,1:3,1:ngrids(1),1) = ralloc_HNFs(1:3,1:3,1:ngrids(1),1)
                grids(1:3,1:3,1:ngrids(1),1) = ralloc_grids(1:3,1:3,1:ngrids(1),1)
                best_HNFs(1:3,1:3,1:ngrids(2),2) = ralloc_HNFs(1:3,1:3,1:ngrids(2),2)
                grids(1:3,1:3,1:ngrids(2),2) = ralloc_grids(1:3,1:3,1:ngrids(2),2)
             end if
             ngrids(1) = ngrids(1) + 1
             grids(:,:,ngrids(1),1) = reduced_grid
             best_HNFs(:,:,ngrids(1),1) = HNF
          else if (equal(temp_rmin, rmin(2), eps)) then
             if (size(grids,3) == ngrids(2)) then
                allocate(ralloc_HNFs(3,3,ngrids(2),2), ralloc_grids(3,3,ngrids(2),2))
                ralloc_HNFs = best_HNFs
                ralloc_grids = grids
                deallocate(best_HNFs, grids)
                allocate(best_HNFs(3,3,2*ngrids(2),2), grids(3,3,2*ngrids(2),2))
                best_HNFs = 0
                grids = 0.0_dp
                best_HNFs(1:3,1:3,1:ngrids(1),1) = ralloc_HNFs(1:3,1:3,1:ngrids(1),1)
                grids(1:3,1:3,1:ngrids(1),1) = ralloc_grids(1:3,1:3,1:ngrids(1),1)
                best_HNFs(1:3,1:3,1:ngrids(2),2) = ralloc_HNFs(1:3,1:3,1:ngrids(2),2)
                grids(1:3,1:3,1:ngrids(2),2) = ralloc_grids(1:3,1:3,1:ngrids(2),2)
             end if
             ngrids(2) = ngrids(2) + 1
             grids(:,:,ngrids(2),2) = reduced_grid
             best_HNFs(:,:,ngrids(2),2) = HNF
          end if
       end if
    end if

  end SUBROUTINE compare_grids

  !!<summary>Selects the best grid from the list of grids.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="cand_grids" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="best_grid" regular="true">The best grid given the
  !!criteria.</parameter>
  !!<parameter name="ngrids" regular="true">The number of
  !!grids.</parameter>
  !!<parameter name="cand_HNFs" regular="true">The cand HNFs.</parameter>
  !!<parameter name="offset" regular="true">The offset for the
  !!k-points grid.</parameter>
  !!<parameter name="find_offset" regular="true">'True' if the offset
  !!needs to be determined by the algorithm.</parameter>
  !!<parameter name="lat_id" regular="true">The niggli lattice
  !!id.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that gives
  !!the best reduction for the best grid.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible kpoints
  !!of best_grid</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE grid_selection(lat_vecs, B_vecs, at, cand_grids, cand_HNFs, ngrids, &
       offset, find_offset, lat_id, best_grid, best_HNF, best_offset, n_irr, &
       symm_flag, eps_)
    real(dp), intent(in) :: lat_vecs(3,3), offset(3)
    real(dp), allocatable, intent(in) :: cand_grids(:,:,:,:)
    real(dp), optional, intent(in) :: eps_
    real(dp), intent(out) :: best_grid(3,3), best_offset(3)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in), allocatable :: cand_HNFs(:,:,:,:)
    integer, intent(out) :: best_HNF(3,3), n_irr
    integer, intent(in) :: ngrids(2), symm_flag, lat_id
    logical, intent(in) :: find_offset

    real(dp) :: temp_grid(3,3), grid_offsets(sum(ngrids),3)
    integer :: i, n_irreducible(sum(ngrids)), n_ir_min(1), count, j
    real(dp) :: eps

    real(dp)              :: R(3,3), invLat(3,3), temp_R(3,3)
    real(dp), pointer     :: rdKlist(:,:)
    integer, pointer      :: weights(:)
    integer :: err
    real(dp), allocatable :: offsets(:,:)
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    call matrix_inverse(lat_vecs, invLat)
    temp_R = transpose(invLat)
    call minkowski_reduce_basis(temp_R, R, eps)

    n_irreducible = 0
    count = 0
    do i=1,ngrids(1)
       temp_grid = cand_grids(:,:,i,1)
       if (find_offset .eqv. .False.) then
          allocate(offsets(1,3))
          offsets(1,:) = offset
       else
          call get_offsets(lat_id, temp_grid, eps, offsets)
       end if
       
       do j=1, size(offsets,1)
          call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, &
               offsets(j,:), rdKlist, weights, reps_=eps, symm_=symm_flag, &
               err_=err)
          if (err==0) then
             if (n_irreducible(i) == 0) then
                n_irreducible(i) = size(rdKlist,1)
                grid_offsets(i,:) = offsets(j,:)
             elseif (size(rdKlist,1) < n_irreducible(i)) then
                n_irreducible(i) = size(rdKlist,1)
                grid_offsets(i,:) = offsets(j,:)
             end if
          end if
       end do
       deallocate(offsets)
       count = count + 1
    end do

    do i=1,ngrids(2)
       temp_grid = cand_grids(:,:,i,2)
       if (find_offset .eqv. .False.) then
          allocate(offsets(1,3))
          offsets(1,:) = offset
       else
          call get_offsets(lat_id, temp_grid, eps, offsets)
       end if
       do j=1, size(offsets, 1)
          call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, &
               offsets(j,:), rdKlist, weights, reps_=eps, symm_=symm_flag, &
               err_=err)
          if (err==0) then
             if (n_irreducible(i) == 0) then
                n_irreducible(count+i) = size(rdKlist,1)
                grid_offsets(count+i,:) = offsets(j,:)
             elseif (size(rdKlist,1) < n_irreducible(i)) then
                n_irreducible(count+i) = size(rdKlist,1)
                grid_offsets(count+i,:) = offsets(j,:)
             end if
          end if
       end do
       deallocate(offsets)
    end do
    if (any(n_irreducible > 0)) then
       n_ir_min = minloc(n_irreducible, n_irreducible>0)
       n_irr = n_irreducible(n_ir_min(1))
       best_offset = grid_offsets(n_ir_min(1),:)
       if (n_ir_min(1) <= count) then
          best_grid = cand_grids(:,:, n_ir_min(1),1)
          best_HNF = cand_HNFs(:,:, n_ir_min(1),1)
       else
          best_grid = cand_grids(:,:, n_ir_min(1)-count,2)
          best_HNF = cand_HNFs(:,:, n_ir_min(1)-count,2)
       end if
    else
       best_HNF = 0
       best_grid = 0.0_dp
       best_offset = 0.0_dp
       n_irr = 0
    end if

  end SUBROUTINE grid_selection

end Module grid_utils
