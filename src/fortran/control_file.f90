Module control_file
  use num_types
  use vector_matrix_utilities, only: matrix_inverse, volume
  implicit none
CONTAINS

  !!<summary>Read the POSCAR for the lattice and atomic
  !!basis.</summary>
  !!<parameter name="lattice" regular="true">The real space
  !!lattice.</parameter>
  !!<parameter name="atom_base">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="atom_type" regular="true">The type of each atom
  !!in the basis.</parameter>
  SUBROUTINE read_POSCAR(lattice, atom_base, atom_type)
    integer, allocatable, intent(out) :: atom_type(:)
    real(dp), pointer :: atom_base(:,:)
    real(dp), intent(out) :: lattice(3,3)

    character(300) :: line
    integer :: count, j, i, z
    integer, allocatable :: concs(:)

    open(1,file="POSCAR",status="old")

    read(1,'(a300)') line
    read(1,'(a300)') line
    do i=1,3
       read(1,*) lattice(:,i)
    end do
    read(1,'(a300)') line
    call parse_line_for_numbers(line,count,concs)
    allocate(atom_type(sum(concs)), atom_base(3,sum(concs)))
    read(1,'(a300)') line
    do i=1,sum(concs)
       read(1,*) atom_base(:,i)
    end do
    close(1)

    z = 1
    do i=1,count
       if (concs(i) > 0) then
          do j=1,concs(i)
             atom_type(z) = i-1
             z = z+1
          end do
       end if
    end do
  end SUBROUTINE read_POSCAR
  
  !!<summary>Reads the inputs from the 'KPGEN' control file and the
  !!POSCAR, all parameters are output.</summary>
  !!<parameter name="nkpts" regular="true">The number of k-points to
  !!search for.</parameter>
  !!<parameter name="lattice" regular="true">The real space lattice
  !!vectors.</parameter>
  !!<parameter name="atom_base">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="atom_type" regular="true">The type of each atom
  !!in the basis.</parameter>
  !!<parameter name="offset" regular="true">The offset to use on the
  !!k-points.</parameter>
  !!<parameter name="find_offset" regular="true">'True' if the use
  !!didn't provide an offset.</parameter>
  !!<parameter name="eps" regular="true">The floating point
  !!tollerance.</parameter>
  SUBROUTINE get_inputs(nkpts, lattice, atom_type, atom_base, offset, find_offset, eps)
    real(dp), intent(out) :: lattice(3,3), offset(3)
    real(dp), pointer :: atom_base(:,:)
    integer, allocatable, intent(out) :: atom_type(:)
    integer, intent(out) :: nkpts
    logical, intent(out) :: find_offset
    
    ! Input related variables
    character(len=100) :: buffer, label
    integer, parameter :: fh = 15
    integer :: ios, line_count, pos
    
    ! Control file variables
    real(dp) :: pi, lkpd, kpd, r_vecs(3,3), r_vol, eps
    integer :: kppra
    logical :: def_eps, nkpts_set

    ios = 0
    line_count = 0
    pi = 3.1415926535897932385_dp
    find_offset = .True.
    def_eps = .True.
    nkpts_set = .False.
    open(fh, file='KPGEN')

    call read_POSCAR(lattice, atom_base, atom_type)

    call matrix_inverse(transpose(lattice),r_vecs)
    r_vecs = 2*pi*r_vecs
    r_vol = volume(r_vecs(:,1), r_vecs(:,2), r_vecs(:,3))
    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    
    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line_count = line_count + 1
          
          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '=')
          label = buffer(1:pos-1)
          buffer = buffer(pos+1:)
          select case (label)
          case ('SHIFT')
             read(buffer, *, iostat=ios) offset
             find_offset = .False.             
          case ('NKPTS')
             read(buffer, *, iostat=ios) nkpts
             nkpts_set = .True.
          case ('KPDENSITY')
             read(buffer, *, iostat=ios) kpd
             nkpts = NINT(kpd*r_vol)
             nkpts_set = .True.
          case ('KSPACING')
             read(buffer, *, iostat=ios) lkpd
             nkpts = NINT(((1.0_dp/lkpd)**3)*r_vol)
             nkpts_set = .True.
          case ('KPPRA')
             read(buffer, *, iostat=ios) kppra
             nkpts = size(atom_base,2)*kppra
             nkpts_set = .True.
          case ('EPS')
             read(buffer, *, iostat=ios) eps
             def_eps = .False.
             nkpts_set = .True.
          case default
             print *, 'Skipping invalid label at line', line_count, label
          end select
       end if
    end do

    if (def_eps .eqv. .True.) then
       eps = 1E-3
    end if

    if (nkpts_set .eqv. .False.) stop "Number of kpoints not set. Exiting."

  end SUBROUTINE get_inputs

  !!<summary>Splits a string on empty spaces.</summary>
  !!<parameter name="instr" regular="true">The input
  !!string.</parameter>
  !!<parameter name="out1" regular="true">The first item in the
  !!string.</parameter>
  subroutine split_str(instr,out1)
    character(300), intent(inout) :: instr
    character(300), intent(out) :: out1
    integer :: index

    instr = adjustl(TRIM(instr))

    index = SCAN(instr," ")
    out1 = adjustl(TRIM(instr(1:index-1)))
    instr = adjustl(TRIM(instr(index+1:)))
  end subroutine split_str

  !!<summary>Searches a string for numbers.</summary>
  !!<parameter name="line" regular="true">The line to
  !!search.</parameter>
  !!<parameter name="readEntries" regular="true">The number of entries
  !!read in.</parameter>
  !!<parameter name="values" regular="true">The output values
  !!found.</parameter>
  subroutine parse_line_for_numbers(line,readEntries,values)
    character(300), intent(inout):: line
    integer, intent(out)          :: readEntries
    integer, allocatable, intent(out) :: values(:)

    character(300) :: temp
    integer iE, ierr
    logical :: done

    done = .False.
    if (allocated(values)) deallocate(values)
    allocate(values(len(line)))
    values = 0

    readEntries=0
    iE = 0

    line = trim(adjustl(line)) ! Remove preceding blanks
    if (trim(line)=="") then ! make sure we didn't get a blank line
       readEntries = 0
       return
    endif

    do while (.not. done)
       iE = iE+1
       call split_str(line,temp)
       if (temp=="") then
          done = .True.
          exit
       end if
       read(temp,*,iostat=ierr) values(iE)
       if (ierr /=0) then
          readEntries=iE-1
          return
          !stop "ERROR: line parsing for number failed. Maybe format mismatch?"
       endif
       if (line=="") then
          done = .True.
       end if
    enddo

    readEntries=iE
  end subroutine parse_line_for_numbers

end Module control_file
