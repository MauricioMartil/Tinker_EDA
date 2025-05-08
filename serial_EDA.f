program analyze_pairs_in_memory
  use atoms , only : n, x, y, z
  use files , only : filename, aname, pname, rname, sname, &
 &                   basefile, suffix, version, openend, abort
  use iounit, only : input, iout, iprm, iseq, ixyz, iarc, icrt, &
 &                   freeunit
  use inform, only : debgu, debgp, debgc, debgr, debgs, debgt, &
 &                   verbose, debug
  use output, only : prterr, die
  ! Add other modules if needed for residue/parameter info (e.g., resdue, kvdws, kchrge)
  implicit none

  integer :: n_total_atoms, num_frames
  integer :: num_residues ! Needs to be set by user code
  integer, allocatable :: residue_atom_indices(:,:) ! (max_atoms_per_res, num_residues)
  integer, allocatable :: atoms_per_residue(:) !(num_residues)
  real*8, allocatable :: full_trajectory(:,:,:) ! (3, n_total_atoms, num_frames)
  character*240 :: arcfile

  ! --- Initialization (Tinker & Program Specific) ---
  call initial
  write (iout,*) 'Tinker Pairwise Analysis In Memory Program'

  ! --- Get Input Archive Filename ---
  call get_arc_filename(arcfile)
  if (arcfile == ' ') call die

  ! --- Read Entire Trajectory into Memory ---
  call read_trajectory_to_memory(arcfile, n_total_atoms, num_frames, full_trajectory)
  if (.not. allocated(full_trajectory)) then
     write (iout,*) 'FATAL ERROR: Failed to read trajectory into memory.'
     call final
     stop
  end if
  write (iout,*) 'Trajectory Read: ', n_total_atoms, ' Atoms, ', num_frames, ' Frames.'

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! TODO: USER CODE SECTION - GET RESIDUE and PARAMETER INFORMATION
  !       - Read Coordinate File (e.g., XYZ) to get structure/residue definitions.
  !       - Read Parameter File (e.g., PRM) if needed for analysis parameters.
  !       - Populate:
  !          - num_residues
  !          - atoms_per_residue(:)
  !          - residue_atom_indices(:,:)
  !       - Load necessary parameters for analyze_nonbonding_pair.

  write(iout, *) "WARNING: Residue/Parameter loading not implemented!"
  ! Example placeholder values (REPLACE WITH REAL DATA LOADING):
  num_residues = 1 ! Set a dummy value
  ! allocate(atoms_per_residue(num_residues))
  ! allocate(residue_atom_indices(1, num_residues)) ! Needs proper sizing
  ! ... load actual data ...

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  ! --- Perform Pairwise Analysis Using In-Memory Trajectory ---
  if (num_residues > 1) then
     call analyze_pairs_from_memory(full_trajectory, n_total_atoms, num_frames, &
 &                                  num_residues, residue_atom_indices, atoms_per_residue)
  else
     write(iout,*) "Skipping pairwise analysis (only 1 residue defined)."
  end if

  ! --- Cleanup ---
  write (iout,*) 'Analysis Complete.'
  if (allocated(full_trajectory)) deallocate (full_trajectory)
  ! TODO: Deallocate residue/parameter arrays if allocated
  if (allocated(atoms_per_residue)) deallocate(atoms_per_residue)
  if (allocated(residue_atom_indices)) deallocate(residue_atom_indices)
  call final

contains

  ! ==================================================================
  subroutine get_arc_filename(arcfile_out)
    character*(*), intent(out) :: arcfile_out
    character*240 :: string
    logical :: exist

    arcfile_out = ' '
    call nextarg (arcfile_out, exist) ! Check command line first
    if (exist) then
       call basefile (arcfile_out)
       call suffix (arcfile_out, 'arc', 'old')
       inquire (file=arcfile_out, exist=exist)
    else
       exist = .false.
    end if

    do while (.not. exist)
       write (iout, 10)
 10    format (/,' Enter Coordinate Archive File Name :  ',$)
       read (input, 20) arcfile_out
 20    format (a240)
       call basefile (arcfile_out)
       call suffix (arcfile_out, 'arc', 'old')
       inquire (file=arcfile_out, exist=exist)
       if (.not. exist) then
           write(iout,*) ' ERROR: Cannot find file: ', trim(arcfile_out)
           call nextarg(string, exist) ! Allow empty input to exit
           if (.not. exist) then
               arcfile_out = ' '
               return
           end if
       end if
    end do
  end subroutine get_arc_filename

  ! ==================================================================
  subroutine read_trajectory_to_memory(arc_filename, n_atoms_out, n_frames_out, trajectory_data)
    use atoms , only : n, x, y, z ! Import n,x,y,z from module
    character(len=*), intent(in) :: arc_filename
    integer, intent(out) :: n_atoms_out, n_frames_out
    real*8, allocatable, intent(out) :: trajectory_data(:,:,:) ! (3, n_atoms, n_frames)

    integer :: iarc_unit, frame_count, stat
    logical :: first, opened
    character*1 :: letter ! For checking format (though assumed ARC)

    n_atoms_out = 0
    n_frames_out = 0
    if (allocated(trajectory_data)) deallocate(trajectory_data)

    ! --- Open File and Get Atom Count ---
    iarc_unit = freeunit()
    open (unit=iarc_unit, file=trim(arc_filename), status='old', iostat=stat)
    if (stat /= 0) then
       write(iout,*) 'ERROR: Could not open ARC file: ', trim(arc_filename)
       call prterr('Read Trajectory')
       return
    end if
    inquire(unit=iarc_unit, opened=opened)
    if (.not. opened) return ! Double check

    ! Check format briefly - assumes ARC format primarily
    read (iarc_unit,'(a1)',iostat=stat) letter
    if (stat /= 0) then
       write(iout,*) 'ERROR: Could not read first character of ARC file.'
       close(iarc_unit)
       return
    end if
    if (letter /= ' ' .and. (letter < '0' .or. letter > '9')) then
        write(iout,*) 'WARNING: File does not appear to be Tinker ARC format.'
        ! Proceed anyway for now, readxyz/readcart should handle errors.
    end if
    rewind(unit=iarc_unit)

    ! Read first frame to get atom count (n is populated in 'atoms' module)
    first = .true.
    call readxyz(iarc_unit) ! Uses ixyz internally, which we set to iarc_unit
    if (abort) then
        write(iout,*) 'ERROR: Could not read header/first frame from ARC file.'
        close(iarc_unit)
        return
    end if
    n_atoms_out = n ! Get atom count from module
    if (n_atoms_out <= 0) then
        write(iout,*) 'ERROR: Invalid atom count read from file (', n_atoms_out, ').'
        close(iarc_unit)
        return
    end if
    rewind(unit=iarc_unit)

    ! --- First Pass: Count Frames ---
    write(iout,*) 'Counting frames in ARC file...'
    frame_count = 0
    first = .true.
    abort = .false.
    do while (.not. abort)
        call readcart(iarc_unit, first) ! readcart reads next frame
        if (.not. abort) then
            frame_count = frame_count + 1
            first = .false.
        end if
    end do
    n_frames_out = frame_count
    write(iout, *) 'Found ', n_frames_out, ' frames.'
    if (n_frames_out <= 0) then
        write(iout,*) 'ERROR: No frames found or read error during counting.'
        close(iarc_unit)
        return
    end if
    rewind(unit=iarc_unit)
    abort = .false. ! Reset abort flag

    ! --- Allocate Memory ---
    allocate(trajectory_data(3, n_atoms_out, n_frames_out), stat=stat)
    if (stat /= 0) then
        write(iout,*) 'ERROR: Failed to allocate memory for trajectory!'
        write(iout,*) 'Required: (3, ', n_atoms_out, ', ', n_frames_out, ') * 8 bytes'
        close(iarc_unit)
        return
    end if

    ! --- Second Pass: Read Frames into Memory ---
    write(iout,*) 'Reading frames into memory...'
    first = .true.
    do frame_count = 1, n_frames_out
       call readcart(iarc_unit, first)
       if (abort) then
          write(iout,*) 'ERROR: Read error on frame ', frame_count, ' during storage pass.'
          deallocate(trajectory_data)
          close(iarc_unit)
          n_frames_out = 0 ! Indicate failure
          n_atoms_out = 0
          return
       end if
       first = .false.
       ! Copy coordinates from atoms module to our array
       trajectory_data(1, 1:n_atoms_out, frame_count) = x(1:n_atoms_out)
       trajectory_data(2, 1:n_atoms_out, frame_count) = y(1:n_atoms_out)
       trajectory_data(3, 1:n_atoms_out, frame_count) = z(1:n_atoms_out)
    end do

    close (unit=iarc_unit)
    write(iout,*) 'Finished reading trajectory to memory.'

  end subroutine read_trajectory_to_memory

  ! ==================================================================
  subroutine analyze_pairs_from_memory(full_traj, n_tot, n_frm, n_res, &
 &                                     res_indices, res_counts)
    real*8, intent(in) :: full_traj(3, n_tot, n_frm)
    integer, intent(in) :: n_tot, n_frm, n_res
    integer, intent(in) :: res_counts(n_res)
    integer, intent(in) :: res_indices(size(res_indices,1), n_res) ! Assumes (max_atoms, n_res)

    integer :: res_i, res_j, frame, i, k, atom_idx_orig, pair_idx
    integer :: nuse_pair, max_atoms_pair
    integer, allocatable :: pair_indices(:) ! Original indices for current pair
    real*8, allocatable :: pair_coords(:,:) ! (3, nuse_pair) - Coords for one frame
    ! TODO: Add allocatable arrays for parameters needed by analyze_nonbonding_pair
    !       e.g., pair_atom_types(:), pair_charges(:), pair_vdw_radii(:), ...

    write(iout,*)
    write(iout,*) 'Starting Pairwise Analysis from Memory...'

    do res_i = 1, n_res - 1
      do res_j = res_i + 1, n_res

         ! --- Determine atoms and allocate for this pair ---
         nuse_pair = res_counts(res_i) + res_counts(res_j)
         allocate(pair_indices(nuse_pair))
         pair_idx = 0
         do k = 1, res_counts(res_i)
             pair_idx = pair_idx + 1
             pair_indices(pair_idx) = res_indices(k, res_i)
         end do
         do k = 1, res_counts(res_j)
             pair_idx = pair_idx + 1
             pair_indices(pair_idx) = res_indices(k, res_j)
         end do

         allocate(pair_coords(3, nuse_pair))
         ! TODO: Allocate pair parameter arrays (size nuse_pair)

         write(iout,'(a,i0,a,i0,a,i0,a)') ' Analyzing Pair: (', res_i, ', ', res_j, ') - ', nuse_pair, ' Atoms'

         ! --- Loop over frames ---
         do frame = 1, n_frm

            ! --- Extract coordinates for the pair for this frame ---
            do k = 1, nuse_pair
               atom_idx_orig = pair_indices(k)
               pair_coords(1:3, k) = full_traj(1:3, atom_idx_orig, frame)
            end do

            ! --- Extract parameters for the pair ---
            ! TODO: User implementation needed to get parameters for the
            !       atoms listed in pair_indices(:) and store them in
            !       the allocated parameter arrays (e.g., pair_atom_types).

            ! --- Call adapted analyze routine (PLACEHOLDER) ---
            !       Pass the extracted coordinates and parameters
            call analyze_nonbonding_pair(pair_coords, nuse_pair) !, pair_atom_types, pair_charges, ...)
            !       Handle/accumulate results from the analysis call for this frame

         end do ! End frame loop

         ! --- Deallocate temporary arrays ---
         deallocate(pair_coords)
         deallocate(pair_indices)
         ! TODO: Deallocate pair parameter arrays

         ! --- Process/Output overall results for pair (i, j) ---
         ! E.g., calculate average energy, print time series, etc.
         write(iout,'(a,i0,a,i0,a)') '  Finished analysis for pair (', res_i, ', ', res_j, ')'


      end do ! End res_j loop
    end do ! End res_i loop

  contains

    subroutine analyze_nonbonding_pair(coords, n_pair_atoms) !, atom_types, charges, vdw_params, ...)
       ! --- Arguments ---
       integer, intent(in) :: n_pair_atoms
       real*8, intent(in) :: coords(3, n_pair_atoms)
       ! TODO: Add intent(in) arguments for necessary parameters:
       !       e.g., integer, intent(in) :: atom_types(n_pair_atoms)
       !       e.g., real*8, intent(in) :: charges(n_pair_atoms)
       !       e.g., real*8, intent(in) :: vdw_radii(n_pair_atoms), vdw_eps(n_pair_atoms)

       ! --- Local Variables ---
       ! Declare variables needed for the energy calculation (e.g., energy terms, loop indices)
       real*8 :: evdw, echarge ! Example energy terms

       ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ! TODO: USER CODE - Analyze Non-Bonding Energy for this Pair/Frame
       !
       ! This is where you insert the core logic from analyze.f for
       ! calculating non-bonding interactions (like VdW, Electrostatics).
       !
       ! Key changes needed from original analyze.f logic:
       !   1. Input: Use the 'coords' array (size n_pair_atoms) passed in.
       !   2. Parameters: Use the parameter arrays passed in (atom_types, charges, etc.).
       !   3. Loops: Loop only over the 'n_pair_atoms' in this subset.
       !   4. Output: Calculate the desired energy components (e.g., evdw, echarge)
       !              and potentially return them or store them in a module variable.
       !
       ! Example sketch (replace with actual calculations):
       ! evdw = 0.0d0
       ! echarge = 0.0d0
       ! do i = 1, n_pair_atoms - 1
       !    do j = i + 1, n_pair_atoms
       !       ! Calculate distance between atom i and j using coords(:,i) and coords(:,j)
       !       ! Calculate VdW energy using parameters for i and j
       !       ! Calculate Charge energy using parameters for i and j
       !       ! Add to evdw, echarge
       !    end do
       ! end do

       ! For now, just print a message:
       if (mod(frame, 100) == 0 .or. frame == n_frm) then ! Print progress occasionally
           write (iout,'(a,i0,a,i0)') '      Placeholder: Analyzing frame ', frame, ' for ', n_pair_atoms, ' atoms.'
       end if
       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    end subroutine analyze_nonbonding_pair

  end subroutine analyze_pairs_from_memory

end program analyze_pairs_in_memory