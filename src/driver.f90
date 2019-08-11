PROGRAM lat_id_driver
  use find_kgrids, only: find_grid
  use kpointgeneration, only: generateIrredKpointList, mapKptsIntoBZ
  use control_file, only: get_inputs
  use vector_matrix_utilities, only: matrix_inverse, determinant, minkowski_reduce_basis
  use numerical_utilities, only: equal
  use num_types

  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3), reduced_R(3,3)
  real(dp) :: Rinv(3,3), point(3), eps, best_offset(3)
  logical :: find_offset
  integer :: nkpts, i, hnf(3,3)
  integer, allocatable :: at(:)
  real(dp), pointer :: IRKps(:,:)
  real(dp), allocatable :: B_vecs(:,:)
  integer, pointer :: weights(:)
  integer :: startTime(8), endTime(8)
  real(dp) :: exeTime, startTimeSec, endTimeSec
  character(10) :: date, time, zone

  call get_inputs(nkpts, lat_vecs, at, B_vecs, offset, find_offset, eps)
  call matrix_inverse(transpose(lat_vecs),r_vecs)

  call date_and_time(date, time, zone, startTime)
  call minkowski_reduce_basis(r_vecs, reduced_R, eps)
  call find_grid(lat_vecs, nkpts, B_vecs, at, offset, find_offset, grid, best_offset, hnf, &
       eps_=eps)
  call date_and_time(date, time, zone, endTime)

  exeTime = 0; startTimeSec = 0; endTimeSec = 0;
  startTimeSec = startTime(5) * 3600 + startTime(6) * 60 + startTime(7) + startTime(8)/1000.0
  if ( endTime(5) < startTime(5) ) then ! add an extra 24 hour
    endTimeSec = (endTime(5) + 24) * 3600 + endTime(6) * 60 + endTime(7) + endTime(8)/1000.0
  else
    endTimeSec = endTime(5) * 3600 + endTime(6) * 60 + endTime(7) + endTime(8)/1000.0
  end if
  exeTime = endTimeSec - startTimeSec

  call generateIrredKpointList(lat_vecs, B_vecs, at, grid, reduced_R, best_offset, &
       IRKps, weights, eps)

  call mapKptsIntoBZ(r_vecs, IRKps, eps)

  open(4,file="KPOINTS")
  write(4,'(A69, I12, A1, I6)')"Dynamic K-point generation: number of irreducible / number reducible: ", size(IRKps,1), '/', sum(weights)
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"

  call matrix_inverse(r_vecs,Rinv)
  do i = 1,size(IRKps,1)
     point = matmul(Rinv,IRKps(i,:))
     write(4,'(F16.14,A1,F16.14,A1,F16.14,A1,I5.1)') point(1), " ",point(2), " ",point(3), " ", weights(i)
  end do
  close(4)

  !print "(4I5)", startTime(5:8)
  !print "(4I5)", endTime(5:8)
  print "(f12.5)", exeTime
  open(10, file="hnf")
  do i = 1, 3
    write(10, '(3I8)') hnf(i,1), hnf(i,2), hnf(i,3)
  end do
  ! write(10, *) best_offset(1), best_offset(2), best_offset(3)
  close(10)

end PROGRAM lat_id_driver
