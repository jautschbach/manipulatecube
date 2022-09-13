program cubcub
  
  ! manipulates Gaussian cube files
  ! (c) copyright Jochen Autschbach (2020-22)
  ! all rights reserved
  
  
  use types
  
  implicit none
  
  character*(LCHARS) :: cubfile1, cubfile2, arg1, arg2, arg3, arg4, &
    res1, res2
  
  character(len=3) :: math
  
  integer(KINT), parameter :: inp1=7, inp2=8, out=6, &
    cub1=9, cub2=10, err=0
  
  integer(KINT) :: nnuc, npts(3),  nnuc2, npts2(3), i, j, k, &
    ngrid, info
  
  real(KREAL) :: vectr(3,3), startp(3), facmul, alpha, &
    rint, dV, thresh, vectr2(3,3), startp2(3)
  
  real(KREAL), parameter :: one=1.0d0,sq2=sqrt(2.0d0),sq2m1=one/sq2, &
    eight=8.0d0,twopi=eight*atan(one),radians=twopi/360.0d0,zero=0.0d0,&
    p005=0.05d0, small=1d-3
  
  real(KREAL), allocatable :: xyznuc(:,:), cubvalue1(:,:,:), &
    cubvalue2(:,:,:), qtch(:)
  real(KREAL), allocatable :: output(:,:,:), points(:)
  integer(KINT), allocatable :: nuctyp(:)
  
  integer :: narg, ios
  
!!$  integer(KINT), parameter :: nelmnt=118
!!$  character(3) :: symbol(0:nelmnt)
  
  ! ==========================================================================
  !Read Input from user 
  narg = iargc()
  
  if ( narg .lt. 2 .or. narg .gt. 4 ) then
    stop &
      'Usage: ./manipulatecube cubfile1 add/sub/mix/rot/mul/fac cubfile2|factor [angle]' 
  end if
  
  call getarg(1,arg1)
  call getarg(2,arg2)
  call getarg(3,arg3)
  if (narg.eq.4) call getarg(4,arg4)
  
  ! check the option for what we're supposed to do with the cubes
  
  ios = 0
  read(arg2,*,iostat=ios) math
  if (ios /= 0) stop 'error reading arg2'
  !write (out,*) math
  if ( (math .eq. 'add') .or. (math .eq. 'sub') .or. &
    (math .eq. 'mix') .or. (math .eq. 'mul') .or. &
    (math .eq. 'fac') .or. (math .eq. 'rot') .or. &
    (math .eq. 'dmp')) then
    write (out,*) 'will perform operation: ',math
  else
    stop 'arg2: unrecognized option'
  end if
  
  if (math.eq.'rot' .and. narg.eq.4) then
    ios = 0
    read(arg4,*,iostat=ios) alpha
    if (ios /= 0) stop 'error reading rotation angle in argument 4'
  elseif (math.eq.'rot' .and. narg.ne.4) then
    stop 'rot procedure requires rotation angle in argument 4'
  end if
  
  
  ! the first cube file in the list is special in the sense that
  ! we'll read the grid specs from that file
  
  ios = 0
  read(arg1,*,iostat=ios) cubfile1
  if (ios /= 0) stop 'error reading arg1'
  
  call open_cube_read_header(inp1,cubfile1,nnuc,startp,npts,vectr)
  
  allocate (cubvalue1(npts(1),npts(2),npts(3)))
  
  cubvalue1 = 0
  
  ios=0
  do i = 1, npts(1)
    do j = 1, npts(2)
      read (inp1, *, iostat=ios) (cubvalue1(i,j,k), k=1,npts(3))
      if (ios /= 0) then
        write (err,*) 'i, j =', i,j
        write (err,*) 'error reading data from file '//trim(cubfile1)
      end if
    end do
  end do
  
  close (inp1)
  
  
  ! check the options and see if we need to read a second cube file.
  ! if so, red header 
  
  if (math /= 'fac' .and. math /= 'dmp') then
    
    ios = 0
    read(arg3,*,iostat=ios) cubfile2
    if (ios /= 0) stop 'error reading arg3'
    
    call open_cube_read_header(inp2,cubfile2,nnuc2,startp2,npts2,vectr2)
    
    ! check that the cube files have the same specs as far as
    ! number of nuclei and number of data points are concerned
    
    if (nnuc2.ne.nnuc) then
      write (err,*) 'incompatible number of atoms in the cubes'
      stop 'error termination'
    end if
    do i =1,3
      if (npts2(i).ne.npts(i)) then
        write(err,*) npts(i), npts2(i)
        write (err,*) 'incompatible number of grid points in the cubes'
        stop 'error termination'
      end if
    end do
    
    allocate (cubvalue2(npts(1),npts(2),npts(3)))
    cubvalue2 = 0
    
    do i = 1, npts(1)
      do j = 1, npts(2)
        read (inp2, *, iostat=ios) (cubvalue2(i,j,k), k=1,npts(3))
        if (ios /= 0) then
          write (err,*) 'i, j =', i,j
          write (err,*) 'error reading data from file '//trim(cubfile2)
        end if
      end do
    end do
    
    close (inp2)
    
  else if (math .eq. 'fac') then
    
    ! the 3rd argument is a numerical value by which to multiply the cube
    
    ios = 0
    read(arg3,*,iostat=ios) facmul
    if (ios /= 0) stop 'error reading arg3'
    
  end if
  
  ! create output cube and write header
  
  if (math .eq. 'mix' .or. math .eq. 'rot') then
    res1 = 'results1.cube'
    res2 = 'results2.cube'
  else if (math .eq. 'fac') then
    res1 = 'results.cube'
    res2 = ''
  else
    res1 = 'data.txt'
    res2 = ''
  end if
  
  
  ios = 0
  open (cub1, file=trim(res1), status='unknown', iostat=ios)
  if (ios /= 0) then
    write (err,*) 'error creating or openinng file '//trim(res1)
    stop 'error termination'
  end if
  
  if (math /= 'dmp') call write_cube_header (cub1)
  
  ! create second output cube if needed
  
  if (math .eq. 'mix' .or. math .eq. 'rot') then
    
    ios = 0
    open (cub2, file=trim(res2), status='unknown', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'error creating or openinng file '//trim(res2)
      stop 'error termination'
    end if
    
    call write_cube_header (cub2)  
    
  end if
  
  
  ! now do the requested operation
  
  allocate (output(npts(1),npts(2),npts(3)))
  
  output=0
  
  select case (math)
  case ('add')
    !Deal with Addition
    write (out,*) 'Adding cube files'
    output = cubvalue1 + cubvalue2
    
    call write_data(cub1,err,res1,npts,output)
    
  case ('sub')
    !Deal with Subtraction
    write (out,*) 'Subtracting cube files'
    output = cubvalue1 - cubvalue2
    
    call write_data(cub1,err,res1,npts,output)
    
  case ('mul')
    !Deal with multiplication
    write (out,*) 'Multiplying cube files'
    output = cubvalue1 * cubvalue2
    
    call write_data(cub1,err,res1,npts,output)
    
  case ('mix')
    !Deal with mixing
    write (out,*) 'Mixing cube files'
    
    output = cubvalue1 + cubvalue2
    output = sq2m1 * output
    
    call write_data(cub1,err,res1,npts,output)
    
    output = cubvalue1 - cubvalue2
    output = sq2m1 * output
    
    call write_data(cub2,err,res2,npts,output)
    
  case ('rot')
    !Deal with rotation
    write (out,'(1x,a,f8.2,a,f8.2,a)') 'Rotating cube files by angle ',alpha, &
      ' degrees or ', alpha*radians,' radians'
    
    !write(out,*) twopi, radians
    alpha = alpha * radians
    !write(out,*) cos(alpha), sin(alpha),sq2m1
    
    output = cos(alpha)*cubvalue1 + sin(alpha)*cubvalue2
    
    call write_data(cub1,err,res1,npts,output)
    
    output = -one*sin(alpha)*cubvalue1 +cos(alpha)*cubvalue2
    
    call write_data(cub2,err,res2,npts,output)
    
  case ('fac')
    !Deal with multiplication by a factor
    write (out,*) 'Multiplying cube file with a factor of ',facmul
    
    output = cubvalue1 * facmul
    
    call write_data(cub1,err,res1,npts,output)
    
    ! let's also take the integral of the cube, since it's easy to do here
    call volume_element(vectr(:,1), vectr(:,2), vectr(:,3), dV)
    rint = sum(cubvalue1) * dV
    write (out,*) 'Volume Integral of the cube: ',rint
    
    ! ... and another hack: if the cube integral is not zero, we assume it is
    ! a density and determine what isosurface value contains what fraction of
    ! density
    
    if (rint>small) then
      
      ngrid = npts(1) * npts(2) * npts(3)
      allocate(points(ngrid))
      
      write(out,*) 'iso, threshold, cumulated integral'
      
      points = reshape(cubvalue1,[ngrid])

      info = 0
      call dlasrt('I', ngrid, points, info)
      if (info.ne.0) write (out,*) &
        'WARNING: DLASRT Failed with INFO =',info
      
      ! make sure array 'points' is ordered lowest to highest
      do i=2,ngrid
        if (points(i)<points(i-1)) then
          write(out,*) 'sort error at ',i,' value ', points(i)
          if (ngrid.le.100) then
            write(out,*) points(:)
          else
            write(out,*) points(i-3:i+3)
          endif
          stop 'sorting error in quicksort'
        endif
      enddo
      
      
      thresh = p005
      rint = zero
      do i = ngrid, 1, -1
        rint = rint + points(i) * dV
        if (rint > thresh) then
          write (out,'(1x,f15.6,3x,f5.2,f18.10)') points(i), thresh, rint
          thresh = thresh + p005
        end if
      end do
      write (out,*) 'final: ', rint
      
      deallocate(points)
      
    end if  ! cube integral > 0
    
  case ('dmp')
    ! dump the data to a file to be read by Mathematica
    write (out,*) 'Dumping data set'
    output = cubvalue1
    call dump_data(cub1,err,res1,npts,output,startp, &
      vectr(:,1), vectr(:,2), vectr(:,3))
    
  end select
  
  
  ! all done. Clean up, and exit gracefully
  
  close (cub1, status='keep')
  if (math.eq.'mix') close (cub2, status='keep')
  
  deallocate(cubvalue1, output, nuctyp, xyznuc, qtch)
  if (math /= 'fac' .and. math /= 'dmp') deallocate (cubvalue2)
  
  stop 'normal termination'
  
  ! =======================================================================
  
contains
  
  ! =======================================================================
  
  subroutine open_cube_read_header(iu,cub,nnuc,startp,npts,vectr)
    
    integer(KINT), intent(in) :: iu
    character(LCHARS), intent(in) :: cub
    integer(KINT), intent(out) :: nnuc, npts(3)
    real(KREAL) :: vectr(3,3), startp(3)
    integer(KINT) :: ii
    
    open (iu,file=trim(cub), status='old', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'error: file '//trim(cub)//' does not exist'
      stop 'error termination'
    end if
    
    ! read grid section and molecule for cubfile1
    
    read (iu, *, err=666, end=666)! Molecule title
    read (iu, *, err=666, end=666)! Field title
    read (iu, *, err=666, end=666) nnuc, startp(1), startp(2), startp(3)
    read (iu, *, err=666, end=666) npts(1), vectr(1,1), vectr(2,1), vectr(3,1)
    read (iu, *, err=666, end=666) npts(2), vectr(1,2), vectr(2,2), vectr(3,2)
    read (iu, *, err=666, end=666) npts(3), vectr(1,3), vectr(2,3), vectr(3,3)
    
    if (iu .eq. inp1) then
      allocate(nuctyp(nnuc))
      allocate(xyznuc(3,nnuc),qtch(nnuc))
      
      do ii = 1, nnuc
        read (iu, *, err=666, end=666) nuctyp(ii), qtch(ii), xyznuc(1,ii), &
          xyznuc(2,ii), xyznuc(3,ii)
      end do ! ii
      
    else
      do ii = 1, nnuc
        read (iu, *, err=666, end=666)
      end do
    end if
    
    return
    
666 write (out,*) 'error reading grid or XYZ from cube file '//trim(cub)
    stop 'aborting'
  end subroutine open_cube_read_header
  
  ! -------------------------------------------------------------------------
  
  subroutine write_cube_header (iu)
    
    integer(KINT), intent(in) :: iu
    
    character(LCHARS) :: fmt
    integer(KINT) :: ii
    
    fmt='(i5,3(f12.6))'
    !Write Header
    write (iu,*) '1'
    write (iu,*) 'Total Density'
    write (iu,trim(fmt)) nnuc, startp(1), startp(2), startp(3)
    write (iu,trim(fmt)) npts(1), vectr(1,1), vectr(2,1), vectr(3,1)
    write (iu,trim(fmt)) npts(2), vectr(1,2), vectr(2,2), vectr(3,2)
    write (iu,trim(fmt)) npts(3), vectr(1,3), vectr(2,3), vectr(3,3)
    
    !Write geometry
    fmt='(i5,4(f12.6))'
    do ii = 1, nnuc
      write (iu,trim(fmt)) nuctyp(ii), qtch(ii), xyznuc(1,ii), &
        xyznuc(2,ii), xyznuc(3,ii)
    end do
    
    return
  end subroutine write_cube_header
  
  
  ! -------------------------------------------------------------------------
  
end program cubcub

! =============================================================================

subroutine write_data (iu, err, cub, npts, data)

  use types
  
  implicit none
  
  integer(KINT), intent(in) :: iu, err, npts(3)
  real(KREAL), intent(in) :: data(npts(1), npts(2), npts(3))
  character(LCHARS), intent(in) :: cub
  
  integer(KINT) :: ios, i, j, k
  character(LCHARS) :: fmt
  
  fmt='(6(e13.5))'
  
  do i = 1, npts(1)
    do j = 1, npts(2)
      write (iu, trim(fmt), iostat=ios) (data(i,j,k), k= 1,npts(3))
      if (ios /= 0) then
        write (err,*) 'i, j = ',i,j
        write (err,*) 'error writing data to file '//trim(cub)
        stop 'error termination'
      end if
    end do
  end do
  
  return
end subroutine write_data

! -------------------------------------------------------------------------

subroutine dump_data (iu, err, cub, npts, data, startp, v1, v2, v3)

  use types

  implicit none
  
  integer(KINT), intent(in) :: iu, err, npts(3)
  real(KREAL), intent(in) :: data(npts(1), npts(2), npts(3))
  real(KREAL), intent(in) :: v1(3), v2(3), v3(3), startp(3)
  character(LCHARS), intent(in) :: cub
  
  integer(KINT) :: ios, i, j, k
  real(KREAL) :: x,y,z
  !character(LCHARS) :: fmt
  
  !fmt='('{',e13.5,'},{',e13.5))'
  
  do i = 1, npts(1)
    do j = 1, npts(2)
      do k = 1, npts(3)
        x = startp(1) + (i-1)*v1(1) + (j-1)*v2(1) + (k-1)*v3(1)
        y = startp(2) + (i-1)*v1(2) + (j-1)*v2(2) + (k-1)*v3(2)
        z = startp(3) + (i-1)*v1(3) + (j-1)*v2(3) + (k-1)*v3(3)
        write (iu, '(4e13.5)', iostat=ios) &
          x,y,z,data(i,j,k)
        if (ios /= 0) then
          write (err,*) 'i, j, k = ',i,j,k
          write (err,*) 'error writing data to file '//trim(cub)
          stop 'error termination'
        end if
      end do
    end do
  end do
  
  return
end subroutine dump_data

! -------------------------------------------------------------------------

subroutine volume_element(v1, v2, v3, dV)

  use types

  implicit none
  
  real(KREAL), intent(in) :: v1(3), v2(3), v3(3)
  real(KREAL), intent(out) :: dV
  real(KREAL) :: vx(3)
  real(KREAL), parameter :: zero=0.0d0
  
  dV = zero
  vx = zero
  
  vx(1) = v2(2) * v3(3) - v2(3) * v3(2)
  vx(2) = v2(3) * v3(1) - v2(1) * v3(3)
  vx(3) = v2(1) * v3(2) - v2(2) * v3(1)
  
  dV = v1(1) * vx(1) + v1(2) * vx(2) + v1(3) * vx(3)
  
  return
  
end subroutine volume_element


! =============================================================================

subroutine dumbsortr(n, a)

  ! a very simple sorting routine for a real array a(n)
  ! the code is dumb and scales as n**2
  ! NOT TESTED

  use types

  integer(KINT), intent(in) :: n
  real(KREAL), intent(inout) :: a(n)
  integer(KINT) :: i, j, k
  real(KREAL) :: rmin, swap
  
  do i = 1,n-1
    write (6,*) 'sort: i =',i
    rmin = a(i)
    k = 0
    do j = i+1, n
      if (a(j).lt.rmin) then
        k = j
        rmin = a(j)
      end if
    end do ! j = i+1,n
    ! if k.ne.0 we found a value a(k>i) < a(i) so we swap them
    if (k.ne.0) then
      swap = a(i)
      a(i) = a(k)
      a(k) = swap
    end if
  end do ! i = 1,n-1

end subroutine dumbsortr

     

       


   
  
