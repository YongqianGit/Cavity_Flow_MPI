
!##############################################################################
! This code solves for the viscous flow in a lid-driven cavity
!##############################################################################
module set_precision ! Sets the default precision for floating point numbers

  implicit none
  
  integer, parameter :: Sngl = Selected_Real_Kind(6 , 37)  
  integer, parameter :: Dbl  = Selected_Real_Kind(15, 307)
  integer, parameter :: Quad = Selected_Real_Kind(33, 4931)
  integer, parameter :: Prec = Dbl ! Set default precision to double precision

end module set_precision

  
!##############################################################################
module set_constants

  use set_precision, only : Prec

  implicit none

  real(Prec), parameter :: zero   = 0.0_Prec
  real(Prec), parameter :: tenth  = 0.1_Prec
  real(Prec), parameter :: sixth  = 1.0_Prec/6.0_Prec
  real(Prec), parameter :: fifth  = 0.2_Prec
  real(Prec), parameter :: fourth = 0.25_Prec
  real(Prec), parameter :: third  = 1.0_Prec/3.0_Prec
  real(Prec), parameter :: half   = 0.5_Prec
  real(Prec), parameter :: one    = 1.0_Prec
  real(Prec), parameter :: two    = 2.0_Prec
  real(Prec), parameter :: three  = 3.0_Prec
  real(Prec), parameter :: four   = 4.0_Prec
  real(Prec), parameter :: six    = 6.0_Prec
!$$$$$$   real(Prec), parameter :: big    = huge(zero) ! Largest resolvable #
!$$$$$$   real(Prec), parameter :: small  = tiny(zero) ! Smallest resolvable #

end module set_constants

  
!##############################################################################
module set_inputs ! Sets the input variables for the code

  use set_precision, only : Prec
  use set_constants, only : zero, one
  
  implicit none
  
  ! Set input values
  Integer, Parameter ::   imax_g = 100            ! Number of points in the x-direction (use odd numbers only)
  Integer, Parameter ::   jmax_g = 100             ! Number of points in the y-direction (use odd numbers only)
  Integer, Parameter ::   neq = 3               ! Number of equation to be solved ( = 3: mass, x-mtm, y-mtm)
  Integer ::              nmax = 100000000         ! Maximum number of iterations
  Integer ::              iterout = 5000        ! Number of time steps between solution output
  Integer ::              imms = 0              ! Manufactured solution flag: = 1 for manuf. sol., = 0 otherwise
  Integer ::              isgs = 0              ! Symmetric Gauss-Seidel  flag: = 1 for SGS, = 0 for point Jacobi
  Integer ::              irstr = 0             ! Restart flag: = 1 for restart (file 'restart.in', = 0 for initial run
  Integer ::              ipgorder = 0          ! Order of pressure gradient: 0 = 2nd, 1 = 3rd (not needed)
  Integer ::              lim = 1               ! variable to be used as the limiter sensor (= 1 for pressure)

  Real(kind=Prec) ::      cfl  = 0.5_Prec       ! CFL number used to determine time step
  Real(kind=Prec) ::      Cx = 0.01_Prec        ! Parameter for 4th order artificial viscosity in x
  Real(kind=Prec) ::      Cy = 0.01_Prec        ! Parameter for 4th order artificial viscosity in y
  Real(kind=Prec) ::      toler = 1.e-10_Prec   ! Tolerance for iterative residual convergence
  Real(kind=Prec) ::      rkappa = 0.1_Prec     ! Time derivative preconditioning constant
  Real(kind=Prec) ::      Re = 100.0_Prec       ! Reynolds number = rho*Uinf*L/rmu
  Real(kind=Prec) ::      pinf = 0.801333844662_Prec ! Initial pressure (N/m^2) -> from MMS value at cavity center
  Real(kind=Prec) ::      uinf = 1.0_Prec            ! Lid velocity (m/s)
  Real(kind=Prec) ::      rho = one             ! Density (kg/m^3)
  Real(kind=Prec) ::      xmin = zero          ! Cavity dimensions...: minimum x location (m)
  Real(kind=Prec) ::      xmax = 0.05_Prec      !                       maximum x location (m)
  Real(kind=Prec) ::      ymin = zero         !                       maximum y location (m)
  Real(kind=Prec) ::      ymax = 0.05_Prec   !                       maximum y location (m)
  Real(kind=Prec) ::      Cx2 = 0.0_Prec        ! Coefficient for 2nd order damping (not required)
  Real(kind=Prec) ::      Cy2 = 0.0_Prec        ! Coefficient for 2nd order damping (not required)
  Real(kind=Prec) ::      fsmall = 1.e-20_Prec  ! small parameter

  ! Derived input quantities (set in subroutine 'set_derived_inputs' below)
  Real(kind=Prec) ::      rhoinv =  -99.9_Prec  ! Inverse density, 1/rho (m^3/kg)
  Real(kind=Prec) ::      rlength = -99.9_Prec  ! Characteristic length (m) [cavity width]
  Real(kind=Prec) ::      rmu =     -99.9_Prec  ! Viscosity (N*s/m^2)
  Real(kind=Prec) ::      vel2ref = -99.9_Prec  ! Reference velocity squared (m^2/s^2)
  Real(kind=Prec) ::      dx =      -99.9_Prec  ! Delta x (m)
  Real(kind=Prec) ::      dy =      -99.9_Prec  ! Delta y (m)
  real(Prec) ::           rpi =     -99.9_Prec  ! Pi = 3.14159... (defined below)
  
  contains

  !#######
  subroutine set_derived_inputs

    implicit none
  
    rhoinv = one/rho                            ! Inverse density, 1/rho (m^3/kg)
    rlength = xmax - xmin                       ! Characteristic length (m) [cavity width]
    rmu = rho*uinf*rlength/Re                   ! Viscosity (N*s/m^2)
    vel2ref = uinf*uinf                         ! Reference velocity squared (m^2/s^2)
    dx = (xmax - xmin)/float(imax_g - 1)          ! Delta x (m)
    dy = (ymax - ymin)/float(jmax_g - 1)          ! Delta y (m)
    rpi = acos(-one)                            ! Pi = 3.14159... 
  
    write(*,*) 'rho,V,L,mu,Re: ',rho,uinf,rlength,rmu,Re

  end subroutine set_derived_inputs

end module set_inputs


!##############################################################################
module variables

  use set_precision, only : Prec
  use set_inputs, only : imax_g, jmax_g, neq
  
  implicit none

  Integer ::  ninit                                                  ! Initial iteration number (used for restart file)

  Real(kind=Prec),Dimension(imax_g,jmax_g,neq) ::  u_g = -99.9_Prec        ! Solution vector [p, u, v]^T at each node
  Real(kind=Prec),Dimension(imax_g,jmax_g,neq) ::  uold_g = -99.9_Prec      ! Previous (old) solution vector
  Real(kind=Prec),Dimension(imax_g,jmax_g,neq) ::  s_g = -99.9_Prec         ! Source term
  Real(kind=Prec),Dimension(imax_g,jmax_g) ::     artviscx_g = -99.9_Prec  ! Artificial viscosity in x-direction
  Real(kind=Prec),Dimension(imax_g,jmax_g) ::     artviscy_g = -99.9_Prec  ! Artificial viscosity in y-direction
  Real(kind=Prec) ::                          rtime = -99.9_Prec     ! Variable to estimate simulation time
  Real(kind=Prec) ::                          dtmin_g = +1.e99_Prec    ! Minimum time step for a given iteration (initialized large)
end module variables

module var_mpi
  use set_precision, only : Prec
  use set_inputs, only : imax_g, jmax_g, neq
  Real(kind=Prec),Allocatable :: u(:,:,:),uold(:,:,:),s(:,:,:),dt(:,:),&
  & artviscx(:,:),artviscy(:,:),&
  & dtmin_gather(:)
  Real(kind=Prec),Dimension(neq) ::res = -99.9_Prec,res_sum=-99.9_Prec     ! Iterative residual for each equation
  Real(kind=Prec),Dimension(neq) ::resinit = -99.9_Prec,resinit_sum=-99.9_Prec ! Initial iterative residual for each equation (from iteration 1)

  Integer :: sz_x,sz_y,sz_r,imax,jmax
  Real(kind=Prec) ::                          dtmin = +1.e99_Prec    ! Minimum time step for a given iteration (initialized large)
  Integer,Allocatable :: disp(:),recvcount(:)
  Real(kind=Prec),Allocatable :: sb_p(:,:),sb_u(:,:),sb_v(:,:)
  Real(kind=Prec),Dimension(imax_g*jmax_g) :: p_gather=-99.9_Prec, &
  & u_gather=-99.9_Prec, v_gather=-99.9_Prec
end module var_mpi


  
!##############################################################################


!##############################################################################
module subroutines

  implicit none

  contains

  subroutine allocate_mpi(rank,numtasks)
    use set_precision, only : Prec
    use set_inputs, only : imax_g,jmax_g,neq
    use variables, only : u_g,uold_g,s_g
    use var_mpi,only : u,uold,s,dt,artviscx,artviscy,&
    & dtmin_gather,imax,jmax,sz_x,sz_y,sz_r,disp,recvcount,&
    & sb_p,sb_u,sb_v,dtmin_gather
    Integer :: rank, numtasks
    Integer :: i
    sz_x=int(imax_g/numtasks)
    sz_r=imax_g-(numtasks-1)*sz_x
    sz_y=jmax_g
    Allocate(dtmin_gather(numtasks))
    dtmin_gather= +1.e99_Prec
    if (rank .eq. 0) then
        Allocate(u(sz_x+2,sz_y,neq), uold(sz_x+2,sz_y,neq),s(sz_x+2,sz_y,neq),&
        & artviscx(sz_x+2,sz_y),artviscy(sz_x+2,sz_y),dt(sz_x+2,sz_y))
        Allocate(sb_p(sz_x,sz_y),sb_u(sz_x,sz_y),sb_v(sz_x,sz_y))
        imax=sz_x+2
        jmax=sz_y
    elseif (rank .eq. numtasks-1) then
        Allocate(u(sz_r+2,sz_y,neq), uold(sz_r+2,sz_y,neq),s(sz_r+2,sz_y,neq),&
        & artviscx(sz_r+2,sz_y),artviscy(sz_r+2,sz_y),dt(sz_r+2,sz_y))
        Allocate(sb_p(sz_r,sz_y),sb_u(sz_r,sz_y),sb_v(sz_r,sz_y))
        imax=sz_r+2
        jmax=sz_y
    else
        Allocate(u(sz_x+4,sz_y,neq), uold(sz_x+4,sz_y,neq),s(sz_x+4,sz_y,neq),&
        & artviscx(sz_x+4,sz_y),artviscy(sz_x+4,sz_y),dt(sz_x+4,sz_y))
        Allocate(sb_p(sz_x,sz_y),sb_u(sz_x,sz_y),sb_v(sz_x,sz_y))
        imax=sz_x+4
        jmax=sz_y
    endif
    dt= -99.9_Prec
    Allocate(disp(numtasks),recvcount(numtasks))
    do i=1,numtasks
      disp(i)=(i-1)*sz_x*sz_y
      if (i .lt. numtasks) then
            recvcount(i)=sz_x*sz_y
      elseif (i .eq. numtasks) then
            recvcount(i)=sz_r*sz_y
      endif   
    enddo
  end subroutine allocate_mpi

  subroutine distribute_mpi(rank,numtasks)
    use set_precision, only : Prec
    use set_inputs, only : imax_g,jmax_g,neq
    use variables, only : u_g,uold_g,s_g
    use var_mpi,only : u,uold,s,dt,artviscx,artviscy,&
    & dtmin_gather,imax,jmax,sz_x,sz_y,sz_r,disp,recvcount
    Integer :: rank, numtasks
    Integer :: i
    sz_x=int(imax_g/numtasks)
    sz_r=imax_g-(numtasks-1)*sz_x
    sz_y=jmax_g
    if (rank .eq. 0) then
        u(1:sz_x+2,:,:)=u_g((1+(rank)*sz_x):(rank+1)*sz_x+2, :,:)
        uold(1:sz_x+2,:,:)=uold_g((1+(rank)*sz_x):(rank+1)*sz_x+2, :,:)
        s(1:sz_x+2,:,:)=s_g((1+(rank)*sz_x):(rank+1)*sz_x+2, :,:)
    elseif (rank .eq. numtasks-1) then
        u(1:sz_r+2,:,:)=u_g((1+(rank)*sz_x-2):imax_g, :,:)
        uold(1:sz_r+2,:,:)=uold_g((1+(rank)*sz_x-2):imax_g, :,:)
        s(1:sz_r+2,:,:)=s_g((1+(rank)*sz_x-2):imax_g, :,:)
    else
        u(1:sz_x+4,:,:)=u_g((1+(rank)*sz_x-2):(rank+1)*sz_x+2, :,:)
        uold(1:sz_x+4,:,:)=uold_g((1+(rank)*sz_x-2):(rank+1)*sz_x+2, :,:)
        s(1:sz_x+4,:,:)=s_g((1+(rank)*sz_x-2):(rank+1)*sz_x+2, :,:)
    endif

  end subroutine distribute_mpi


  subroutine combine_mpi(rank,numtasks,n)
    use set_precision, only : Prec
    use var_mpi, only : sz_x, sz_y, sz_r, disp, recvcount,sb_p,sb_u,sb_v,&
    & p_gather, u_gather, v_gather,u
    use variables, only : u_g
    use set_inputs, only : imax_g, jmax_g
    include 'mpif.h'
    Integer ::i, j,rank,numtasks,n
    integer stat(MPI_STATUS_SIZE),ierr

    if (rank .eq. 0) then
      sb_p=u(1:sz_x,1:sz_y,1)
      sb_u=u(1:sz_x,1:sz_y,2)
      sb_v=u(1:sz_x,1:sz_y,3)
    elseif (rank .eq. numtasks-1) then
      sb_p=u(3:sz_r+2,1:sz_y,1)
      sb_u=u(3:sz_r+2,1:sz_y,2)
      sb_v=u(3:sz_r+2,1:sz_y,3)
    else
      sb_p=u(3:sz_x+2,1:sz_y,1)
      sb_u=u(3:sz_x+2,1:sz_y,2)
      sb_v=u(3:sz_x+2,1:sz_y,3)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    
    if (rank .lt. numtasks-1) then    
      call MPI_ALLGATHERV(sb_p,sz_x*sz_y,MPI_DOUBLE,&
     & p_gather,recvcount,disp,MPI_DOUBLE,MPI_COMM_WORLD,stat,ierr)
      call MPI_ALLGATHERV(sb_u,sz_x*sz_y,MPI_DOUBLE,&
     & u_gather,recvcount,disp,MPI_DOUBLE,MPI_COMM_WORLD,stat,ierr)
      call MPI_ALLGATHERV(sb_v,sz_x*sz_y,MPI_DOUBLE,&
     & v_gather,recvcount,disp,MPI_DOUBLE,MPI_COMM_WORLD,stat,ierr)
  elseif (rank .eq. numtasks-1) then
      call MPI_ALLGATHERV(sb_p,sz_r*sz_y,MPI_DOUBLE,&
     & p_gather,recvcount,disp,MPI_DOUBLE,MPI_COMM_WORLD,stat,ierr)
      call MPI_ALLGATHERV(sb_u,sz_r*sz_y,MPI_DOUBLE,&
     & u_gather,recvcount,disp,MPI_DOUBLE,MPI_COMM_WORLD,stat,ierr)
      call MPI_ALLGATHERV(sb_v,sz_r*sz_y,MPI_DOUBLE,&
     & v_gather,recvcount,disp,MPI_DOUBLE,MPI_COMM_WORLD,stat,ierr) 
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    
    do j=1,numtasks !reshape message from all procs
        if (j .lt. numtasks) then
            u_g(1+(j-1)*sz_x:j*sz_x,:,1)=&
            & reshape(p_gather(1+sz_x*sz_y*(j-1):j*sz_x*sz_y),(/sz_x,sz_y/))
            u_g(1+(j-1)*sz_x:j*sz_x,:,2)=&
            & reshape(u_gather(1+sz_x*sz_y*(j-1):j*sz_x*sz_y),(/sz_x,sz_y/))
            u_g(1+(j-1)*sz_x:j*sz_x,:,3)=&
            & reshape(v_gather(1+sz_x*sz_y*(j-1):j*sz_x*sz_y),(/sz_x,sz_y/))
        else 
            u_g(1+(j-1)*sz_x:imax_g,:,1)=&
            & reshape(p_gather(1+sz_x*sz_y*(j-1):imax_g*jmax_g),(/sz_r,sz_y/))
            u_g(1+(j-1)*sz_x:imax_g,:,2)=&
            & reshape(u_gather(1+sz_x*sz_y*(j-1):imax_g*jmax_g),(/sz_r,sz_y/))
            u_g(1+(j-1)*sz_x:imax_g,:,3)=&
            & reshape(v_gather(1+sz_x*sz_y*(j-1):imax_g*jmax_g),(/sz_r,sz_y/))
        endif
    enddo
  end subroutine combine_mpi


  !#######
  subroutine output_file_headers 

  use set_inputs, only : imms

  ! Note: The vector of primitive variables is: 
  !               u = [p, u, v]^T
  ! Set up output files (history and solution)      
  open(30,file='history_mpi.dat',status='unknown')
  write(30,*) 'TITLE = "Cavity Iterative Residual History"'
  write(30,*) 'variables="Iteration""Time(s)""Res1""Res2""Res3"'
  
  open(40,file='cavity_mpi.dat',status='unknown')
  write(40,*) 'TITLE = "Cavity Field Data"'

    write(40,*) 'variables="x(m)""y(m)""p(N/m^2)""u(m/s)""v(m/s)"'

  
  ! Header for Screen Output
  write(*,*) '   Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum'

  end subroutine output_file_headers

!#######
!  Subroutine initial

  subroutine initial
  
  use set_precision, only : Prec
  use set_constants, only : zero, one
  use set_inputs, only : irstr, imax_g, jmax_g, neq, uinf, pinf
  use variables, only : ninit, rtime, u_g, s_g, uold_g !YY add uold
  use var_mpi, only : u,s,uold,resinit,res

  implicit none

  Integer ::              i                    ! i index (x direction)
  Integer ::              j                    ! j index (y direction)
  Integer ::              k                    ! k index (# of equations)

  Real(kind=Prec) ::      x = -99.9_Prec       ! Temporary variable for x location
  Real(kind=Prec) ::      y = -99.9_Prec       ! Temporary variable for y location

  ! This subroutine sets inital conditions in the cavity 

  ! Note: The vector of primitive variables is: 
  !               u = [p, u, v]^T

  if(irstr.eq.0) then  ! Starting run from scratch
    
    ninit = 1          ! set initial iteration to one
    rtime = zero       ! set initial time to zero
    do k = 1, neq
      resinit(k) = one
    enddo
    do i = 1, imax_g
    do j = 1, jmax_g
      uold_g(i,j,1)=pinf !YY
      uold_g(i,j,2)=zero !YY
      uold_g(i,j,3)=zero !YY
      u_g(i,j,1) = pinf
      u_g(i,j,2) = zero
      u_g(i,j,3) = zero
      s_g(i,j,1) = zero
      s_g(i,j,2) = zero
      s_g(i,j,3) = zero
    enddo
    u_g(i,jmax_g,2) = uinf ! Initialize lid (top) to freestream velocity
    enddo
    
  else
    
    write(*,*) 'ERROR: irstr must equal 0 or 1!'
    stop
    
  endif

  return

  end subroutine initial
 
!#######
!  Subroutine set_boundary_conditions

  subroutine set_boundary_conditions
  
  use set_inputs, only : imms
      
  implicit none

  ! This subroutine determines the appropriate BC routines to call

      call bndry


  end subroutine set_boundary_conditions
     
!#######
!  Subroutine bndry

  subroutine bndry
  
  use set_precision, only : Prec
  use set_constants, only : zero, one, two, half
  use set_inputs, only : imax_g, jmax_g, neq, xmax,&
  & xmin, ymax, ymin, uinf,& !!YY
  & imax_g,jmax_g,xmax,xmin,ymax,ymin
  use variables, only : u_g
      
  implicit none

  Integer ::              i                    ! i index (x direction)
  Integer ::              j                    ! j index (y direction)
  Integer ::              k                    ! k index (# of equations)

  Real(kind=Prec) ::      x = -99.9_Prec       ! Temporary variable for x location
  Real(kind=Prec) ::      y = -99.9_Prec       ! Temporary variable for y location

! This applies the cavity boundary conditions
 

!**************************************************************
! Side Walls
  do j = 2, jmax_g-1
    y = (ymax - ymin)*float(j - 1)/float(jmax_g - 1)
    i = 1
    x = xmin !left wall
    do k = 1,neq
      u_g(i,j,k) = zero
    enddo
    u_g(1,j,1) = two*u_g(2,j,1) - u_g(3,j,1)   ! 2nd Order BC (pressure term)
!$$$$$$     u(1,j,1) = u(2,j,1)                  ! 1st Order BC

    i=imax_g
    x = xmax
    do k = 1,neq
      u_g(i,j,k) = zero
    enddo
    u_g(imax_g,j,1) = two*u_g(imax_g-1,j,1) - u_g(imax_g-2,j,1)   ! 2nd Order BC(pressure)
!$$$$$$     u(imax,j,1) = u(imax-1,j,1)                       ! 1st Order BC
  enddo

! Top/Bottom Walls
  do i = 1, imax_g 
    x = (xmax - xmin)*float(i - 1)/float(imax_g - 1)
    j = 1
    y = ymin
    do k = 1,neq
      u_g(i,j,k) = zero
    enddo
    u_g(i,1,1) = two*u_g(i,2,1) - u_g(i,3,1)   ! 2nd Order BC(pressure)
!$$$$$$     u(i,1,1) = u(i,2,1)                  ! 1st Order BC

    j = jmax_g
    y = ymax
      u_g(i,j,2) = uinf
      u_g(i,j,3) = zero
    u_g(i,jmax_g,1) = two*u_g(i,jmax_g-1,1) - u_g(i,jmax_g-2,1)   ! 2nd Order BC(presure)
!$$$$$$     u(i,jmax,1) = u(i,jmax-1,1)                       ! 1st Order
!**************************************************************
  enddo
  return
  end subroutine bndry


!#######
  subroutine compute_time_step(n)

  use set_precision, only : Prec
  use set_constants, only : one, two, four, half, fourth
  use set_inputs, only : vel2ref, rmu, rho, dx, dy, cfl, rkappa
  use var_mpi, only : u,imax,jmax,dtmin,dt

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              n
  Real(kind=Prec) ::      dtvisc = -99.9_Prec      ! Viscous time step stability criteria (constant over domain)
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Local velocity squared
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared paramete for time derivative preconditioning
  Real(kind=Prec) ::      lambda_x = -99.9_Prec    ! Max absolute value eigenvalue in (x,t)
  Real(kind=Prec) ::      lambda_y = -99.9_Prec    ! Max absolute value eigenvalue in (y,t)
  Real(kind=Prec) ::      lambda_max = -99.9_Prec  ! Max absolute value eigenvalue (used in convective time step computation)
  Real(kind=Prec) ::      dtconv = -99.9_Prec      ! Local convective time step restriction

  dt=-99.9_Prec
!**************************************************************
!  print*,u(:,:,2)
  do i=1,imax
    do j=1,jmax
      beta2=max((u(i,j,2)**2+u(i,j,3)**2),vel2ref*rkappa)
      lambda_x=half*(dabs(u(i,j,2))+dsqrt(u(i,j,2)**2+four*beta2))
      lambda_y=half*(dabs(u(i,j,3))+dsqrt(u(i,j,3)**2+four*beta2))
      lambda_max=max(lambda_x,lambda_y)
      dtconv=min(dx,dy)/lambda_max
      dtvisc=dx*dy*rho/(four*rmu)
      dt(i,j)=min(dtconv,dtvisc)*cfl
!      dt(i,j)=dtconv*cfl
    end do
  end do 
          dtmin=minval(dt(:,:))
!**************************************************************
  end subroutine compute_time_step

  !#######
  subroutine compute_source_terms

  use set_precision, only : Prec
!  use set_constants, only : zero, one, two, four, half, fourth
  use set_inputs, only : imax_g, jmax_g, imms, xmax, xmin, ymax, ymin
  use variables, only : s_g

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  s_g=0.0_Prec

  end subroutine compute_source_terms

  !#######
  subroutine Compute_Artificial_Viscosity(rank,numtasks,n)

  use set_precision, only : Prec
  use set_constants, only : zero, one, two, four, six, half, fourth
  use set_inputs, only : imax_g, jmax_g, lim, rho, dx, dy, Cx, Cy, Cx2, Cy2, &
           &             fsmall, vel2ref, rkappa
  use var_mpi, only : u, artviscx, artviscy,imax,jmax

  implicit none
  include 'mpif.h'
  Integer stat(MPI_STATUS_SIZE),ierr
  Integer ::              rank,numtasks,n
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Local velocity squared
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared paramete for time derivative preconditioning
  Real(kind=Prec) ::      lambda_x = -99.9_Prec    ! Max absolute value e-value in (x,t)
  Real(kind=Prec) ::      lambda_y = -99.9_Prec    ! Max absolute value e-value in (y,t)
  Real(kind=Prec) ::      d4pdx4 = -99.9_Prec      ! 4th derivative of pressure w.r.t. x
  Real(kind=Prec) ::      d4pdy4 = -99.9_Prec      ! 4th derivative of pressure w.r.t. y
  Real(kind=Prec),Dimension(imax,jmax) :: lambda_x_tmp, lambda_y_tmp
  Real(kind=Prec),Dimension(numtasks) :: lambda_x_gather, lambda_y_gather
!**************************************************************
  lambda_x_tmp=zero
  lambda_y_tmp=zero
  artviscx=zero
  artviscy=zero
  do i=3,imax-2
    do j=3,jmax-2
      beta2=max((u(i,j,2)**2+u(i,j,3)**2),vel2ref*rkappa)
      lambda_x_tmp(i,j)=half*(dabs(u(i,j,2))+dsqrt(u(i,j,2)**2+four*beta2))
      lambda_y_tmp(i,j)=half*(dabs(u(i,j,3))+dsqrt(u(i,j,3)**2+four*beta2))
            d4pdx4=u(i+2,j,1)-four*u(i+1,j,1)+six*u(i,j,1)-four*u(i-1,j,1)+u(i-2,j,1)
            d4pdy4=u(i,j+2,1)-four*u(i,j+1,1)+six*u(i,j,1)-four*u(i,j-1,1)+u(i,j-2,1)
      artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
      artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4
    end do
  end do

  do j=3,jmax-2
      i=2 ! left wall -- YY
            d4pdx4=u(i,j,1)-four*u(i+1,j,1)+six*u(i+2,j,1)-four*u(i+3,j,1)+u(i+4,j,1) !one side for d4pdx--left wall--YY
            d4pdy4=u(i,j+2,1)-four*u(i,j+1,1)+six*u(i,j,1)-four*u(i,j-1,1)+u(i,j-2,1) !center for d4pdy--left wall--YY
            artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
            artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4
      
      i=imax-1 !right wall -- YY
            d4pdx4=u(i,j,1)-four*u(i-1,j,1)+six*u(i-2,j,1)&
                  & -four*u(i-3,j,1)+u(i-4,j,1)   !one side for right wall--right wal--YY
            d4pdy4=u(i,j+2,1)-four*u(i,j+1,1)+six*u(i,j,1)-four*u(i,j-1,1)+u(i,j-2,1) !center for d4pdy--right wall--YY
            artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
            artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4
  end do

  do i=3,imax-2
      j=2 ! bottom wall -- YY
            d4pdx4=u(i+2,j,1)-four*u(i+1,j,1)+six*u(i,j,1)-four*u(i-1,j,1)+u(i-2,j,1) !center for d4pdx--bottom wall--YY
            d4pdy4=u(i,j,1)-four*u(i,j+1,1)+six*u(i,j+2,1)-four*u(i,j+3,1)+u(i,j+4,1) !one side for d4pdy--bottom wall--YY
            artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
            artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4

      j=jmax-1
            d4pdx4=u(i+2,j,1)-four*u(i+1,j,1)+six*u(i,j,1)-four*u(i-1,j,1)+u(i-2,j,1) !center for d4pdx--top wall--YY
            d4pdy4=u(i,j,1)-four*u(i,j-1,1)+six*u(i,j-2,1)&
                  & -four*u(i,j-3,1)+u(i,j-4,1)  !one side for d4pdy--top wall--YY
            artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
            artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4
  end do

  i=2;j=2 !left bottom corner--one side for both x and y --YY
  d4pdx4=u(i,j,1)-four*u(i+1,j,1)+six*u(i+2,j,1)-four*u(i+3,j,1)+u(i+4,j,1)
  d4pdy4=u(i,j,1)-four*u(i,j+1,1)+six*u(i,j+2,1)-four*u(i,j+3,1)+u(i,j+4,1)
      artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
      artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4

  i=imax-1;j=2 !right bottom corner--one side for both x and y --YY
  d4pdx4=u(i,j,1)-four*u(i-1,j,1)+six*u(i-2,j,1)-four*u(i-3,j,1)+u(i-4,j,1)
  d4pdy4=u(i,j,1)-four*u(i,j+1,1)+six*u(i,j+2,1)-four*u(i,j+3,1)+u(i,j+4,1)
      artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
      artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4

  i=2;j=jmax-1 !left top corner--one side for both x and y--YY
  d4pdx4=u(i,j,1)-four*u(i+1,j,1)+six*u(i+2,j,1)-four*u(i+3,j,1)+u(i+4,j,1)
  d4pdy4=u(i,j,1)-four*u(i,j-1,1)+six*u(i,j-2,1)-four*u(i,j-3,1)+u(i,j-4,1)
      artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
      artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4

  i=imax-1;j=jmax-1 !right top corner--one side for both x and y --YY
  d4pdx4=u(i,j,1)-four*u(i-1,j,1)+six*u(i-2,j,1)-four*u(i-3,j,1)+u(i-4,j,1)
  d4pdy4=u(i,j,1)-four*u(i,j-1,1)+six*u(i,j-2,1)-four*u(i,j-3,1)+u(i,j-4,1)
      artviscx(i,j)=-Cx/(beta2*dx)*d4pdx4
      artviscy(i,j)=-Cy/(beta2*dy)*d4pdy4

  lambda_x=maxval(lambda_x_tmp)
  lambda_y=maxval(lambda_y_tmp)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)  
  call MPI_ALLGATHER(lambda_x,1,MPI_DOUBLE,lambda_x_gather,1,MPI_DOUBLE, &
      & MPI_COMM_WORLD,stat,ierr)
  call MPI_ALLGATHER(lambda_y,1,MPI_DOUBLE,lambda_y_gather,1,MPI_DOUBLE, &
      & MPI_COMM_WORLD,stat,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)  
  artviscx=artviscx*maxval(lambda_x_gather)
  artviscy=artviscy*maxval(lambda_y_gather)
!**************************************************************



  end subroutine Compute_Artificial_Viscosity

  !#######


  subroutine point_Jacobi

  use set_precision, only : Prec
  use set_constants, only :  two, three, six, half
  use set_inputs, only : imax_g, jmax_g, ipgorder, rho, rhoinv, dx, dy, rkappa, &
                       & xmax, xmin, ymax, ymin, rmu, vel2ref
  use var_mpi, only : u, uold, artviscx, artviscy, dt, s, imax, jmax

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  integer :: n

  Real(kind=Prec) ::      dpdx = -99.9_Prec        ! First derivative of pressure w.r.t. x
  Real(kind=Prec) ::      dudx = -99.9_Prec        ! First derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      dvdx = -99.9_Prec        ! First derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      dpdy = -99.9_Prec        ! First derivative of pressure w.r.t. y
  Real(kind=Prec) ::      dudy = -99.9_Prec        ! First derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      dvdy = -99.9_Prec        ! First derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      d2udx2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      d2vdx2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      d2udy2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      d2vdy2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared parameter for time derivative preconditioning
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Velocity squared 


  ! Point Jacobi method
!**************************************************************
  do i=2,imax-1
    do j=2,jmax-1
      beta2=max((u(i,j,2)**2+u(i,j,3)**2),rkappa*vel2ref)
      u(i,j,1)=u(i,j,1)-beta2*dt(i,j)*( &
    &   rho*(u(i+1,j,2)-u(i-1,j,2))/(two*dx) &
    &  +rho*(u(i,j+1,3)-u(i,j-1,3))/(two*dy) &
    &  -s(i,j,1) -artviscx(i,j) - artviscy(i,j))

      u(i,j,2)=u(i,j,2)-dt(i,j)/rho*( &
    &  rho*u(i,j,2)*(u(i+1,j,2)-u(i-1,j,2))/(two*dx)  &
    & +rho*u(i,j,3)*(u(i,j+1,2)-u(i,j-1,2))/(two*dy)  &
    & +(u(i+1,j,1)-u(i-1,j,1))/(two*dx) &
    & -rmu*(u(i+1,j,2)-two*u(i,j,2)+u(i-1,j,2))/(dx**2) &
    & -rmu*(u(i,j+1,2)-two*u(i,j,2)+u(i,j-1,2))/(dy**2) &
    & - s(i,j,2)  )

      u(i,j,3)=u(i,j,3)-dt(i,j)/rho*( &
    &  rho*u(i,j,2)*(u(i+1,j,3)-u(i-1,j,3))/(two*dx)  &
    & +rho*u(i,j,3)*(u(i,j+1,3)-u(i,j-1,3))/(two*dy)  &
    & +(u(i,j+1,1)-u(i,j-1,1))/(two*dy) &
    & -rmu*(u(i+1,j,3)-two*u(i,j,3)+u(i-1,j,3))/(dx**2) &
    & -rmu*(u(i,j+1,3)-two*u(i,j,3)+u(i,j-1,3))/(dy**2) &
    & - s(i,j,3)  )
    end do
  end do 
!**************************************************************
  end subroutine point_Jacobi

  !#######

  
  subroutine SGS_forward_sweep

  use set_precision, only : Prec
  use set_constants, only :  two, three, six, half
  use set_inputs, only : imax_g, jmax_g, ipgorder, rho, rhoinv, dx, dy, rkappa, &
                       & xmax, xmin, ymax, ymin, rmu, vel2ref
  use var_mpi, only : u, artviscx, artviscy, dt, s, imax, jmax

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      dpdx = -99.9_Prec        ! First derivative of pressure w.r.t. x
  Real(kind=Prec) ::      dudx = -99.9_Prec        ! First derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      dvdx = -99.9_Prec        ! First derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      dpdy = -99.9_Prec        ! First derivative of pressure w.r.t. y
  Real(kind=Prec) ::      dudy = -99.9_Prec        ! First derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      dvdy = -99.9_Prec        ! First derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      d2udx2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      d2vdx2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      d2udy2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      d2vdy2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared parameter for time derivative preconditioning
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Velocity squared 


  ! Symmetric Gauss-Siedel: Forward Sweep

!**************************************************************
  do i=2,imax-1
    do j=2,jmax-1
      beta2=max((u(i,j,2)**2+u(i,j,3)**2),rkappa*vel2ref)
      u(i,j,1)=u(i,j,1)-beta2*(half*dt(i,j))*( &
    &   rho*(u(i+1,j,2)-u(i-1,j,2))/(two*dx) &
    &  +rho*(u(i,j+1,3)-u(i,j-1,3))/(two*dy) &
    &  -s(i,j,1) -artviscx(i,j) - artviscy(i,j))

      u(i,j,2)=u(i,j,2)-(half*dt(i,j))/rho*( &
    &  rho*u(i,j,2)*(u(i+1,j,2)-u(i-1,j,2))/(two*dx)  &
    & +rho*u(i,j,3)*(u(i,j+1,2)-u(i,j-1,2))/(two*dy)  &
    & +(u(i+1,j,1)-u(i-1,j,1))/(two*dx) &
    & -rmu*(u(i+1,j,2)-two*u(i,j,2)+u(i-1,j,2))/(dx**2) &
    & -rmu*(u(i,j+1,2)-two*u(i,j,2)+u(i,j-1,2))/(dy**2) &
    & - s(i,j,2)  )

      u(i,j,3)=u(i,j,3)-(half*dt(i,j))/rho*( &
    &  rho*u(i,j,2)*(u(i+1,j,3)-u(i-1,j,3))/(two*dx)  &
    & +rho*u(i,j,3)*(u(i,j+1,3)-u(i,j-1,3))/(two*dy)  &
    & +(u(i,j+1,1)-u(i,j-1,1))/(two*dy) &
    & -rmu*(u(i+1,j,3)-two*u(i,j,3)+u(i-1,j,3))/(dx**2) &
    & -rmu*(u(i,j+1,3)-two*u(i,j,3)+u(i,j-1,3))/(dy**2) &
    & - s(i,j,3)  )
    end do
  end do
!**************************************************************



  end subroutine SGS_forward_sweep

  !#######
  subroutine SGS_backward_sweep(rank,numtasks,n)

  use set_precision, only : Prec
  use set_constants, only :  two, three, six, half
  use set_inputs, only : imax_g, jmax_g, ipgorder, rho, rhoinv, dx, dy, rkappa, &
                       &  xmax, xmin, ymax, ymin, rmu, vel2ref
  use var_mpi, only : u, artviscx, artviscy, dt, s, imax, jmax

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer :: rank,numtasks,n
  Real(kind=Prec) ::      dpdx = -99.9_Prec        ! First derivative of pressure w.r.t. x
  Real(kind=Prec) ::      dudx = -99.9_Prec        ! First derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      dvdx = -99.9_Prec        ! First derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      dpdy = -99.9_Prec        ! First derivative of pressure w.r.t. y
  Real(kind=Prec) ::      dudy = -99.9_Prec        ! First derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      dvdy = -99.9_Prec        ! First derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      d2udx2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      d2vdx2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      d2udy2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      d2vdy2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared parameter for time derivative preconditioning
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Velocity squared 

  ! Symmetric Gauss-Siedel: Backward Sweep

!**************************************************************
  do i=imax-1,2,-1
    do j=jmax-1,2,-1
      beta2=max((u(i,j,2)**2+u(i,j,3)**2),rkappa*vel2ref)
      u(i,j,1)=u(i,j,1)-beta2*(half*dt(i,j))*( &
    &   rho*(u(i+1,j,2)-u(i-1,j,2))/(two*dx) &
    &  +rho*(u(i,j+1,3)-u(i,j-1,3))/(two*dy) &
    &  -s(i,j,1) -artviscx(i,j) - artviscy(i,j))

      u(i,j,2)=u(i,j,2)-(half*dt(i,j))/rho*( &
    &  rho*u(i,j,2)*(u(i+1,j,2)-u(i-1,j,2))/(two*dx)  &
    & +rho*u(i,j,3)*(u(i,j+1,2)-u(i,j-1,2))/(two*dy)  &
    & +(u(i+1,j,1)-u(i-1,j,1))/(two*dx) &
    & -rmu*(u(i+1,j,2)-two*u(i,j,2)+u(i-1,j,2))/(dx**2) &
    & -rmu*(u(i,j+1,2)-two*u(i,j,2)+u(i,j-1,2))/(dy**2) &
    & - s(i,j,2)  )

      u(i,j,3)=u(i,j,3)-(half*dt(i,j))/rho*( &
    &  rho*u(i,j,2)*(u(i+1,j,3)-u(i-1,j,3))/(two*dx)  &
    & +rho*u(i,j,3)*(u(i,j+1,3)-u(i,j-1,3))/(two*dy)  &
    & +(u(i,j+1,1)-u(i,j-1,1))/(two*dy) &
    & -rmu*(u(i+1,j,3)-two*u(i,j,3)+u(i-1,j,3))/(dx**2) &
    & -rmu*(u(i,j+1,3)-two*u(i,j,3)+u(i,j-1,3))/(dy**2) &
    & - s(i,j,3)  )
    end do
  end do

!**************************************************************



  end subroutine SGS_backward_sweep

  !#######
  subroutine pressure_rescaling

  use set_precision, only : Prec
  use set_inputs, only : imax_g, jmax_g, imms, xmax, xmin, ymax, ymin, rlength, pinf
  use variables, only : u_g,uold_g 

  implicit none
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              iref                     ! i index location of pressure rescaling point
  Integer ::              jref                     ! j index location of pressure rescaling point

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location
  Real(kind=Prec) ::      deltap = -99.9_Prec ,deltap_old=-99.9_Prec     ! Temporary variable: y location

  iref = (imax_g-1)/2+1     ! Set reference pressure to center of cavity
  jref = (jmax_g-1)/2+1
 
    deltap = u_g(iref,jref,1) - pinf ! Reference pressure 
    deltap_old=uold_g(iref,jref,1) - pinf
  do i = 1, imax_g
  do j = 1, jmax_g
    u_g(i,j,1) = u_g(i,j,1) - deltap
    uold_g(i,j,1)=uold_g(i,j,1) - deltap_old
  enddo
  enddo

  end subroutine pressure_rescaling

  !#######
  subroutine check_iterative_convergence(n,conv,rank,numtasks)

  use set_precision, only : Prec
  use set_constants, only : zero ,&
        & two !YY
  use set_inputs, only : imax_g, jmax_g, neq, fsmall, &
        & dx, dy, rmu, rho, rkappa, vel2ref, toler,iterout !YY 
  use variables, only : rtime,ninit
  use var_mpi, only : u, uold, dt, res, resinit, dtmin, &
        & artviscx, artviscy, s, imax, jmax, res_sum,resinit_sum ! YY

  implicit none
  include 'mpif.h'
  integer stat(MPI_STATUS_SIZE),ierr
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                        ! k index (over # of equations)
  Integer ::              n                        ! Current iteration number
  Integer ::              rank
  Real(kind=Prec) ::      conv                     ! Minimum of iterative residual norms from three equations
  real(kind=Prec) :: beta2
  Integer :: numtasks
  ! Compute iterative residuals to monitor iterative convergence

!**************************************************************
  conv=zero
  res=zero
  if (n .eq. ninit)then
        resinit=zero  
  endif
  do i=2,imax-1
    do j=2,jmax-1
      beta2=max((uold(i,j,2)**2+uold(i,j,3)**2),rkappa*vel2ref)
      if (n .eq. ninit)then
            resinit(1)=((u(i,j,1)-uold(i,j,1))/(beta2*dt(i,j)))**2+resinit(1)
            resinit(2)=((u(i,j,2)-uold(i,j,2))*rho/dt(i,j))**2+resinit(2)
            resinit(3)=((u(i,j,3)-uold(i,j,3))*rho/dt(i,j))**2+resinit(3)
      endif 

      
      res(1)=((u(i,j,1)-uold(i,j,1))/(beta2*dt(i,j)))**2+res(1)
      res(2)=((u(i,j,2)-uold(i,j,2))*rho/dt(i,j))**2+res(2)
      res(3)=((u(i,j,3)-uold(i,j,3))*rho/dt(i,j))**2+res(3)
    end do
  end do
      if (n .eq. ninit ) then
        resinit=sqrt(resinit/float((jmax-2)*(imax-2)))
        call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        call MPI_ALLREDUCE(resinit(1),resinit_sum(1),1,MPI_DOUBLE,MPI_SUM, &
        & MPI_COMM_WORLD,stat,ierr)
        call MPI_ALLREDUCE(resinit(2),resinit_sum(2),1,MPI_DOUBLE,MPI_SUM, &
        & MPI_COMM_WORLD,stat,ierr)
        call MPI_ALLREDUCE(resinit(3),resinit_sum(3),1,MPI_DOUBLE,MPI_SUM, &
        & MPI_COMM_WORLD,stat,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      endif
     
      do k=1,3 !if resinit is too close to 0, set it to res
        if (dabs(resinit(k)) .le. 1000.0_Prec*toler)  then
            resinit(k)=res(k)
            resinit(k)=sqrt(resinit(k)/float((jmax-2)*(imax-2)))
            call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
            call MPI_ALLREDUCE(resinit(1),resinit_sum(1),1,MPI_DOUBLE,MPI_SUM, &
            & MPI_COMM_WORLD,stat,ierr)
            call MPI_ALLREDUCE(resinit(2),resinit_sum(2),1,MPI_DOUBLE,MPI_SUM, &
            & MPI_COMM_WORLD,stat,ierr)
            call MPI_ALLREDUCE(resinit(3),resinit_sum(3),1,MPI_DOUBLE,MPI_SUM, &
            & MPI_COMM_WORLD,stat,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        endif
      enddo
      
      res=sqrt(res/float((jmax-2)*(imax-2)))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      call MPI_ALLREDUCE(res(1),res_sum(1),1,MPI_DOUBLE,MPI_SUM, &
      & MPI_COMM_WORLD,stat,ierr)
      call MPI_ALLREDUCE(res(2),res_sum(2),1,MPI_DOUBLE,MPI_SUM, &
      & MPI_COMM_WORLD,stat,ierr)
      call MPI_ALLREDUCE(res(3),res_sum(3),1,MPI_DOUBLE,MPI_SUM, &
      & MPI_COMM_WORLD,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 



      res=res_sum/resinit_sum
!**************************************************************

        conv = maxval(res) ! Place current maximum scaled iterative residual into variable 'conv'
        
        ! Write iterative residuals every 500 iterations
        if( (mod(n,iterout).eq.0).or.n.eq.ninit ) then
          if (rank .eq. 0) then
            write (30,800) n,rtime,(res(k),k=1,3)
            write(*,300) n,rtime,dtmin,(res(k),k=1,3)
            write(*,300) n,rtime,dtmin,(resinit_sum(k)/numtasks,k=1,3)
      300   format(i8,2(e11.3),3(e14.5))   
      800   format(I10,E14.6,3(E18.10))       
          write(*,*) '   Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum'
          endif
        endif  
  end subroutine check_iterative_convergence

  !#######
  subroutine write_output(n)

  use set_precision, only : Prec
  use set_inputs, only : imax_g, jmax_g, neq, xmax, xmin, ymax, ymin, &
      &                  rlength, imms
  use variables, only : u_g

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                        ! k index (over # of equations)
  Integer ::              n                        ! Current iteration number

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location
   
  ! Field output
  write(40,*) 'zone T="n = ',n,' " '
  write(40,*) 'I=',imax_g,' J=',jmax_g
  write(40,*) 'DATAPACKING=POINT'

    do j = 1, jmax_g
    do i = 1, imax_g
      x = (xmax - xmin)*float(i - 1)/float(imax_g - 1)
      y = (ymax - ymin)*float(j - 1)/float(jmax_g - 1)
      write(40,500)x,y,(u_g(i,j,k),k=1,neq)
  500 format(2(E14.6),3(E22.14))     
    enddo
    enddo
 
  ! Restart file: overwrites every 'iterout' iteration
  end subroutine write_output

 
  subroutine output_final

  use set_precision, only : Prec
  use set_inputs, only : imax_g, jmax_g, neq, xmax, xmin, ymax, ymin, &
      &                  rlength, imms
  use variables, only : u_g
  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                        ! k index (over # of equations)
  Integer ::              n                        ! Current iteration number

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location
  open(520,file='final_results_mpi.dat',status='unknown')

    write(520,*)'%x     y     p     u     v'
    do j = 1, jmax_g
    do i = 1, imax_g
      x = (xmax - xmin)*float(i - 1)/float(imax_g - 1)
      y = (ymax - ymin)*float(j - 1)/float(jmax_g - 1)
      write(520,400)x,y,(u_g(i,j,k),k=1,neq)
  400 format(2(E14.6),3(E22.14))    
    enddo
    enddo

close(520)
  end subroutine output_final

  !#######


end module subroutines


!##############################################################################
!####################### P R O G R A M  C A V I T Y ###########################
!##############################################################################
Program Cavity


  ! Note: The vector of primitive variables is: 
  !               u = [p, u, v]^T
 
  use set_precision, only : Prec
  use set_constants, only : zero
  use set_inputs, only : imax_g, jmax_g, xmax, xmin, ymax, ymin, rlength, nmax, &
                  &      neq, toler, iterout, isgs, set_derived_inputs
  use variables, only : ninit, u_g, uold_g, artviscx_g, artviscy_g, &
                  &      rtime, dtmin_g
  use var_mpi, only : dtmin,dtmin_gather,u,uold,res,resinit,artviscx,artviscy
  use subroutines  
      
  implicit none
  
include 'mpif.h'

  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                    ! k index (# of equations)
  Integer ::              n                    ! Iteration number index
 
  Real(kind=Prec) ::      conv = -99.9_Prec    ! Minimum of iterative residual norms from three equations
  integer numtasks, rank, dest, source, count, tag, ierr
  integer stat(MPI_STATUS_SIZE)
  Real(kind=Prec) :: start_t, end_t

  call cpu_time(start_t)
  
!===================MPI 2015/5=======================  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
!===================MPI 2015/5=======================  
  ! Initialize Artificial Viscosity arrays to zero (note: artviscx(i,j) and artviscy(i,j)
      artviscx_g=zero
      artviscy_g=zero
  ! Set derived input quantities
  call set_derived_inputs

  ! Set Initial Profile for u vector
  call initial
  ! Set Boundary Conditions for u
  call set_boundary_conditions
  if (rank .eq. 0) then
  ! Set up headers for output files
      call output_file_headers
  ! Write out inital conditions to solution file
      call write_output(ninit)
  endif
      
  ! Evaluate Source Terms Once at Beginning (only interior points; will be zero for standard cavity)
  call compute_source_terms


  
  call allocate_mpi(rank,numtasks)! start proc distribution
  ! ========== Main Loop ==========
  do n = ninit, nmax
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    
    call set_boundary_conditions
    call distribute_mpi(rank,numtasks)! start proc distribution
    artviscx=zero
    artviscy=zero
    ! Calculate time step   
    call compute_time_step(n)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

    call MPI_ALLGATHER(dtmin,1,MPI_DOUBLE,dtmin_gather,1,MPI_DOUBLE, &
      & MPI_COMM_WORLD,stat,ierr)  
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    dtmin_g=minval(dtmin_gather)


    ! Save u values at time level n (u and uold are 2D arrays)
    uold=u
    uold_g=u_g

    if (isgs.eq.1) then ! ==Symmetric Gauss Seidel==
      
      ! Artificial Viscosity
      call Compute_Artificial_Viscosity(rank,numtasks,n)
      
      ! Symmetric Gauss-Siedel: Forward Sweep
      call SGS_forward_sweep

      call combine_mpi(rank,numtasks,n)
   
      ! Set Boundary Conditions for u
      call set_boundary_conditions !renew BC every step

      call distribute_mpi(rank,numtasks)! start proc distribution

      ! Artificial Viscosity
      call Compute_Artificial_Viscosity(rank,numtasks,n)

   
      ! Symmetric Gauss-Siedel: Backward Sweep
      call SGS_backward_sweep(rank,numtasks,n)


      ! Set Boundary Conditions for u
 
    elseif (isgs.eq.0) then ! ==Point Jacobi==

      ! Artificial Viscosity
      call Compute_Artificial_Viscosity(rank,numtasks,n)
      
      ! Symmetric Gauss-Siedel: Forward Sweep
      call point_Jacobi
      
    else
      write(*,*) 'ERROR: isgs must be 0 or 1!!!'
      stop
    endif
 
    call combine_mpi(rank,numtasks,n)

    ! Pressure Rescaling (based on center point)
    call pressure_rescaling

    call distribute_mpi(rank,numtasks)

    ! Update the time
    rtime = rtime + dtmin_g
    ! Check iterative convergence using L2 norms of iterative residuals
    call check_iterative_convergence(n,conv,rank,numtasks)
          
          if(conv.le.toler) then
            goto 200  !exit iteration loop if converged
          else
            continue
          endif
                      
          ! Output solution and restart file every 'iterout' steps
!    if (rank .eq. 0) then                 
!          if( (mod(n,iterout).eq.0) ) then
!            call write_output(n)
!          endif
!    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
  enddo  ! ========== End Main Loop ==========
  
        write(*,*) 'Solution failed to converge in ',nmax,' iterations!!!'
          
        goto 201
          
  200  continue  ! go here once solution is converged
      
  if (rank .eq. 0) then  
        write(30,700) n,rtime,(res(k),k=1,3)
        write(*,*) 'Solution converged in ',n,' iterations!!!'
  700   format(I10,E14.6,3(E18.10))
  endif   
  201  continue
      
  if (rank .eq. 0) then  
        ! Output solution and restart file 
        call output_final
        call write_output(n)
      
        ! Close open files
        close(50)
        close(30)
        close(40)
  endif    
  call MPI_FINALIZE(ierr) !MPI finalization 2015
  if (rank .eq. 0) then
      call cpu_time(end_t)

      print '("Time = ",f6.3," seconds.")',end_t-start_t
      open(1111,file='runtime_mpi.dat',status='unknown')
      write(1111,*),start_t,end_t,end_t-start_t
      close(1111)
  endif

end Program Cavity
