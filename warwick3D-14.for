
C     https://warwick.ac.uk/research/rtp/sc/rse
C     /training/advancedmpi/02_case_study_1.pdf
C     Chris Brady & Heather Ratcliffe (warwick university)
C     compile: mpif77 -o warwick warwick.for
C     or:      mpif90 -o warwick warwick.for
C     run:     mpirun -np 16 ./warwick

      MODULE allocated

      USE MPI

      IMPLICIT NONE

      INTEGER, PARAMETER :: im = 400, jm = 100,  km = 100
      INTEGER, PARAMETER :: ndim = 3, nvec = 15, nrow = 5
      INTEGER, PARAMETER :: tag = 100

      INTEGER  :: obs(im,jm,km)
      REAL     :: vel(im,jm,km,ndim)
      REAL     :: rho(im,jm,km)

      REAL,    ALLOCATABLE :: vel_l(:,:,:,:)
      REAL,    ALLOCATABLE :: fin_l(:,:,:,:)
      REAL,    ALLOCATABLE :: fout_l(:,:,:,:)
      REAL,    ALLOCATABLE :: feq_l(:,:,:,:)
      REAL,    ALLOCATABLE :: rho_l(:,:,:)

      INTEGER :: im_l, jm_l, km_l
      INTEGER :: imin_l, imax_l
      INTEGER :: jmin_l, jmax_l
      INTEGER :: kmin_l, kmax_l

      INTEGER :: nproc, rank, cart_comm
      INTEGER :: nprocs(3), coords(3)
      INTEGER :: xmin_rank, xmax_rank
      INTEGER :: ymin_rank, ymax_rank
      INTEGER :: zmin_rank, zmax_rank

      END MODULE allocated




      SUBROUTINE gather_to_zero
      USE allocated

      INTEGER :: ierr
      REAL, DIMENSION(im,jm,km,ndim) :: red
      REAL, DIMENSION(im,jm,km) :: blu

      vel = 0.0
      red = 0.0

      rho = 0.0
      blu = 0.0

      vel(imin_l:imax_l,jmin_l:jmax_l,kmin_l:kmax_l,1:3) =
     &    vel_l(1:im_l,1:jm_l,1:km_l,1:3)

      rho(imin_l:imax_l,jmin_l:jmax_l,kmin_l:kmax_l) = 
     &    rho_l(1:im_l,1:jm_l,1:km_l)

      CALL MPI_Reduce(vel,red,im*jm*km*3,
     &     MPI_REAL,MPI_SUM,0,cart_comm,ierr)

      CALL MPI_Reduce(rho,blu,im*jm*km,
     &     MPI_REAL,MPI_SUM,0,cart_comm,ierr)

      vel = red
      rho = blu

      END SUBROUTINE




      SUBROUTINE bcs(array)
      USE allocated

      REAL,DIMENSION(0:im_l+1, 0:jm_l+1, 0:km_l+1, 1:nvec),
     &  INTENT(INOUT)::array
      INTEGER :: ierr

C     Send left strip of cells left and receive into left ghost cells
      CALL MPI_Sendrecv(
     & array(1, 1:jm_l, 1:km_l, 1:nvec),jm_l*km_l*nvec,
     & MPI_REAL,xmin_rank,tag,
     & array(0, 1:jm_l, 1:km_l, 1:nvec),jm_l*km_l*nvec,
     & MPI_REAL,xmin_rank,tag,
     & cart_comm, MPI_STATUS_IGNORE, ierr)

C     Send right strip of cells right and receive into right ghost cells
      CALL MPI_Sendrecv(
     & array(im_l,   1:jm_l, 1:km_l, 1:nvec), jm_l*km_l*nvec,
     & MPI_REAL,xmax_rank,tag,
     & array(im_l+1, 1:jm_l, 1:km_l, 1:nvec), jm_l*km_l*nvec,
     & MPI_REAL,xmax_rank,tag,
     & cart_comm, MPI_STATUS_IGNORE, ierr)

      END SUBROUTINE bcs




      SUBROUTINE setup_mpi
      USE allocated

      LOGICAL, DIMENSION(3) :: periods
      INTEGER :: ierr 

      periods = .FALSE.
      ndimpar = 1 ! num of dims to decompose

      CALL MPI_Init(ierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      CALL MPI_Dims_create(nproc, ndimpar, nprocs, ierr)

      IF (rank == 0) THEN
        print *,'Processor decomposition is ', nprocs
      ENDIF

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

      im_l = im / nprocs(1)
      jm_l = jm ! / nprocs(2)
      km_l = km ! / nprocs(3)

      CALL MPI_Cart_create(MPI_COMM_WORLD,ndimpar,
     &     nprocs,periods,.TRUE.,cart_comm,ierr)

      CALL MPI_Comm_rank(cart_comm,rank,ierr)
      CALL MPI_Cart_shift(cart_comm,0,1,xmin_rank,xmax_rank,ierr)
C      CALL MPI_Cart_shift(cart_comm,1,1,ymin_rank,ymax_rank,ierr)
C      CALL MPI_Cart_shift(cart_comm,2,1,zmin_rank,zmax_rank,ierr)
      CALL MPI_Cart_coords(cart_comm,rank,ndimpar,coords,ierr)

      imin_l = im_l *  coords(1) + 1
      imax_l = im_l * (coords(1) + 1)
      jmin_l = jm_l *  coords(2) + 1
      jmax_l = jm_l * (coords(2) + 1)
      kmin_l = km_l *  coords(3) + 1
      kmax_l = km_l * (coords(3) + 1)

      END SUBROUTINE setup_mpi



C-----start ascii VTK file
      SUBROUTINE write_ascii_vtk
      USE allocated

      INTEGER i,j,k
      WRITE(6,*) 'write_ascii_vtk()...'
      OPEN(UNIT=20, FILE='warwick.vtk')
      WRITE(20,10)'# vtk DataFile Version 2.0'
      WRITE(20,10)'sample rectilinear grid'
      WRITE(20,10)'ASCII'
      WRITE(20,10)'DATASET RECTILINEAR_GRID'
      WRITE(20,20)'DIMENSIONS ', im+1, jm+1, km+1
      WRITE(20,30)'X_COORDINATES ', im+1, ' float'
      WRITE(20,*)( real(i-1),   i=1,im+1)
      WRITE(20,30)'Y_COORDINATES ', jm+1, ' float'
      WRITE(20,*)( real(j-1),   j=1,jm+1)
      WRITE(20,30)'Z_COORDINATES ', km+1, ' float'
      WRITE(20,*)( real(k-1),   k=1,km+1)
      WRITE(20,40)'CELL_DATA ', im*jm*km
C-----scalar obs
      WRITE(20,10)'SCALARS obs int'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)(((obs(i,j,k),i=1,im),j=1,jm),k=1,km)
C-----scalar rho
      WRITE(20,10)'SCALARS rho float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)(((rho(i,j,k),i=1,im),j=1,jm),k=1,km)
C-----scalar vel
      WRITE(20,10)'VECTORS vel float'
      WRITE(20,*)(((vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),
     &              i=1,im),j=1,jm),k=1,km)
      CLOSE(20)
   10 FORMAT(A)
   20 FORMAT(A,3I5)
   30 FORMAT(A,I5,A)
   40 FORMAT(A,I9)
      END SUBROUTINE




C-----start binary VTK file
      SUBROUTINE write_binary_vtk
      USE allocated

      INTEGER i,j,k
      CHARACTER(LEN=1)  :: lf
      CHARACTER(LEN=10) :: str1,str2,str3,str4
      lf = char(10)
      WRITE(6,*) 'write_binary_vtk()...'
      OPEN(unit=20, file='warwick.vtk', form='unformatted',
     &  access='stream',status='replace',convert="big_endian")
      write(20)'# vtk DataFile Version 3.0'//lf
      write(20)'vtk output'//lf
      write(20)'BINARY'//lf
      write(20)'DATASET RECTILINEAR_GRID'//lf
      write(str1(1:10),'(i10)') im+1
      write(str2(1:10),'(i10)') jm+1
      write(str3(1:10),'(i10)') km+1
      write(str4(1:10),'(i10)') im*jm*km
      write(20)'DIMENSIONS '//str1//str2//str3//lf
      write(20)'X_COORDINATES '//str1//' float'//lf
      write(20)(real(i-1),i=1,im+1)
      write(20)'Y_COORDINATES '//str2//' float'//lf
      write(20)(real(j-1),j=1,jm+1)
      write(20)'Z_COORDINATES '//str3//' float'//lf
      write(20)(real(k-1),k=1,km+1)
      write(20)'CELL_DATA '//str4//lf
C-----obstacle
      write(20)'SCALARS obs int'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((obs(i,j,k),i=1,im),j=1,jm),k=1,km)
C-----obstacle
      write(20)'SCALARS rho float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((rho(i,j,k),i=1,im),j=1,jm),k=1,km)
C-----velocity
      write(20)'VECTORS vel float'//lf
      write(20)(((vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),
     &         i=1,im),j=1,jm),k=1,km)
      close(20)
C-----finish binary VTK file
      END SUBROUTINE





      PROGRAM serial

      USE allocated

      IMPLICIT NONE

      integer vec(nvec,ndim),col1(nrow),col2(nrow),col3(nrow)
      integer iter,ierr,i,j,k,ii,jj,kk,n,c,nexti,nextj,nextk
      integer nits,nout

      real lref,re,uLB,nuLB,omega,cu,usqr,sum2,sum3
      real cx,cy,cz,cr,dx,dy,dz,dist
      real wt(nvec),vin(ndim)

      CALL setup_mpi

      ALLOCATE( rho_l(0:im_l+1,0:jm_l+1,0:km_l+1))
      ALLOCATE( vel_l(0:im_l+1,0:jm_l+1,0:km_l+1,ndim))
      ALLOCATE( fin_l(0:im_l+1,0:jm_l+1,0:km_l+1,nvec)) 
      ALLOCATE(fout_l(0:im_l+1,0:jm_l+1,0:km_l+1,nvec))
      ALLOCATE( feq_l(0:im_l+1,0:jm_l+1,0:km_l+1,nvec)) 

      Re = 100.0         ! reynolds number
      Lref  = 10.0       ! length scale (prev Lref = cylinder radius)
      uLB   = 0.04       ! velocity in lattice units
      nuLB  = uLB*Lref/Re            ! viscosity in lattice units
      omega = 1.0 / (3.0*nuLB + 0.5) ! relaxation parameter
      vin = [uLB,0.0,0.0]   
      nits = 10000
      nout = 10

C=====initialise vectors
      vec(1, :) = (/ 1, 0, 0 /)
      vec(2, :) = (/ 0, 1, 0 /)
      vec(3, :) = (/ 0, 0, 1 /)
      vec(4, :) = (/ 1, 1, 1 /) 
      vec(5, :) = (/ 1, 1,-1 /)
      vec(6, :) = (/ 1,-1, 1 /) 
      vec(7, :) = (/ 1,-1,-1 /) 
      vec(8, :) = (/ 0, 0, 0 /)  
      vec(9, :) = (/-1, 1, 1 /)
      vec(10,:) = (/-1, 1,-1 /)
      vec(11,:) = (/-1,-1, 1 /)
      vec(12,:) = (/-1,-1,-1 /)
      vec(13,:) = (/ 0, 0,-1 /)
      vec(14,:) = (/ 0,-1, 0 /) 
      vec(15,:) = (/-1, 0, 0 /)

C=====initialise weights

      wt(1)  = 1./9.
      wt(2)  = 1./9.
      wt(3)  = 1./9.
      wt(4)  = 1./72.
      wt(5)  = 1./72.
      wt(6)  = 1./72.
      wt(7)  = 1./72.
      wt(8)  = 2./9.
      wt(9)  = 1./72.
      wt(10) = 1./72.
      wt(11) = 1./72.
      wt(12) = 1./72.
      wt(13) = 1./9.
      wt(14) = 1./9.
      wt(15) = 1./9.

      col1 = (/ 1,  4,  5,  6,   7 /)
      col2 = (/ 2,  3,  8,  13, 14 /)
      col3 = (/ 9,  10, 11, 12, 15 /)

C-----outer boundaries
c      if (imin_l.eq.1) vel_l(0,:,:,1) = 1.0 ! xmin
c      if (jmin_l.eq.1) vel_l(:,0,:,1) = 1.0 ! ymin
c      if (kmin_l.eq.1) vel_l(:,:,0,1) = 1.0 ! zmin
c      if (imax_l.eq.im) vel_l(im_l+1,:,:,1) = 1.0 ! xmax
c      if (jmax_l.eq.jm) vel_l(:,jm_l+1,:,1) = 1.0 ! ymax
c      if (kmax_l.eq.km) vel_l(:,:,km_l+1,1) = 1.0 ! zmax

C=====obstacle size / location
      cr = real(jm)/5.0 ! radius
      cx = real(im)/4.0
      cy = real(jm)/2.0
      cz = real(km)/2.0

      obs = 0

      DO k = 1, km
      DO j = 1, jm
      DO i = 1, im

        dx = (real(i)-0.5) - cx
        dy = (real(j)-0.5) - cy
        dz = (real(k)-0.5) - cz
        dist = sqrt(dx*dx + dy*dy + dz*dz)
        if (dist.LE.cr) obs(i,j,k) = 1

C-------put obstacle boundary on side walls
c        if ((j.eq.1).OR.(j.eq.jm)) obs(i,j,k) = 1
c        if ((k.eq.1).OR.(k.eq.km)) obs(i,j,k) = 1

      END DO
      END DO
      END DO


C=====initialise local velocity field
      DO k = 1, km_l
      DO j = 1, jm_l
      DO i = 1, im_l

        ii = imin_l + i-1
        jj = jmin_l + j-1
        kk = kmin_l + k-1

        if (obs(ii,jj,kk).eq.0) vel_l(i,j,k,:) = vin(:)
        if (obs(ii,jj,kk).eq.1) vel_l(i,j,k,:) = 0.0

      END DO
      END DO
      END DO


C=====equilibrium distribution function
C-----fin = equilibrium(1,u)
      DO k = 1, km_l
      DO j = 1, jm_l
      DO i = 1, im_l

          usqr = vel_l(i,j,k,1)**2
     &         + vel_l(i,j,k,2)**2
     &         + vel_l(i,j,k,3)**2

          do n = 1,nvec
            cu = vec(n,1)*vel_l(i,j,k,1)
     &         + vec(n,2)*vel_l(i,j,k,2)
     &         + vec(n,3)*vel_l(i,j,k,3)
            fin_l(i,j,k,n) = 1.0*wt(n)*
     &      (1.0 + 3.0*cu + 4.5*cu**2 - 1.5*usqr)
          enddo

      enddo
      enddo
      enddo


      write(6,*)rank,imin_l,imax_l,jmin_l,jmax_l,kmin_l,kmax_l


C==============
      DO iter = 1, nits
C==============

C=======right wall outflow condition
        if (imax_l.eq.im) then
          do k = 1, km_l
          do j = 1, jm_l
            do c = 1,nrow
              fin_l(im_l,j,k,col3(c)) = fin_l(im_l-1,j,k,col3(c))
            enddo
          enddo
          enddo
        endif

C=======compute macroscopic variables rho and u
        rho_l = 0.0
        vel_l = 0.0

        do k = 1, km_l
        do j = 1, jm_l
        do i = 1, im_l

          do n = 1,nvec
            rho_l(i,j,k) = rho_l(i,j,k) + fin_l(i,j,k,n)
            vel_l(i,j,k,1) = vel_l(i,j,k,1) + vec(n,1)*fin_l(i,j,k,n)
            vel_l(i,j,k,2) = vel_l(i,j,k,2) + vec(n,2)*fin_l(i,j,k,n)
            vel_l(i,j,k,3) = vel_l(i,j,k,3) + vec(n,3)*fin_l(i,j,k,n)
          enddo

          vel_l(i,j,k,1) = vel_l(i,j,k,1) / rho_l(i,j,k)
          vel_l(i,j,k,2) = vel_l(i,j,k,2) / rho_l(i,j,k)
          vel_l(i,j,k,3) = vel_l(i,j,k,3) / rho_l(i,j,k)

        enddo
        enddo
        enddo


C=======left wall inflow condition
        if (imin_l.eq.1) then
        do k = 1, km_l
        do j = 1, jm_l

          vel_l(1,j,k,1) = vin(1)
          vel_l(1,j,k,2) = vin(2)
          vel_l(1,j,k,3) = vin(3)

          sum2 = 0.0
          sum3 = 0.0

          do c = 1,nrow
            sum2 = sum2 + fin_l(1,j,k,col2(c))
            sum3 = sum3 + fin_l(1,j,k,col3(c))
            rho_l(1,j,k) = (sum2 + 2.0*sum3) / (1.0-vel_l(1,j,k,1))
          enddo

        enddo
        enddo
        endif



C=======compute equilibrium (rho, vel)
        do k = 1,km_l
        do j = 1,jm_l
        do i = 1,im_l

          usqr = vel_l(i,j,k,1)**2
     &         + vel_l(i,j,k,2)**2
     &         + vel_l(i,j,k,3)**2

          do n = 1,nvec
            cu = vec(n,1)*vel_l(i,j,k,1)
     &         + vec(n,2)*vel_l(i,j,k,2)
     &         + vec(n,3)*vel_l(i,j,k,3)

            feq_l(i,j,k,n) = rho_l(i,j,k)*wt(n)*
     &     (1.0 + 3.0*cu + 4.5*cu**2 - 1.5*usqr)
          enddo

        enddo
        enddo
        enddo


C=======calculate populations (at inlet)       
        if (imin_l.eq.1) then
        do j = 1, jm_l
        do k = 1, km_l
          do c = 1,nrow
            fin_l(1,j,k,col1(c)) = feq_l(1,j,k,col1(c))
     &                           + fin_l(1,j,k,col3(6-c))
     &                           - feq_l(1,j,k,col3(6-c))
          enddo
        enddo
        enddo
        endif



C=======collision step
        do k = 1, km_l
        do j = 1, jm_l
        do i = 1, im_l
          do n = 1,nvec
            fout_l(i,j,k,n) = fin_l(i,j,k,n)
     &               - omega*(fin_l(i,j,k,n) - feq_l(i,j,k,n))
          enddo
        enddo
        enddo
        enddo



C=======bounce-back condition
        do k = 1, km_l
        do j = 1, jm_l
        do i = 1, im_l

        ii = imin_l + (i-1)
        jj = jmin_l + (j-1)
        kk = kmin_l + (k-1)

        if (obs(ii,jj,kk).eq.1) then ! obstacle
          do n = 1,nvec
            fout_l(i,j,k,n) = fin_l(i,j,k,16-n)
          enddo
        endif

        enddo
        enddo
        enddo


C=======communicate before streaming step
        CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
        CALL bcs(fout_l)
        CALL MPI_Barrier(MPI_COMM_WORLD, ierr)


C=======streaming step
        do k = 1, km_l
        do j = 1, jm_l
        do i = 1, im_l
        do n = 1, nvec
           nexti = i + vec(n,1)
           nextj = j + vec(n,2)
           nextk = k + vec(n,3)
           fin_l(nexti,nextj,nextk,n) = fout_l(i,j,k,n)
        enddo
        enddo
        enddo
        enddo


C=======streaming at xmin/xmax block boundaries
        do j = 1, jm_l
        do k = 1, km_l
        do n = 1, nvec

          if (imin_l.ne.1) then ! rhs block inlet
            i = 0
            nexti = i + vec(n,1)
            nextj = j + vec(n,2)
            nextk = k + vec(n,3)
            if  (nexti.eq.1) then
              fin_l(nexti,nextj,nextk,n) = fout_l(i,j,k,n)
            endif
          endif

          if (imax_l.ne.im) then ! lhs block outlet
            i = im_l+1
            nexti = i + vec(n,1)
            nextj = j + vec(n,2)
            nextk = k + vec(n,3)
            if  (nexti.eq.im_l) then
              fin_l(nexti,nextj,nextk,n) = fout_l(i,j,k,n)
            endif
          endif

        enddo
        enddo
        enddo





        IF (MOD(iter,nout) == 0) THEN

          CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
          CALL gather_to_zero
          CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
          IF (rank == 0) write(6,*) "iter = ",iter,
     &               "vel = ",vel(im/2,jm/2,km/2,1)

        END IF


C===========
      END DO ! iterations
C===========
 
      CALL MPI_Finalize(ierr)
      IF (rank == 0) write(6,*) "solver finished"
      IF (rank == 0) CALL write_binary_vtk()

      END


