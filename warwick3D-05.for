
C     https://warwick.ac.uk/research/rtp/sc/rse/training/advancedmpi/02_case_study_1.pdf
C     Chris Brady & Heather Ratcliffe (warwick university)
C     compile: mpif77 -o warwick warwick.for
C     or:      mpif90 -o warwick warwick.for
C     run:     mpirun -np 16 ./warwick

      MODULE display

      USE mpi

      IMPLICIT NONE

      INTEGER, PARAMETER :: nx = 200, ny = 200, nz = 200
      INTEGER, PARAMETER :: tag = 100
      REAL, DIMENSION(0:nx+1, 0:ny+1, 0:nz+1) :: values
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: values_local, temp_local

      INTEGER :: nx_local, ny_local, nz_local
      INTEGER :: x_cell_min_local, x_cell_max_local
      INTEGER :: y_cell_min_local, y_cell_max_local
      INTEGER :: z_cell_min_local, z_cell_max_local

      INTEGER :: nproc, rank, cart_comm
      INTEGER, DIMENSION(3) :: nprocs, coordinates
      INTEGER :: x_min_rank, x_max_rank
      INTEGER :: y_min_rank, y_max_rank
      INTEGER :: z_min_rank, z_max_rank

      CONTAINS




      SUBROUTINE gather_to_zero

      INTEGER :: ierr
      REAL, DIMENSION(0:nx+1, 0:ny+1, 0:nz+1) :: red

      values = 0.0
      red = 0.0

      values(x_cell_min_local:x_cell_max_local,
     &       y_cell_min_local:y_cell_max_local,
     &       z_cell_min_local:z_cell_max_local) = 
     &       values_local(1:nx_local, 1:ny_local, 1:nz_local)

      CALL MPI_Reduce(values,red,(nx+2)*(ny+2)*(nz+2),
     &     MPI_REAL,MPI_SUM,0,cart_comm,ierr)

      values = red

      END SUBROUTINE gather_to_zero




      SUBROUTINE bcs(array)

C     REAL, DIMENSION(0:,0:,0:), INTENT(INOUT) :: array
      REAL, DIMENSION(0:nx_local+1, 0:ny_local+1, 0:nz_local+1),
     &  INTENT(INOUT) :: array
      INTEGER :: ierr

C     write(6,*) "got here 11"
C     Send left most strip of cells left and receive into right guard cells
      CALL MPI_Sendrecv(
     & array(1,         1:ny_local,1:nz_local),
     & ny_local*nz_local, MPI_REAL,x_min_rank,tag, 
     & array(nx_local+1,1:ny_local,1:nz_local),
     & ny_local*nz_local, MPI_REAL,x_max_rank,tag,
     & cart_comm, MPI_STATUS_IGNORE, ierr)

C     write(6,*) "got here 22"
C     Send right most strip of cells right and receive into left guard cells
      CALL MPI_Sendrecv(
     & array(nx_local,1:ny_local,1:nz_local),
     & ny_local*nz_local, MPI_REAL,x_max_rank,tag,
     & array(0,       1:ny_local,1:nz_local),
     & ny_local*nz_local, MPI_REAL,x_min_rank,tag, 
     & cart_comm, MPI_STATUS_IGNORE, ierr)

C     write(6,*) "got here 33"
C     Send lower most strip of cells down and receive into upper guard cells
      CALL MPI_Sendrecv(
     & array(1:nx_local,         1,1:nz_local),
     & nx_local*nz_local, MPI_REAL,y_min_rank,tag, 
     & array(1:nx_local,ny_local+1,1:nz_local),
     & nx_local*nz_local, MPI_REAL,y_max_rank,tag, 
     & cart_comm, MPI_STATUS_IGNORE, ierr)

C     Send upper most strip of cells up and receive into lower guard cells
      CALL MPI_Sendrecv(
     & array(1:nx_local,  ny_local,1:nz_local),
     & nx_local*nz_local, MPI_REAL,y_max_rank,tag, 
     & array(1:nx_local,         0,1:nz_local),
     & nx_local*nz_local, MPI_REAL,y_min_rank,tag, 
     & cart_comm, MPI_STATUS_IGNORE, ierr)

C     Send back most strip of cells back and receive into front guard cells
      CALL MPI_Sendrecv(
     & array(1:nx_local,1:ny_local,         1),
     & nx_local*ny_local, MPI_REAL,z_min_rank,tag, 
     & array(1:nx_local,1:ny_local,nz_local+1),
     & nx_local*ny_local, MPI_REAL,z_max_rank,tag, 
     & cart_comm, MPI_STATUS_IGNORE, ierr)

C     Send front most strip of cells forward and receive into rear guard cells
      CALL MPI_Sendrecv(
     & array(1:nx_local,1:ny_local, nz_local),
     & nx_local*ny_local, MPI_REAL,z_max_rank,tag, 
     & array(1:nx_local,1:ny_local,        0),
     & nx_local*ny_local, MPI_REAL,z_min_rank,tag, 
     & cart_comm, MPI_STATUS_IGNORE, ierr)

      END SUBROUTINE bcs




      SUBROUTINE setup_mpi


      LOGICAL, DIMENSION(3) :: periods
      INTEGER :: ierr 

      periods = .FALSE.

      CALL MPI_Init(ierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      CALL MPI_Dims_create(nproc, 3, nprocs, ierr)

      IF (rank == 0) THEN
        PRINT *,'Processor decomposition is ', nprocs
      ENDIF

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

      nx_local = nx / nprocs(1)
      ny_local = ny / nprocs(2)
      nz_local = nz / nprocs(3)

      CALL MPI_Cart_create(MPI_COMM_WORLD,3,nprocs,periods,.TRUE.,
     &     cart_comm,ierr)

      CALL MPI_Comm_rank(cart_comm, rank, ierr)
      CALL MPI_Cart_shift(cart_comm, 0, 1, x_min_rank, x_max_rank, ierr)
      CALL MPI_Cart_shift(cart_comm, 1, 1, y_min_rank, y_max_rank, ierr)
      CALL MPI_Cart_shift(cart_comm, 2, 1, z_min_rank, z_max_rank, ierr)
      CALL MPI_Cart_coords(cart_comm, rank, 3, coordinates, ierr)

      x_cell_min_local = nx_local *  coordinates(1) + 1
      x_cell_max_local = nx_local * (coordinates(1) + 1)
      y_cell_min_local = ny_local *  coordinates(2) + 1
      y_cell_max_local = ny_local * (coordinates(2) + 1)
      z_cell_min_local = nz_local *  coordinates(3) + 1
      z_cell_max_local = nz_local * (coordinates(3) + 1)


      END SUBROUTINE setup_mpi




      SUBROUTINE write_vtk
      INTEGER i,j,k,im,jm,km
      CHARACTER(20) outfile

      WRITE(6,*) 'write_vtk()...'

!-----filename for vtk output
      outfile = 'warwick.vtk'
      OPEN(UNIT=20, FILE=outfile)

      im = nx
      jm = ny
      km = nz

!-----header
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

!-------scalar temperature
      WRITE(6,*)  '... write temperature'
      WRITE(20,10)'SCALARS temperature float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((( values(i,j,k), i=1,im), j=1,jm), k=1,km)

      WRITE(6,*) '... write FINISHED'
      CLOSE(20)

   10 FORMAT(A)
   20 FORMAT(A,3I5)
   30 FORMAT(A,I5,A)
   40 FORMAT(A,I9)

      END SUBROUTINE write_vtk

!-----------------------------------------------------------------------

      END MODULE display




      PROGRAM serial

      USE display

      IMPLICIT NONE

      INTEGER :: ix, iy, iz, icycle, ierr

      CALL setup_mpi

      ALLOCATE(values_local(0:nx_local+1, 0:ny_local+1, 0:nz_local+1))
      ALLOCATE(  temp_local(0:nx_local+1, 0:ny_local+1, 0:nz_local+1))

C     write(6,*) "got here 1"

      values_local = 0.0

      values_local(0,:,:) = 10.0 ! xmin
      values_local(:,0,:) = 10.0 ! ymin
      values_local(:,:,0) = 10.0 ! zmin

      values_local(nx_local+1,:,:) = 10.0 ! xmax
      values_local(:,ny_local+1,:) = 10.0 ! ymax
      values_local(:,:,nz_local+1) = 10.0 ! zmax

      values = 0.0

C     write(6,*) "got here 2"

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      CALL gather_to_zero
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

C     write(6,*) "got here 3"

      DO icycle = 1, 1000

      DO iz = 1, nz_local
       DO iy = 1, ny_local
        DO ix = 1, nx_local
          temp_local(ix,iy,iz) =
     &       (values_local(ix-1,iy,iz) + values_local(ix+1,iy,iz)
     &      + values_local(ix,iy-1,iz) + values_local(ix,iy+1,iz)
     &      + values_local(ix,iy,iz-1) + values_local(ix,iy,iz+1))/6.0
        END DO
       END DO
      END DO

C     write(6,*) "got here 4"

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

      values_local(1:nx_local,1:ny_local,1:nz_local) = 
     &  temp_local(1:nx_local,1:ny_local,1:nz_local)

C     write(6,*) "got here 5"

      CALL bcs(values_local)

C     write(6,*) "got here 6"

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

      IF (MOD(icycle,100) == 0) THEN
        CALL gather_to_zero
        IF (rank == 0) write(6,*) "icycle = ",icycle
      END IF

C     write(6,*) "got here 7"

      END DO ! iterations
 
      CALL MPI_Finalize(ierr)

      IF (rank == 0) write(6,*) "solver finished"

      IF (rank == 0) CALL write_vtk()

      END PROGRAM serial


