PROGRAM TEST
! Driver code for computing a Delaunay simplex containing a given point.
! Reads a file containing a list of N vertices in R^D and constructs the
! umbrella neighbourhood of some vertex in that file.
! 
! Author: Tyler Chang                                   
! Last Update: August, 2017                               
!
! Usage :$ ./delfan D N V FILENAME
!
! D is an integer specifying the dimension of problem.
! N is an integer specifying the number of lines (vertices) in the file.
! V is an integer specifying the index of vertex to compute the umbrella
!    neighbourhood about.
! FILENAME gives the input file path (32 chars max).
USE DELAUNAYFAN_MOD
IMPLICIT NONE

! Variables and input data.
INTEGER :: D, N, V, ERROR
REAL(KIND=R8), ALLOCATABLE :: PTS(:,:), WEIGHTS(:,:), WORK(:)
REAL(KIND=R8) :: START, FINISH
INTEGER, ALLOCATABLE :: FAN(:,:)
INTEGER :: I, J
CHARACTER(LEN=80) :: ARG

! Read the values of D, N, V, and FILENAME from the command line.
CALL GET_COMMAND_ARGUMENT(1, ARG)
READ(ARG, *) D
CALL GET_COMMAND_ARGUMENT(2, ARG)
READ(ARG, *) N
CALL GET_COMMAND_ARGUMENT(3, ARG)
READ(ARG, *) V
CALL GET_COMMAND_ARGUMENT(4, ARG)

! Allocate memory to store inputs.
ALLOCATE(PTS(D,N), WORK(5*D), STAT=I)
IF (I .NE. 0) PRINT *, 'Warning: Malloc Fail!'

! Read data from input file.
OPEN(1, FILE=TRIM(ARG))
DO I = 1, N
   READ(1, *) PTS(:, I)
END DO
CLOSE(1)

! Get the results and record the computation time.
CALL CPU_TIME(START)
CALL DELAUNAYFAN(D, N, PTS, V, FAN, ERROR, IBUDGET=10000)
CALL CPU_TIME(FINISH)
FINISH = FINISH - START
! Print the timing results.
IF(ERROR .LE. 1) THEN
   PRINT *, 'Fan: '
   DO I = 1, SIZE(FAN,2)
      PRINT *, FAN(:,I)
   END DO
   PRINT *, ''
   IF (ERROR .EQ. 1) PRINT *, 'Vertex ', V, ' is a vertex of the convex hull.'
   PRINT *, ''
ELSE
   PRINT *, 'Error. Code = ', ERROR
END IF
PRINT *, 'Fan of ', SIZE(FAN,2), ' simps built in ', FINISH, ' seconds.'

! Free the heap memory.
DEALLOCATE(PTS, FAN, WORK, STAT=ERROR)
IF (ERROR .NE. 0) PRINT *, 'There was an error freeing memory.'
RETURN
END PROGRAM TEST
