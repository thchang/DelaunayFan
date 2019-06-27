MODULE REAL_PRECISION ! HOMPACK90 module for 64-bit arithmetic.
INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

MODULE DELAUNAYFAN_MOD
! This module contains the subroutine VTDELAUNAYFAN for computing the
! umbrella neighbourhood of a vertex V in R^D from the Delaunay
! triangulation of data points PTS.
USE REAL_PRECISION

PUBLIC

CONTAINS

SUBROUTINE DELAUNAYFAN( D, N, PTS, V, FAN, IERR, EPS, IBUDGET )
! This is a serial implementation of an algorithm for computing the umbrella
! neighborhood of a single vertex in R^D in the Delaunay triangulation.
! The algorithm is fully described and analyzed in
!
! T. H. Chang, L. T. Watson, T. C.H. Lux, S. Raghvendra, B. Li, L. Xu,
! A. R. Butt, K. W. Cameron, and Y. Hong. Computing the umbrella
! neighbourhood of a vertex in the Delaunay triangulation and a single
! Voronoi cell in arbitrary dimension. In Proceedings of IEEE Southeast
! Conference 2018 (SEC 2018). IEEE, St. Petersburg, FL, 2018.
!
! Note that this algorithm is NOT robust for degenerate or near degenerate
! data. Such inputs will result in an error.
!
! On input:
!
! D is the dimension of the space for PTS.
!
! N is the number of data points in PTS.
!
! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
!    coordinates of a single data point in R^D.
!
! V is the index in PTS of the vertex about which to construct the umbrella
!    neighbourhood.
!
! FAN(:,:) is an unallocated ALLOCATABLE array of type INTEGER that
!    will be allocated on output.
!
!
! On output:
!
! PTS has been rescaled and shifted. All the data points in PTS are now
!    contained in the unit hyperball in R^D.
!
! FAN(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
!    in PTS) for the D+1 vertices of each Delaunay simplex in the umbrella
!    neighbourhood of PTS(:,V).
!
! IERR is an integer valued error flag. The error codes are:
!
! 00 : Succesfully computed the entire umbrella neighbourhood.
!
! 10 : The input data sizes don't match or contained illegal values.
! 11 : Too few data points to construct a triangulation (i.e., N < D+1).
! 12 : Two or more points in the data set PTS are too close together with
!      respect to the working precision (EPS), which would result in a
!      numerically degenerate simplex.
! 13 : The input V does not contain a valid index in PTS.
!
! 20 : A budget violation has occurred. All computed simplices were still
!      returned, but the complete neighbourhood was not computed. Note,
!      it is possible that this error was caused by a degeneracy that
!      went undetected, resulting in an infinite loop.
! 21 : The supplied data was in some way degenerate or nearly degenerate,
!      causing the algorithm to fail.
!
! 31 : All the data points in PTS lie in some lower dimensional linear
!      manifold (up to the working precision), and no valid triangulation
!      exists.
!
! 50 : An error occurred while managing the dynamic arrays containing the
!      FaceList and SimplexLists.
! 51 : A memory allocation error occurred while allocating the output
!      array FAN.
!
! 90 : A LAPACK subroutine has reported an illegal value.
! 91 : The LAPACK subroutine DGESVD failed to converge while performing a
!      singular value decomposition.
!
!
! Optional arguments:
!
! EPS contains the working precision for the problem on input. By default,
!    EPS is assigned \sqrt{\mu} where \mu denotes the unit roundoff for the
!    machine. In general, any values that differ by less than EPS are judged
!    as equal, and any weights that are greater than -EPS are judged as
!    nonnegative.  EPS cannot take a value less than the default value of
!    \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied, the default
!    value will be used instead automatically. 
! 
! IBUDGET contains the integer valued budget for performing flips while
!    iterating toward the simplex containing each interpolation point in Q.
!    This prevents DelaunayFan from falling into an infinite loop when
!    supplied with degenerate or near degenerate data.  By default,
!    IBUDGET=50000. However, for extremely high-dimensional problems and
!    pathological data sets, the default value may be insufficient. 
!
!
! Subroutines and functions directly referenced from BLAS are
!      DDOT, DNRM2, DTRSM,
! and from LAPACK are
!      DGELS, DGEQP3, DGESV, DGESVD.
! 
! Primary Author: Tyler H. Chang
! Last Update: June, 2019

USE AFL
IMPLICIT NONE

! Input arguments.
INTEGER, INTENT(IN) :: D, N
REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
INTEGER, INTENT(IN) :: V
! Output arguments.
INTEGER, ALLOCATABLE, INTENT(OUT) :: FAN(:,:)
INTEGER, INTENT(OUT) :: IERR
! Optional arguments.
REAL(KIND=R8), INTENT(IN), OPTIONAL:: EPS
INTEGER, INTENT(IN), OPTIONAL :: IBUDGET

! Local variables.
INTEGER :: IBUDGETL ! Local copy of IBUDGET.
INTEGER :: ITMP ! Temporary value for swapping.
INTEGER :: I, J, K ! Loop iteration variables.
INTEGER :: LWORK ! Size of WORK array (5*D).
REAL(KIND=R8) :: CURRRAD ! Radius of the circumball.
REAL(KIND=R8) :: EPSL ! Local copy of EPS.
REAL(KIND=R8) :: MINRAD ! Smallest radius found.
REAL(KIND=R8) :: SIDE1 ! Side of the hyperplane to flip toward.
REAL(KIND=R8) :: SIDE2 ! Side of the hyperplane for a point.

! Local arrays requiring O(d^2) extra memory.
INTEGER :: IPIV(D) ! Pivot indices.
INTEGER :: SIMP(D+1) ! Current simplex.
REAL(KIND=R8) :: AT(D,D) ! The transpose of A, the linear coefficient matrix.
REAL(KIND=R8) :: B(D) ! The RHS of a linear system.
REAL(KIND=R8) :: CENTER(D) ! The circumcenter of a simplex.
REAL(KIND=R8) :: LQ(D,D) ! Hold LU or QR of AT.
REAL(KIND=R8) :: PLANE(D+1) ! The hyperplane containing a facet.
REAL(KIND=R8) :: TAU(D) ! Householder reflector constants.
REAL(KIND=R8) :: WORK(5*D) ! Work array for DGEQP3 and DORMQR.

! Dynamic arrays of unknown size from AFL and Vector modules.
TYPE(INTVECTOR) :: SL
TYPE(FACELIST) :: FL

! External functions and subroutines.
REAL(KIND=R8), EXTERNAL :: DDOT ! Inner product (BLAS).
REAL(KIND=R8), EXTERNAL :: DNRM2 ! Euclidean norm (BLAS).
EXTERNAL :: DGEQP3 ! Perform a QR factorization with column pivoting (LAPACK).
EXTERNAL :: DGETRF ! Perform a LU factorization with partial pivoting (LAPACK).
EXTERNAL :: DGETRS ! Use the output of DGETRF to solve a linear system (LAPACK).
EXTERNAL :: DORMQR ! Apply householder reflectors to a matrix (LAPACK).
EXTERNAL :: DTRSM ! Perform a triangular solve (BLAS).

! Check for input size errors.
IF ( D < 1 .OR. N < 1 .OR. SIZE(PTS,1) .NE. D &
 .OR. SIZE(PTS,2) .NE. N) THEN
   IERR = 10
   RETURN
ELSE IF (N < D+1) THEN
   IERR = 11
   RETURN
ELSE IF (V > N .OR. V < 1) THEN
   IERR = 13
   RETURN
END IF
! Compute the machine precision.
EPSL = SQRT(EPSILON(1.0_R8))
IF (PRESENT(EPS)) THEN
   IF(EPSL < EPS) THEN
      EPSL = EPS
   END IF
END IF
! Rescale the points to the unit hypersphere.
CALL RESCALE(MINRAD)
! Check for degeneracies in point spacing.
IF (MINRAD < EPSL) THEN
   IERR = 12
   RETURN
END IF
! Check other optional inputs.
IF (PRESENT(IBUDGET)) THEN
   IBUDGETL = IBUDGET
   IF (IBUDGETL < 1) THEN
      IERR = 10
      RETURN
   END IF
ELSE
   IBUDGETL = 50000
END IF

! Initialize the simplex and face lists.
SL = INTVECTOR(DIM=D+1, ISTAT=I)
IF (I .NE. 0) THEN
   IERR = 50
   RETURN
END IF
FL = NEWFACELIST(DIM=D+1, IND=V, ISTAT=I)
IF (I .NE. 0) THEN
   IERR = 50
   RETURN
END IF
! Initialize LWORK
LWORK = 5*D

! Make the first "seed" simplex.
CALL MAKEFIRSTSIMP()
IF(IERR .NE. 0) RETURN
! Add all the faces to the stack.
CALL SELECTSORTSIMP()
CALL FL%PUSHALL(SIMP, IERR)
IF(IERR .NE. 0) RETURN
! Save the first simplex.
SL%PUSH(SIMP, IERR)
IF(IERR .NE. 0) THEN
   IERR = 50
   RETURN
END IF

! Loop for filling all simplices containing the vertex with index V.
INNER : DO K = 1, IBUDGETL
!! Debugging statements. Uncomment for step-by-step print out of progress.
!PRINT *, 'Iteration: ', K
!PRINT *, 'Simplices:'
!DO I=1, SL%LENGTH()
!PRINT *, SL%ITEM(I)
!END DO
!PRINT *, 'Faces:'
!DO I=1, FL%LENGTH()
!PRINT *, FL%ITEM(I)
!END DO
   ! Get the top face off the stack, but don't delete.
   CALL FL%TOP(SIMP, IERR)
   IF(IERR .NE. 0) THEN
      IERR = 50
      RETURN
   END IF
   IF (SIMP(1) .EQ. 0) EXIT INNER ! Exit condition.
   ! Reconstruct the linear system for the new face.
   DO I = 1, D-1
      AT(:,I) = PTS(:,SIMP(I+1)) - PTS(:,SIMP(1))
      B(I) = DDOT(D, AT(:,I), 1, AT(:,I), 1) / 2.0_R8
   END DO
   ! Compute the next simplex.
   CALL MAKESIMPLEX()
   IF (IERR .NE. 0) RETURN
   ! Check for the extrapolation condition.
   IF (SIMP(D+1) .EQ. 0) THEN
      CALL FL%POP(SIMP, IERR)
      IF (IERR .NE. 0) RETURN
      CYCLE INNER
   END IF
   ! Sort the simplex indices and push all its faces to the stack.
   CALL BUBBLESORTSIMP()
   CALL FL%PUSHSIMP(SIMP, IERR)
   IF(IERR .NE. 0) RETURN
   ! Save the current simplex to the stack.
   CALL SL%PUSH(SIMP, IERR)
   IF(IERR .NE. 0) THEN
      IERR = 50
      RETURN
   END IF
END DO INNER

! Check for budget violation.
IF (K > IBUDGETL) THEN 
   IERR = 20
END IF

! Copy the fan into the output array.
ALLOCATE(FAN(D+1,SL%LENGTH()), STAT=I)
IF(I .NE. 0) THEN
   IERR = 51
   RETURN
END IF
FAN(:,:) = SL%DATA()
CALL SL%FREE()
RETURN

CONTAINS ! Internal subroutines and functions.

SUBROUTINE MAKEFIRSTSIMP()
! Iteratively construct the first simplex by choosing points that
! minimize the radius of the smallest circumball. Let P_1, P_2, ..., P_K
! denote the current set of vertices for the simplex. Let P* denote the
! candidate vertex to be added to the simplex. Let CENTER denote the
! circumcenter of the simplex.  Then
!
! X = CENTER - P_1
!
! is given by the minimum norm solution to the underdetermined linear system 
!
! A X = B, where
!
! A^T = [ P_2 - P_1, P_3 - P_1, ..., P_K - P_1, P* - P_1 ] and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{K.},A_{K.}>/2 ]^T.
!
! Then the radius of the smallest circumsphere is CURRRAD = \| X \|,
! and the next vertex is given by P_{K+1} = argmin_{P*} CURRRAD, where P*
! ranges over points in PTS that are not already a vertex of the simplex.
!
! On output, this subroutine fully populates the matrix A^T and vector B,
! and fills SIMP(:) with the indices of a valid Delaunay simplex.

! Add the first point.
SIMP(:) = 0
SIMP(1) = V
! Find the second point, i.e., the closest point to PTS(:,SIMP(1)).
MINRAD = HUGE(0.0_R8)
DO I = 1, N
   ! Skip repeated vertices.
   IF (I .EQ. SIMP(1)) CYCLE
   ! Check the diameter of the resulting circumsphere.
   CURRRAD = DNRM2(D, PTS(:,I)-PTS(:,SIMP(1)), 1)
   IF (CURRRAD < MINRAD) THEN; MINRAD = CURRRAD; SIMP(2) = I; END IF
END DO
! Set up the first row of the system A X = B.
AT(:,1) = PTS(:,SIMP(2)) - PTS(:,SIMP(1))
B(1) = DDOT(D, AT(:,1), 1, AT(:,1), 1) / 2.0_R8
! Loop to collect the remaining D-1 vertices for the first simplex.
DO I = 2, D
   MINRAD = HUGE(0.0_R8) ! Re-initialize the radius for each iteration.
   ! Check each point P* in PTS.
   DO J = 1, N
      ! Check that this point is not already in the simplex.
      IF (ANY(SIMP(:) .EQ. J)) CYCLE
      ! Add P* to linear system, and compute the minimum norm solution.
      AT(:,I) = PTS(:,J) - PTS(:,SIMP(1))
      B(I) = DDOT(D, AT(:,I), 1, AT(:,I), 1) / 2.0_R8
      LQ(:,1:I) = AT(:,1:I)
      ! Compute A^T P = Q R.
      CALL DGEQP3(D, I, LQ, D, IPIV, TAU, WORK, LWORK, IERR)
      IF(IERR < 0) THEN ! LAPACK illegal input error.
         IERR = 80; RETURN
      ELSE IF (ABS(LQ(I,I)) < EPSL) THEN ! A is rank-deficient.
         CYCLE ! If rank-deficient, skip this point.
      END IF
      ! Set CENTER = P^T B.
      FORALL (ITMP = 1:I) CENTER(ITMP) = B(IPIV(ITMP))
      ! Get the radius using R^T Q^T X = P^T B.
      CALL DTRSM('L', 'U', 'T', 'N', I, 1, 1.0_R8, LQ, D, CENTER, D)
      CENTER(I+1:D) = 0.0_R8
      CALL DORMQR('L', 'N', D, 1, I, LQ, D, TAU, CENTER, D, WORK, LWORK, &
         IERR)
      IF(IERR < 0) THEN ! LAPACK illegal input error.
         IERR = 83; RETURN
      END IF
      ! Calculate the radius and compare it to the current minimum.
      CURRRAD = DNRM2(D, CENTER, 1)
      IF (CURRRAD < MINRAD) THEN; MINRAD = CURRRAD; SIMP(I+1) = J; END IF
   END DO
   ! Check that a point was found. If not, then all the points must lie in a
   ! lower dimensional linear manifold (error case).
   IF (SIMP(I+1) .EQ. 0) THEN; IERR = 31; RETURN; END IF
   ! If all operations were successful, add the best P* to the linear system.
   AT(:,I) = PTS(:,SIMP(I+1)) - PTS(:,SIMP(1))
   B(I) = DDOT(D, AT(:,I), 1, AT(:,I), 1) / 2.0_R8 
END DO
IERR = 0 ! Set error flag to 'success' for a normal return.
RETURN
END SUBROUTINE MAKEFIRSTSIMP

SUBROUTINE MAKESIMPLEX()
! Given a Delaunay facet F and a specified side of its containing hyperplane,
! complete the simplex by adding a point from PTS on the opposite `side' of
! F. Assume SIMP(1:D) contains the vertex indices of F (corresponding to
! data points P_1, P_2, ..., P_D in PTS), and assume the matrix A(1:D-1,:)^T
! and vector B(1:D-1) are filled appropriately (similarly as in
! MAKEFIRSTSIMP()). Then for any P* (not in the hyperplane containing
! F) in PTS, let CENTER denote the circumcenter of the simplex with vertices
! P_1, P_2, ..., P_D, P*. Then
!
! X = CENTER - P_1
!
! is given by the solution to the nonsingular linear system
!
! A X = B where
!
! A^T = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1, P* - P_1 ] and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{D.},A_{D.}>/2 ]^T.
!
! Then CENTER = X + P_1 and RADIUS = \| X \|.  P_{D+1} will be given by the 
! candidate P* that satisfies both of the following:
!
! 1) Let PLANE denote the hyperplane containing F. Then P_{D+1} and
! PTS(SIMP(D+1)) must be on opposite sides of PLANE.
!
! 2) The circumball about CENTER must not contain any points in PTS in its
! interior (Delaunay property).
! 
! The above are necessary and sufficient conditions for flipping the
! Delaunay simplex, given that F is indeed a Delaunay facet.
!
! On input, SIMP(1:D) should contain the vertex indices (column indices
! from PTS) of the facet F and SIMP(D+1) must contain the index of a piont
! on the current side of F (to flip away from). Upon output, SIMP(:) will
! contain the vertex indices of a Delaunay simplex on the opposite side of
! PLANE.  Also, the matrix A^T and vector B will be updated accordingly. If
! SIMP(D+1)=0, then there were no points in PTS on the appropriate side of
! F, meaning that F is a facet of the convex hull of PTS.

! Compute the hyperplane PLANE.
CALL MAKEPLANE()
IF(IERR .NE. 0) RETURN ! Check for errors.
! Compute the sign for the side of the viable points.
SIDE1 =  -1.0_R8 * DDOT(D,PLANE(1:D),1,PTS(:,SIMP(D+1)),1) + PLANE(D+1)
SIDE1 = SIGN(1.0_R8,SIDE1)
! Initialize the center, radius, and simplex.
SIMP(D+1) = 0
CENTER(:) = 0.0_R8
MINRAD = HUGE(0.0_R8)
! Loop through all points P* in PTS.
DO I = 1, N
   ! Check that P* is inside the current ball.
   IF (DNRM2(D, PTS(:,I) - CENTER(:), 1) > MINRAD) CYCLE ! If not, skip.
   ! Check that P* is on the appropriate halfspace.
   SIDE2 = DDOT(D,PLANE(1:D),1,PTS(:,I),1) - PLANE(D+1)
   IF (SIDE1*SIDE2 < EPSL .OR. ANY(SIMP(:) .EQ. I)) CYCLE ! If not, skip.
   ! Add P* to the linear system, and solve to get shifted CENTER.
   AT(:,D) = PTS(:,I) - PTS(:,SIMP(1))
   B(D) = DDOT(D, AT(:,D), 1, AT(:,D), 1) / 2.0_R8
   LQ = AT
   CENTER = B
   ! Compute A^T = LU
   CALL DGETRF(D, D, LQ, D, IPIV, IERR)
   IF (IERR < 0) THEN ! LAPACK illegal input.
      IERR = 81; RETURN
   ELSE IF (IERR > 0) THEN ! Rank-deficiency detected.
      IERR = 61; RETURN
   END IF
   ! Use A^T = LU to solve A X = B, where X = CENTER - P_1.
   CALL DGETRS('T', D, 1, LQ, D, IPIV, CENTER, D, IERR)
   IF (IERR < 0) THEN ! LAPACK illegal input.
      IERR = 82; RETURN
   END IF
   ! Update the new radius, center, and simplex.
   MINRAD = DNRM2(D, CENTER, 1)
   CENTER(:) = CENTER(:) + PTS(:,SIMP(1))
   SIMP(D+1) = I
END DO
IERR = 0 ! Reset the error flag to 'success' code.
! Check for extrapolation condition.
IF(SIMP(D+1) .EQ. 0) RETURN
! Add new point to the linear system.
AT(:,D) = PTS(:,SIMP(D+1)) - PTS(:,SIMP(1))
B(D) = DDOT(D, AT(:,D), 1, AT(:,D), 1) / 2.0_R8
RETURN
END SUBROUTINE MAKESIMPLEX

SUBROUTINE MAKEPLANE()
! Construct a hyperplane c^T x = \alpha containing the first D vertices indexed
! in SIMP(:). The plane is determined by its normal vector c and \alpha.
! Let P_1, P_2, ..., P_D be the vertices indexed in SIMP(1:D). A normal
! vector is any nonzero vector in ker A, where the matrix
!
! A^T = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1 ].
!
! Since rank A = D-1, dim ker A = 1, and ker A can be found from a QR
! factorization of A^T:  A^T P = QR, where P permutes the columns of A^T.
! Then the last column of Q is orthogonal to the range of A^T, and in ker A.
!
! Upon output, PLANE(1:D) contains the normal vector c and PLANE(D+1)
! contains \alpha defining the plane.

IF (D > 1) THEN ! Check that D-1 > 0, otherwise the plane is trivial.
   ! Compute the QR factorization.
   IPIV=0
   LQ = AT
   CALL DGEQP3(D, D-1, LQ, D, IPIV, TAU, WORK, LWORK, IERR)
   IF(IERR < 0) THEN ! LAPACK illegal input error.
      IERR = 80; RETURN
   END IF
   ! The nullspace is given by the last column of Q.
   PLANE(1:D-1) = 0.0_R8
   PLANE(D) = 1.0_R8
   CALL DORMQR('L', 'N', D, 1, D-1, LQ, D, TAU, PLANE, D, WORK, &
    LWORK, IERR)
   IF(IERR < 0) THEN ! LAPACK illegal input error.
      IERR = 83; RETURN
   END IF
   ! Calculate the constant \alpha defining the plane.
   PLANE(D+1) = DDOT(D,PLANE(1:D),1,PTS(:,SIMP(1)),1)
ELSE ! Special case where D=1.
   PLANE(1) = 1.0_R8
   PLANE(2) = PTS(1,SIMP(1))
END IF
RETURN
END SUBROUTINE MAKEPLANE

SUBROUTINE SELECTSORTSIMP()
! Sort the indices in SIMP(:), but keep SIMP(1) fixed. O(D^2) complexity.
DO I = 2, D+1
   J = MINLOC(SIMP(I:D+1),DIM=1) + I - 1
   ITMP = SIMP(J)
   SIMP(J) = SIMP(I)
   SIMP(I) = ITMP
END DO
RETURN
END SUBROUTINE SELECTSORTSIMP

SUBROUTINE BUBBLESORTSIMP()
! Sort the indices in SIMP(:), given that only SIMP(D+1) is out of place and
! SIMP(1) is fixed. O(D) complexity.
DO I = D+1, 3, -1
   IF (SIMP(I) < SIMP(I-1)) THEN
      ITMP = SIMP(I)
      SIMP(I) = SIMP(I-1)
      SIMP(I-1) = ITMP
   ELSE
      RETURN
   END IF
END DO
RETURN
END SUBROUTINE BUBBLESORTSIMP

SUBROUTINE RESCALE(MINDIST)
! Rescale and transform data to be centered at the origin with unit
! radius. This subroutine has O(n^2) complexity.
!
! On output, PTS has been rescaled and shifted. All the data
! points in PTS are centered with unit radius.
!
! MINDIST is a real number containing the (scaled) minimum distance
!    between any two data points in PTS.

! Output arguments.
REAL(KIND=R8), INTENT(OUT) :: MINDIST

! Local variables.
REAL(KIND=R8) :: PTS_CENTER(D) ! The center of the data points PTS.
REAL(KIND=R8) :: DISTANCE ! The current distance.
REAL(KIND=R8) :: SCALE ! The scale factor.

! Initialize local values.
MINDIST = HUGE(MINDIST)

! Compute barycenter of all data points.
PTS_CENTER(:) = SUM(PTS(:,:), DIM=2)/REAL(N, KIND=R8)
! Compute the minimum and maximum distances.
DO I = 1, N ! Cycle through all pairs of points.
   DO J = I + 1, N
      DISTANCE = DNRM2(D, PTS(:,I) - PTS(:,J), 1) ! Compute the distance.
      IF (DISTANCE < MINDIST) THEN ! Compare to the current minimum distance.
         MINDIST = DISTANCE
      END IF
   END DO
END DO
! Center the points.
FORALL (I = 1:N) PTS(:,I) = PTS(:,I) - PTS_CENTER(:)
! Compute the scale factor (for unit radius).
DO I = 1, N ! Cycle through all points again.
   DISTANCE = DNRM2(D, PTS(:,I), 1) ! Compute the distance from the center.
   IF (DISTANCE > SCALE) THEN ! Compare to the current radius.
      SCALE = DISTANCE
   END IF
END DO
! Scale the points to unit radius.
PTS = PTS / SCALE
MINDIST = MINDIST / SCALE
RETURN
END SUBROUTINE RESCALE

END SUBROUTINE DELAUNAYFAN

END MODULE DELAUNAYFAN_MOD
