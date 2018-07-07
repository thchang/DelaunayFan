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

SUBROUTINE DELAUNAYFAN( D, N, PTS, V, FAN, IERR, EPS_OPT, IBUDGET_OPT )
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
! FAN(:,:) is an unallocated ALLOCATABLE array of type REAL(KIND=R8) that
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
! EPS_OPT contains the working precision for the problem on input. By default,
!    EPS_OPT is assigned \sqrt{\mu} where \mu denotes the unit roundoff for the
!    machine. In general, any values that differ by less than EPS_OPT are judged
!    as equal, and any weights that are greater than -EPS_OPT are judged as
!    nonnegative.  EPS_OPT cannot take a value less than the default value of
!    \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied, the default
!    value will be used instead automatically. 
! 
! IBUDGET_OPT contains the integer valued budget for performing flips while
!    iterating toward the simplex containing each interpolation point in Q.
!    This prevents DelaunayFan from falling into an infinite loop when
!    supplied with degenerate or near degenerate data.  By default,
!    IBUDGET_OPT=50000. However, for extremely high-dimensional problems and
!    pathological data sets, the default value may be insufficient. 
!
!
! Subroutines and functions directly referenced from BLAS are
!      DDOT, DNRM2, DTRSM,
! and from LAPACK are
!      DGELS, DGEQP3, DGESV, DGESVD.
! 
! Primary Author: Tyler H. Chang
! Last Update: June, 2018

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
REAL(KIND=R8), INTENT(IN), OPTIONAL:: EPS_OPT
INTEGER, INTENT(IN), OPTIONAL :: IBUDGET_OPT

! Local variables.
INTEGER :: IBUDGET ! Local copy of IBUDGET_OPT.
INTEGER :: ITMP ! Temporary value for swapping.
INTEGER :: I, J, K ! Loop iteration variables.
INTEGER :: LWORK ! Size of WORK array (5*D).
REAL(KIND=R8) :: CURRRAD ! Radius of the circumball.
REAL(KIND=R8) :: EPS ! Local copy of EPS_OPT.
REAL(KIND=R8) :: MINRAD ! Smallest radius found.
REAL(KIND=R8) :: SIDE1 ! Side of the hyperplane to flip toward.
REAL(KIND=R8) :: SIDE2 ! Side of the hyperplane for a point.

! Local arrays requiring O(d^2) extra memory.
INTEGER :: PIV(D)
INTEGER :: SIMP(D+1)
REAL(KIND=R8) :: A(D,D)
REAL(KIND=R8) :: B(D)
REAL(KIND=R8) :: CENTER(D)
REAL(KIND=R8) :: LQ(D,D)
REAL(KIND=R8) :: PLANE(D+1)
REAL(KIND=R8) :: WORK(5*D)

! Dynamic arrays of unknown size from AFL and Vector modules.
TYPE(INTVECTOR) :: SL
TYPE(FACELIST) :: FL

! External functions and subroutines.
REAL(KIND=R8), EXTERNAL :: DDOT ! BLAS inner product.
REAL(KIND=R8), EXTERNAL :: DNRM2 ! BLAS Euclidean norm.
EXTERNAL :: DTRSM ! BLAS triangular solve.
EXTERNAL :: DGELS ! LAPACK linear solve.
EXTERNAL :: DGEQP3 ! LAPACK QR factorization.
EXTERNAL :: DGESV ! LAPACK solve least squares via QR.
EXTERNAL :: DGESVD ! LAPACK singular value decomposition.

! Check for input size errors.
IF ( D < 1 .OR. N < 1 .OR. SIZE(PTS,1) .NE. D &
 .OR. SIZE(PTS,2) .NE. N) THEN
   IERR = 10
   RETURN
ELSE IF (N < D+1) THEN
   IERR = 11
   RETURN
END IF
! Compute the machine precision.
EPS = SQRT(EPSILON(1.0_R8))
IF (PRESENT(EPS_OPT)) THEN
   IF(EPS < EPS_OPT .AND. EPS_OPT < 0.05_R8) THEN
      EPS = EPS_OPT
   END IF
END IF
! Rescale the points to the unit hypersphere.
CALL RESCALE(MINRAD)
! Check for degeneracies in point spacing.
IF (MINRAD < EPS) THEN
   IERR = 12
   RETURN
END IF
! Check other optional inputs.
IF (PRESENT(IBUDGET_OPT)) THEN
   IBUDGET = IBUDGET_OPT
   IF (IBUDGET < 1) THEN
      IERR = 10
      RETURN
   END IF
ELSE
   IBUDGET = 50000
END IF

! Initialize the simplex and face lists.
CALL NEWINTVECTOR(SL, I, DIM=D+1)
IF (I .NE. 0) THEN
   IERR = 50
   RETURN
END IF
FL = NEWFACELIST(V, D, I)
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
CALL PUSHALL(FL, SIMP, IERR)
IF(IERR .NE. 0) RETURN
! Save the first simplex.
CALL INTVECTORPUSH(SL, SIMP, IERR)
IF(IERR .NE. 0) THEN
   IERR = 50
   RETURN
END IF

! Loop for filling all simplices containing the vertex with index V.
INNER : DO K = 1, IBUDGET
!! Debugging statements. Uncomment for step-by-step print out of progress.
!PRINT *, 'Iteration: ', k
!PRINT *, 'Simplices:'
!DO i=1, SL%length
!print *, SL%dat(:,i)
!END DO
!PRINT *, 'Faces:'
!DO i=1, FL%Table%length
!print *, FL%TABLE%dat(1:d,i)
!END DO
   ! Get the top face off the stack.
   CALL POPFACE(FL, SIMP, IERR)
   IF(IERR .NE. 0) THEN
      IERR = 50
      RETURN
   END IF
   IF (SIMP(1) .EQ. 0) EXIT INNER ! Exit condition.
   ! Reconstruct the linear system for the new face.
   DO I=1,D-1
      A(I,:) = PTS(:,SIMP(I+1)) - PTS(:,SIMP(1))
      B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8
   END DO
   ! Compute the next simplex.
   CALL MAKESIMPLEX()
   IF (IERR .NE. 0) RETURN
   ! Check for the extrapolation condition.
   IF (SIMP(D+1) .EQ. 0) THEN
      CALL INTVECTORPOP(FL%TABLE, SIMP, IERR)
      IF (IERR .NE. 0) RETURN
      CYCLE INNER
   END IF
   ! Sort the simplex indices and push all its faces to the stack.
   CALL BUBBLESORTSIMP()
   CALL PUSHALL(FL, SIMP, IERR)
   IF(IERR .NE. 0) RETURN
   ! Save the current simplex to the stack.
   CALL INTVECTORPUSH(SL, SIMP, IERR)
   IF(IERR .NE. 0) THEN
      IERR = 50
      RETURN
   END IF
END DO INNER

! Check for budget violation.
IF (K > IBUDGET) THEN 
   IERR = 20
END IF

! Copy the fan into the output array.
ALLOCATE(FAN(D+1,SL%LENGTH), STAT=I)
IF(I .NE. 0) THEN
   IERR = 51
   RETURN
END IF
FAN(:,:) = SL%DAT(1:D+1,1:SL%LENGTH)
CALL INTVECTORFREE(SL, I)
RETURN

CONTAINS ! Internal subroutines and functions.

SUBROUTINE MAKEFIRSTSIMP()
! Iteratively construct an initial simplex starting with the PTS(:,V) by
! choosing points that minimize the radius of the smallest circumball.
! Let (P_1, P_2, ..., P_K) denote the list of vertices for the simplex
! after K iterations. Let P* denote the candidate vertex to be added
! to SIMP in iteration K+1. Let CENTER denote the circumcenter of the
! resulting simplex.  Then
!
! X = CENTER - P_1
!
! is given by the minimum norm solution to the underdetermined linear system 
!
! AX = B, where
!
! A = [ P_2 - P_1, P_3 - P_1, ..., P_K - P_1, P* - P_1 ]^T and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{K.},A_{K.}>/2 ]^T.
!
! Then the radius of the smallest circumsphere is CURRRAD = \| X \|,
! and the next vertex is given by P_{K+1} = argmin_{P*} CURRRAD, where P*
! ranges over points in PTS that are not already a vertex of the simplex.
!
! On output, this subroutine fully populates the matrix A and vector B, and
! fills SIMP(:) with the indices of a valid Delaunay simplex containing V as
! a vertex. This subroutine returns an error if a degenerate case is
! detected.

! Dummy variables U and VT for storing U and V^T when computing SVD.
! Values are never initialized or used.
REAL(KIND=R8), ALLOCATABLE :: U(:), VT(:)
! Add the first point.
SIMP(:) = 0
IF (V > 0 .AND. V .LE. N) THEN
   SIMP(1) = V
ELSE
   IERR = 13
   RETURN
END IF
! The second point is the closest point to the first.
MINRAD = HUGE(MINRAD)
DO I = 1, N
   ! Avoid repeats.
   IF (I == SIMP(1)) CYCLE
   ! Check for a new minimum.
   CURRRAD = DNRM2(D, PTS(:,I)-PTS(:,SIMP(1)), 1)
   IF (CURRRAD < MINRAD) THEN
      MINRAD = CURRRAD
      SIMP(2) = I
   END IF
END DO
! Set up the first row of the least squares system.
A(1,:) = PTS(:,SIMP(2)) - PTS(:,SIMP(1))
B(1) = DDOT(D, A(1,:), 1, A(1,:), 1) / 2.0_R8
! Loop over the remaining D-1 vertices in the first simplex.
DO I = 2, D
   ! Re-initialize the radius for each iteration.
   MINRAD = HUGE(MINRAD)
   ! Find the next point to add.
   DO J = 1, N
      ! Check that the point is not already in the face.
      IF (ANY(SIMP(:) == J)) CYCLE
      ! Add the current point to the least squares system.
      A(I,:) = PTS(:,J) - PTS(:,SIMP(1))
      B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8
      ! Solve least squares problem using QR/LQ factorization.
      LQ(1:I,:) = A(1:I,:)
      CENTER(1:I) = B(1:I)
      CALL DGELS('N', I, D, 1, LQ, D, CENTER, D, WORK, LWORK, IERR)
      IF (IERR < 0) THEN ! Check for errors.
         IERR = 90
         RETURN
      ELSE IF (IERR > 0) THEN ! Indicates rank deficiency.
         CYCLE
      END IF
      ! Calculate the radius of the circumball.
      CURRRAD = DNRM2(D, CENTER, 1)
      ! Check for a new minimum radius.
      IF (CURRRAD < MINRAD) THEN
         MINRAD = CURRRAD
         SIMP(I+1) = J
      END IF
   END DO
   ! Check that a point was found.
   IF (SIMP(I+1) .EQ. 0) THEN
      EXIT
   END IF
   ! Add the next point to the least squares system.
   A(I,:) = PTS(:,SIMP(I+1)) - PTS(:,SIMP(1))
   B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8 
END DO
! If no rank deficiency was detected, double-check.
IF (I > D) THEN
   ! Compute the singular value decomposition of A.
   LQ = A
   CALL DGESVD('N','N',D,D,LQ,D,PLANE(1:D),U,1,VT,1,WORK,LWORK,IERR)
   IF (IERR < 0) THEN ! Check for errors.
      IERR = 90
      RETURN
   ELSE IF (IERR > 0) THEN ! Failure to converge.
      IERR = 91
      RETURN
   END IF
   ! Check for rank deficiency up to working precision.
   IF (PLANE(D)/PLANE(1) < EPS) THEN
      IERR = 21
      RETURN
   END IF
ELSE ! Otherwise, rank deficiency has already been detected.
   IERR = 21
   RETURN
END IF
IERR = 0
RETURN
END SUBROUTINE MAKEFIRSTSIMP

SUBROUTINE MAKESIMPLEX()
! Given a Delaunay facet F and a specified side of its containing hyperplane,
! complete the simplex by adding a point from PTS on the opposite `side' of
! F. Assume SIMP(1:D) contains the vertex indices of F (corresponding to
! data points P_1, P_2, ..., P_D in PTS) and SIMP(D+1) contains a point on
! the current side (to flip away from), and assume that the matrix A(1:D-1,:)
! and vector B(1:D-1) are filled appropriately (similarly as in
! MAKEFIRSTSIMP()). Then for any P* (not in the hyperplane containing F)
! in PTS, let CENTER denote the circumcenter of the simplex with vertices
! P_1, P_2, ..., P_D, P*. Then
!
! X = CENTER - P_1
!
! is given by the solution to the nonsingular linear system
!
! AX = B where
!
! A = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1, P* - P_1 ]^T and
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
! from PTS) of the facet F and SIMP(D+1) must contain the index of a point
! on the current side of F (to flip away from).  Upon output, SIMP(:) will
! contain the vertex indices of a Delaunay simplex on the opposite side of
! PLANE.  Also, the matrix A and vector B will be updated accordingly. If
! SIMP(D+1)=0, then there were no points in PTS on the appropriate side of
! F, meaning that F is a facet of the convex hull.

! Calculate the hyperplane boundary of the halfspace and store in PLANE.
CALL MAKEPLANE()
IF(IERR .NE. 0) RETURN
! Calculate the side of the viable point.
SIDE1 = -1.0_R8 * DDOT(D,PLANE(1:D),1,PTS(:,SIMP(D+1)),1) + PLANE(D+1)
! Normalize the magnitude of SIDE1.
SIDE1 = SIGN(1.0_R8,SIDE1)
! Initialize the center, radius, and simplex.
SIMP(D+1) = 0
CENTER = 0_R8
MINRAD = HUGE(MINRAD)
! Loop through all the points in PTS.
DO I = 1, N
   ! Check whether PTS(:,I) is on the viable side of halfspace.
   SIDE2 = DDOT(D,PLANE(1:D),1,PTS(:,I),1) - PLANE(D+1)
   IF (SIDE1 * SIDE2 < EPS .OR. ANY(SIMP(:) .EQ. I)) CYCLE
   ! Check for evidence of degeneracy.
   IF (ABS(DNRM2(D, PTS(:,I) - CENTER(:), 1) - MINRAD) < EPS) THEN
      IERR = 21
      RETURN
   END IF
   ! Check whether PTS(:,I) is inside the current ball.
   IF (DNRM2(D, PTS(:,I) - CENTER(:), 1) > MINRAD) CYCLE
   ! Add the point with index I to the linear system.
   A(D,:) = PTS(:,I) - PTS(:,SIMP(1))
   B(D) = DDOT(D, A(D,:), 1, A(D,:), 1) / 2.0_R8
   LQ = A
   CENTER = B
   ! Solve the linear system to get the center.
   CALL DGESV(D, 1, LQ, D, PIV, CENTER, D, IERR)
   IF (IERR < 0) THEN
      IERR = 90
      RETURN
   ELSE IF (IERR > 0) THEN
      IERR = 21
      RETURN
   END IF
   ! Update the radius, center, and simplex.
   CURRRAD = DNRM2(D, CENTER, 1)
   MINRAD = CURRRAD
   CENTER = CENTER + PTS(:,SIMP(1))
   SIMP(D+1) = I
END DO
! Reset the error flag.
IERR = 0
! Check for the extrapolation condition.
IF(SIMP(D+1) .EQ. 0) RETURN
! If there was no extrapolation, add the new point to the linear system.
A(D,:) = PTS(:,SIMP(D+1)) - PTS(:,SIMP(1))
B(D) = DDOT(D, A(D,:), 1, A(D,:), 1) / 2.0_R8
RETURN
END SUBROUTINE MAKESIMPLEX

SUBROUTINE MAKEPLANE()
! Construct a plane c^T x = \alpha containing the first D vertices indexed
! in SIMP(:). The plane is determined by its normal vector c and \alpha.
! Let (P_1, P_2, ..., P_D) be the vertices indexed in SIMP(1:D). A normal
! vector is any nonzero vector in ker(A), where the matrix
! 
! A = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1 ]^T.
! 
! Since rank A = D-1, dim ker(A) = 1, and ker(A) can be found from a QR
! factorization of A:  AP = QR, where P permutes the columns of A. Solving
! AP X = QR X = 0 with X_D =  1 gives a normal vector PX.
! 
! Upon output, PLANE(1:D) contains the normal vector c and PLANE(D+1)
! contains \alpha defining the plane.

! CHECK THAT D-1 > 0.
IF (D > 1) THEN
   ! Compute the QR factorization.
   PIV=0
   LQ = A
   CALL DGEQP3(D-1,D,LQ,D,PIV,PLANE,WORK,LWORK,IERR)
   IF(IERR < 0) THEN
      IERR = 90
      RETURN
   ELSE IF (IERR > 0) THEN
      IERR = 21
      RETURN
   END IF
   ! Perform a triangular solve, fixing the final value to get
   ! a nontrivial solution.
   WORK(1:D-1) = LQ(1:D-1,D)
   CALL DTRSM('L','U','N','N',D-1,1,-1.0_R8,LQ,D,WORK,D)
   WORK(D) = 1.0_R8
   ! Undo the pivots.
   DO I = 1,D
      PLANE(PIV(I)) = WORK(I)
   END DO
   ! Normalize the orthogonal vector.
   PLANE(1:D) = PLANE(1:D) / DNRM2(D,PLANE(1:D),1)
   ! Calculate the intercept.
   PLANE(D+1) = DDOT(D,PLANE(1:D),1,PTS(:,SIMP(1)),1)
ELSE
   ! Compute the hyper plane when D=1.
   PLANE(1) = 1.0_R8
   PLANE(2) = PTS(1,SIMP(1))
END IF
RETURN
END SUBROUTINE MAKEPLANE

SUBROUTINE SelectSortSimp()
! Sort simplex indices, but keep the top item at the bottom.
! O(D^2) complexity.
DO I=1, D+1
   J = MINLOC(SIMP(I:D+1),DIM=1) + I - 1
   ITMP = SIMP(J)
   SIMP(J) = SIMP(I)
   SIMP(I) = ITMP
END DO
RETURN
END SUBROUTINE SELECTSORTSIMP

SUBROUTINE BUBBLESORTSIMP()
! Sort simplex indices given that only the last index is out of place and
! the first index is fixed. O(D) complexity.
DO I=D+1, 1, -1
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
