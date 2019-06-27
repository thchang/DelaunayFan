MODULE AFL
! This module implements the active face list as described in:
!
! T. H. Chang, L. T. Watson, T. C.H. Lux, S. Raghvendra, B. Li, L. Xu,
! A. R. Butt, K. W. Cameron, and Y. Hong. Computing the umbrella
! neighbourhood of a vertex in the Delaunay triangulation and a single
! Voronoi cell in arbitrary dimension. In Proceedings of IEEE Southeast
! Conference 2018 (SEC 2018). IEEE, St. Petersburg, FL, 2018.
!
! It is specifically designed for usage in the VTDELAUNAYFAN module.
!
! Primary author: Tyler H. Chang
! Last Update : June, 2019
USE VECTOR_CLASS, ONLY : INTVECTOR
IMPLICIT NONE

PRIVATE
PUBLIC :: FACELIST, INTVECTOR

! The face list class extends the INTVECTOR class for dynamic memory.
TYPE, EXTENDS(INTVECTOR) :: FACELIST
   INTEGER, PRIVATE :: IND ! store the index of the central vertex
CONTAINS
   PROCEDURE, NOPASS :: CONSTRUCTOR => NEWFACELIST
   PROCEDURE, PASS :: CONTAIN => FACELISTCONTAINS
   PROCEDURE, PASS :: PUSHSIMP => PUSHSIMPLEX
   PROCEDURE, PASS :: TOP => FACELISTTOP
   PROCEDURE, PRIVATE, PASS :: PF => PUSHFACE
END TYPE FACELIST

! Define an interface for the constructor.
INTERFACE FACELIST
   MODULE PROCEDURE NEWFACELIST
END INTERFACE FACELIST

CONTAINS

FUNCTION NEWFACELIST(DIM, IND, ISTAT) RESULT(FL)
! This is the constructor for FACELIST class.
! 
! On input, DIM specifies the dimension of the FACELIST and IND specifies
! the index of the central point. For a D-dimensional space, DIM should
! be D+1 because an additional entry is required to store the side of
! each facet.
!
! The optional argument ISTAT returns an error flag.
!
! The result is a new FACELIST object.

   ! Parameters.
   INTEGER, INTENT(IN) :: DIM, IND
   INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
   TYPE(FACELIST) :: FL
   ! Set metadata.
   FL%N = 0
   FL%MAXLEN = 0
   FL%IND = IND
   ! Initialize ISTAT, if present.
   IF(PRESENT(ISTAT)) ISTAT = 0
   ! Check for input errors and set the dimension of the FACELIST.
   IF (DIM < 1) THEN
      IF(PRESENT(ISTAT)) ISTAT = -1
      RETURN
   END IF
   FL%D = DIM
   RETURN
END FUNCTION NEWFACELIST

SUBROUTINE FACELISTTOP(FL, SIMP, ISTAT)
! Get the top item in FL without deleting it.

   ! Parameters.
   CLASS(FACELIST), INTENT(IN) :: FL
   INTEGER, INTENT(OUT) :: SIMP(:)
   INTEGER, INTENT(OUT) :: ISTAT

   ! Check for illegal size.
   IF(SIZE(SIMP) .NE. FL%D) THEN
      ISTAT = -1
      RETURN
   END IF   
   ! Otherwise, return the data.
   IF(FL%N .EQ. 0) THEN
      SIMP(:) = 0
   ELSE
      SIMP(:) = FL%DAT(:,FL%N)
   END IF
   RETURN
END SUBROUTINE FACELISTTOP

SUBROUTINE PUSHSIMPLEX(FL, SIMP, ISTAT)
! This subroutine `pushes' all the facets of a simplex to FL.
!
! On input, FL is the face list and SIMP is the simplex to push.
!
! On output, FL has been updated following the collision rule in
! Chang et. al., and ISTAT contains an error flag (successful push = 0).

   ! Parameters.
   CLASS(FACELIST), INTENT(INOUT) :: FL
   INTEGER, INTENT(INOUT) :: SIMP(FL%D)
   INTEGER, INTENT(OUT) :: ISTAT
   ! Local variables.
   INTEGER :: I, TMP

   ! Loop over all entries.
   DO I = FL%D, 2, -1
      ! Move the unused vertex to index D.
      TMP = SIMP(FL%D)
      SIMP(FL%D) = SIMP(I)
      SIMP(I) = TMP
      ! Push the facet defined by the first D vertices to FL.
      CALL FL%PF(SIMP, ISTAT)
      IF (ISTAT .NE. 0) RETURN
   END DO
   RETURN
END SUBROUTINE PUSHSIMPLEX

SUBROUTINE PUSHFACE(FL, FACE, ISTAT)
! This subroutine implements the `push' operation.
!
! On input, FL is the face list and FACE is the face to push.
!
! On output, FL has been updated following the collision rule in
! Chang et. al., and ISTAT contains an error flag (successful push = 0).

   ! Parameters.
   CLASS(FACELIST), INTENT(INOUT) :: FL
   INTEGER, INTENT(IN) :: FACE(:)
   INTEGER, INTENT(OUT) :: ISTAT
   ! Local variables.
   INTEGER :: I

   ! Check if FACE(:) is an appropriate dimension.
   IF (SIZE(FACE,1) .NE. FL%D) THEN; ISTAT=-1; RETURN; END IF
   ! Check if FACE(:) contains IND.
   IF (.NOT. ANY(FACE(1:FL%D) .EQ. FL%IND)) RETURN
   ! Look for collisions.
   I = FL%CONTAIN(FACE, ISTAT)
   IF (ISTAT .NE. 0) RETURN
   ! If a collision is detected, delete the entry.
   IF (I > 0) THEN
      CALL FL%DELETE(I, ISTAT)
      IF (ISTAT .NE. 0) RETURN
   ! If there was no collision, push the new entry.
   ELSE
      CALL FL%PUSH(FACE, ISTAT)
      IF (ISTAT .NE. 0) RETURN
   END IF
   RETURN 
END SUBROUTINE PUSHFACE

FUNCTION FACELISTCONTAINS(FL, FACET, ISTAT) RESULT(IND)
! The FACELISTCONTAINS subroutine searches FL for FACET.
!
! Required arguments:
!
! On input, FL is an initialized FACELIST object with leading dimension D,
! and FACET(1:D-1) is an integer array containing the facet to search for.
!
! If FACET is found, FACELISTCONTAINS returns the index. Otherwise, 0 is
! returned.
!
! ISTAT returns an error code associated with execution (0 for success).

! Parameter list.
CLASS(FACELIST), INTENT(IN) :: FL
INTEGER, INTENT(IN) :: FACET(:)
INTEGER, INTENT(OUT) :: ISTAT
! Return value.
INTEGER :: IND
! Local variables.
INTEGER :: I

! Check for an illegal size.
IF(FL%D-1 .NE. SIZE(FACET)) THEN
   ISTAT = -1
   RETURN
END IF
! Search an unsorted list.
IND = 0
DO I=1, FL%N
   IF(ALL(FL%DAT(1:D-1,I) .EQ. FACET(:))) THEN
      IND = I
      EXIT
   END IF
END DO
RETURN
END FUNCTION FACELISTCONTAINS

END MODULE AFL
