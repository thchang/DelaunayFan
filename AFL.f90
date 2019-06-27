MODULE AFL
! This module implements the active facet list as described in:
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
! Last Update : 2018
USE VECTOR
IMPLICIT NONE

! The face list will be stored in dynamic memory.
TYPE FACELIST
   INTEGER :: IND, D
   TYPE(INTVECTOR) :: TABLE
END TYPE FACELIST

CONTAINS

FUNCTION NEWFACELIST(IND, D, IERR) RESULT(FL)
   INTEGER :: IND, D, IERR
   TYPE(FACELIST) :: FL
   CALL NEWINTVECTOR(FL%TABLE, IERR, DIM=D+1)
   IF (IERR .NE. 0) RETURN
   FL%D = D
   FL%IND = IND
   RETURN
END FUNCTION

SUBROUTINE PUSHFACE(FL, FACE, IERR)
! This subroutine implements the `push' operation.
!
! On input, FL is the face list and FACE is the face to push.
!
! On output, FL has been updated following the collision rule in
! Chang et. al., and IERR contains an error flag (successful push = 0).

   TYPE(FACELIST) :: FL
   INTEGER :: FACE(FL%D+1), IERR, I

   ! Check if the face contains the vertex.
   IF (.NOT. ANY(FACE(1:FL%D) .EQ. FL%IND)) RETURN
   ! Look for collisions.
   I = INTVECTORCONTAINS(FL%TABLE, FACE, IERR)
   IF (IERR .NE. 0) RETURN
   ! If a collision is detected, delete the entry.
   IF (I > 0) THEN
      CALL INTVECTORDEL(FL%TABLE, I, IERR)
      IF (IERR .NE. 0) RETURN
   ! If there was no collision, push the new entry.
   ELSE
      CALL INTVECTORPUSH(FL%TABLE, FACE, IERR)
      IF (IERR .NE. 0) RETURN
   END IF
   
   RETURN 
END SUBROUTINE PUSHFACE

SUBROUTINE PUSHALL(FL, SIMP, IERR)
! This subroutine `pushes' all the faces in a simplex.
!
! On input, FL is the face list and SIMP is the simplex to push.
!
! On output, FL has been updated following the collision rule in
! Chang et. al., and IERR contains an error flag (successful push = 0).

   TYPE(FACELIST) :: FL
   INTEGER :: SIMP(FL%D+1), IERR, I, TMP

   DO I = FL%D+1, 2, -1
      ! Move the unused vertex to index D+1.
      TMP = SIMP(FL%D+1)
      SIMP(FL%D+1) = SIMP(I)
      SIMP(I) = TMP
      ! Push the facet defined by the first D vertices to FL.
      CALL PUSHFACE(FL, SIMP, IERR)
      IF (IERR .NE. 0) RETURN
   END DO
END SUBROUTINE PUSHALL

SUBROUTINE POPFACE(FL, FACE, IERR)
! This subroutine `pops' the top facet in FL off and returns it in FACE.
!
! On input, FL is the face list.
!
! On output, FL has had its top item removed and it is contained in FACE,
! and IERR contains an error flag (successful push = 0).

   TYPE(FACELIST) :: FL
   INTEGER :: FACE(FL%D+1), IERR
   IERR = 0
   IF(FL%TABLE%LENGTH .EQ. 0) THEN
      FACE = 0
   ELSE
      FACE(:) = FL%TABLE%DAT(:,FL%TABLE%LENGTH)
   END IF
   RETURN
END SUBROUTINE POPFACE

END MODULE AFL
