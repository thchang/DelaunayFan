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

PRIVATE
PUBLIC :: FACELIST, INTVECTOR

! The INTVECTOR type defines a dynamic 2D INTEGER array.
TYPE INTVECTOR
   ! Private meta data.
   INTEGER, PRIVATE :: D, N, MAXLEN
   ! Private data
   INTEGER, ALLOCATABLE, PRIVATE :: DAT(:,:)
CONTAINS
   PROCEDURE, PUBLIC, PASS :: LENGTH => VECTORSIZE
   PROCEDURE, PUBLIC, PASS :: DIM => VECTORDIM
   PROCEDURE, PUBLIC, PASS :: PUSH => INTVECTORPUSH
   PROCEDURE, PUBLIC, PASS :: POP => INTVECTORPOP
   PROCEDURE, PUBLIC, PASS :: INSERT => INTVECTORINS
   PROCEDURE, PUBLIC, PASS :: DELETE => INTVECTORDEL
   PROCEDURE, PUBLIC, PASS :: SET => INTVECTORSET
   PROCEDURE, PUBLIC, PASS :: ITEM => INTVECTORITEM
   PROCEDURE, PUBLIC, PASS :: DATA => INTVECTORDATA
   PROCEDURE, PUBLIC, PASS :: FREE => INTVECTORFREE
END TYPE INTVECTOR

! The FACELIST class extends the INTVECTOR class for dynamic memory.
TYPE, EXTENDS(INTVECTOR) :: FACELIST
   INTEGER, PRIVATE :: IND ! store the index of the central vertex
CONTAINS
   PROCEDURE, PUBLIC, NOPASS :: CONSTRUCTOR => NEWFACELIST
   PROCEDURE, PUBLIC, PASS :: CONTAIN => FACELISTCONTAINS
   PROCEDURE, PUBLIC, PASS :: PUSHSIMP => PUSHSIMPLEX
   PROCEDURE, PUBLIC, PASS :: TOP => FACELISTTOP
   PROCEDURE, PRIVATE, PASS :: PF => PUSHFACE
END TYPE FACELIST

! Define an interface for the INTVECTOR constructor.
INTERFACE INTVECTOR
   MODULE PROCEDURE NEWINTVECTOR
END INTERFACE INTVECTOR

! Define an interface for the FACELIST constructor.
INTERFACE FACELIST
   MODULE PROCEDURE NEWFACELIST
END INTERFACE FACELIST

CONTAINS

! These are the procedures associated exclusively with the FACELIST class,
! which extends the INTVECTOR class. The INTVECTOR class procedures are
! further below.

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
IMPLICIT NONE
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
IMPLICIT NONE
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
! On input, FL is the face list and SIMP is the simplex to push.
! On output, FL has been updated following the collision rule in
! Chang et. al., and ISTAT contains an error flag (successful push = 0).
IMPLICIT NONE
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

SUBROUTINE PUSHFACE(FL, FACET, ISTAT)
! This subroutine implements the `push' operation.
! On input, FL is the face list and FACET is the facet to push.
! On output, FL has been updated following the collision rule in
! Chang et. al., and ISTAT contains an error flag (successful push = 0).
IMPLICIT NONE
! Parameters.
CLASS(FACELIST), INTENT(INOUT) :: FL
INTEGER, INTENT(IN) :: FACET(:)
INTEGER, INTENT(OUT) :: ISTAT
! Local variables.
INTEGER :: I
! Check if FACET(:) is an appropriate dimension.
IF (SIZE(FACET) .NE. FL%D) THEN; ISTAT=-1; RETURN; END IF
! Check if FACET(:) contains IND.
IF (.NOT. ANY(FACET(1:FL%D) .EQ. FL%IND)) RETURN
! Look for collisions.
I = FL%CONTAIN(FACET(1:FL%D-1), ISTAT)
IF (ISTAT .NE. 0) RETURN
! If a collision is detected, delete the entry.
IF (I > 0) THEN
   CALL FL%DELETE(I, ISTAT)
   IF (ISTAT .NE. 0) RETURN
! If there was no collision, push the new entry.
ELSE
   CALL FL%PUSH(FACET, ISTAT)
   IF (ISTAT .NE. 0) RETURN
END IF
RETURN 
END SUBROUTINE PUSHFACE

FUNCTION FACELISTCONTAINS(FL, FACET, ISTAT) RESULT(IND)
! The FACELISTCONTAINS subroutine searches FL for FACET.
! On input, FL is an initialized FACELIST object with leading dimension D,
! and FACET(1:D-1) is an integer array containing the facet to search for.
! If FACET is found, FACELISTCONTAINS returns the index. Otherwise, 0 is
! returned.
! ISTAT returns an error code associated with execution (0 for success).
IMPLICIT NONE
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
   IF(ALL(FL%DAT(1:FL%D-1,I) .EQ. FACET(:))) THEN
      IND = I
      EXIT
   END IF
END DO
RETURN
END FUNCTION FACELISTCONTAINS

! Define the procedures for the INTVECTOR class below.
! All procedures were lifted from the VECTOR_CLASS module on 6/28/2019.
! See https://github.com/thchang/VectorClass.

FUNCTION NEWINTVECTOR(DIM, ISTAT) RESULT(VEC)
! This is the INTVECTOR constructor.
! DIM (optional) contains the INTEGER value of the first (static)
! dimension of the INTVECTOR object. If not present, DIM is assumed
! to be 1.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
INTEGER, OPTIONAL, INTENT(IN) :: DIM
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
TYPE(INTVECTOR) :: VEC
! Set metadata.
VEC%N = 0
VEC%MAXLEN = 0
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check for optional input.
IF(PRESENT(DIM)) THEN
   IF (DIM < 1) THEN
      IF(PRESENT(ISTAT)) ISTAT = -1
      RETURN
   END IF
   VEC%D = DIM
ELSE
   VEC%D = 1
END IF
RETURN
END FUNCTION NEWINTVECTOR

FUNCTION VECTORSIZE(VEC) RESULT(SIZE)
! This function returns the size of the second (dynamic) dimension
! of an INTVECTOR object. I.e., it's length.
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(IN) :: VEC
INTEGER :: SIZE
! Get the size of the INTVECTOR object.
SIZE = VEC%N
RETURN
END FUNCTION VECTORSIZE

FUNCTION VECTORDIM(VEC) RESULT(DIM)
! This function returns the size of the first (static) dimension
! of an INTVECTOR object.
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(IN) :: VEC
INTEGER :: DIM
! Get the dimension of the INVECTOR object.
DIM = VEC%D
RETURN
END FUNCTION VECTORDIM

SUBROUTINE INTVECTORPUSH(VEC, ITEM, ISTAT)
! This subroutine pushes the contents of ITEM to the next open
! entry in an INTVECTOR object, resizing if necessarry.
! ITEM is an INTEGER array whose dimension matches VEC%DIM(). On
! input, ITEM(:) contains the value to append to the INTVECTOR.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! Local data.
INTEGER :: TMP(VEC%D,VEC%N), LSTAT
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check that ITEM is a legal size.
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! No resizing of the INTVECTOR is required.
IF(VEC%N < VEC%MAXLEN) THEN
   ! Increment the counter and append ITEM(:) to the INTVECTOR.
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
! The INTVECTOR needs to be initialized.
ELSE IF(VEC%MAXLEN .EQ. 0) THEN
   ! Allocate the INTVECTOR.
   ALLOCATE(VEC%DAT(VEC%D,8), STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 8
   ! Increment the counter and append ITEM(:) to the INTVECTOR.
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
! The INTVECTOR needs to be resized.
ELSE
   ! Make a copy of the current data, and reallocate the INTVECTOR.
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ALLOCATE(VEC%DAT(VEC%D,VEC%MAXLEN*2), &
      STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN*2
   ! Restore the copied data.
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
   ! Increment the counter and append ITEM(:) to the INTVECTOR.
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
END IF
RETURN
END SUBROUTINE INTVECTORPUSH

SUBROUTINE INTVECTORPOP(VEC, ITEM, ISTAT)
! This subroutine pops the top entry off of the INTVECTOR object,
! resizing if necessarry.
! ITEM(:) is an INTEGER array whose dimension matches VEC%DIM().
! On output, ITEM(:) contains the data that was popped from the
! INTVECTOR object. If the INTVECTOR was empty, then ITEM(:) = 0.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(OUT) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! Local variables.
INTEGER :: TMP(VEC%D, VEC%N-1), LSTAT
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check that ITEM is a legal size.
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! If the INTVECTOR is empty, then return 0.
IF(VEC%N .EQ. 0) THEN
   ITEM = 0
   RETURN
END IF
! Pop ITEM(:) off of the top of the INTVECTOR.
ITEM(:) = VEC%DAT(:,VEC%N)
VEC%N = VEC%N - 1
! Reclaim memory if the INTVECTOR%SIZE() underflows the MAXLEN.
IF (VEC%N < VEC%MAXLEN / 2 .AND. VEC%N > 6) THEN
   ! Make a copy of the current data, and reallocate the INTVECTOR.
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT=LSTAT
      RETURN
   END IF
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT=LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! Restore the copied data.
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
! If the INTVECTOR is now empty, then reclaim all the memory.
ELSE IF (VEC%N .EQ. 0) THEN
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT=LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE INTVECTORPOP

SUBROUTINE INTVECTORINS(VEC, ITEM, IND, ISTAT)
! This subroutine inserts the contents of ITEM to the next open
! entry in an INTVECTOR object, resizing if necessarry.
! ITEM(:) is an INTEGER array whose dimension matches VEC%DIM(). On
! input, ITEM(:) contains the value to insert into the INTVECTOR.
! IND is the INTEGER index at which to insert ITEM(:).
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: ITEM(:), IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! Local data.
INTEGER :: TMP(VEC%D,VEC%N), LSTAT
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check that ITEM is a legal size.
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! Check that IND is in an appropriate range.
IF(IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -2
   RETURN
END IF
! No resizing of the INTVECTOR is required.
IF(VEC%N < VEC%MAXLEN) THEN
   ! Increment the counter and insert ITEM(:) into the INTVECTOR.
   VEC%DAT(:,IND+1:VEC%N+1) = VEC%DAT(:,IND:VEC%N)
   VEC%DAT(:,IND) = ITEM(:)
   VEC%N = VEC%N + 1
! The INTVECTOR needs to be resized.
ELSE
   ! Make a copy of the current data, and reallocate the INTVECTOR.
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ALLOCATE(VEC%DAT(VEC%D,VEC%MAXLEN*2), STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN=VEC%MAXLEN*2
   ! Restore the copied data.
   VEC%DAT(:,1:IND-1) = TMP(:,1:IND-1)
   VEC%DAT(:,IND+1:VEC%N+1) = TMP(:,IND:VEC%N)
   ! Increment the counter and insert ITEM(:) into the INTVECTOR.
   VEC%N = VEC%N + 1
   VEC%DAT(:,IND) = ITEM(:)
END IF
RETURN
END SUBROUTINE INTVECTORINS

SUBROUTINE INTVECTORDEL(VEC, IND, ISTAT)
! This subroutine removes the element of the INTVECTOR object located
! at the index IND, resizing if necessarry.
! IND is the INTEGER index of the element to delete from the INTVECTOR.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! Local data.
INTEGER :: TMP(VEC%D,VEC%N-1), LSTAT
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check that IND is in an appropriate range.
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! Shift all other elements down to overwrite the data at IND.
VEC%DAT(:,IND:VEC%N-1) = VEC%DAT(:,IND+1:VEC%N)
VEC%N = VEC%N - 1
! Reclaim memory if the INTVECTOR%SIZE() underflows the MAXLEN.
IF (VEC%N < VEC%MAXLEN / 2 .AND. &
  VEC%N > 6) THEN
   ! Make a copy of the current data, and reallocate the INTVECTOR.
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! Restore the copied data.
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
! If the INTVECTOR is now empty, then reclaim all the memory.
ELSE IF (VEC%N .EQ. 0) THEN
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE INTVECTORDEL

SUBROUTINE INTVECTORSET(VEC, IND, ITEM, ISTAT)
! This subroutine sets the value of an INTVECTOR at some index.
! IND is the INTEGER index at which to set the value of INTVECTOR to
! that contained in ITEM(:).
! ITEM(:) is an INTEGER array whose dimension matches VEC%DIM(). On
! input, ITEM(:) contains the value to set in the INTVECTOR.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, INTENT(IN) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check that IND is in an appropriate range.
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! Check that ITEM is a legal size.
IF (VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -2
   RETURN
END IF
! Set the value of VEC%DAT(:,IND) to that of ITEM(:).
VEC%DAT(:,IND) = ITEM(:)
RETURN
END SUBROUTINE INTVECTORSET

FUNCTION INTVECTORITEM(VEC, IND, ISTAT) RESULT(ITEM)
! Function for accessing a single item in the INTVECTOR object.
! The function returns a 1D INTEGER array containing the value(s) of
! the item in the INTVECTOR that was located at the specified index.
! IND is the INTEGER index of the item to be accessed.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(IN) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
INTEGER :: ITEM(VEC%D)
! Initialize ISTAT.
IF(PRESENT(ISTAT)) ISTAT = 0
! Check that IND is in an appropriate range.
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! Get the value of VEC%DAT(:,IND), and return in ITEM(:).
ITEM = VEC%DAT(:,IND)
RETURN
END FUNCTION INTVECTORITEM

FUNCTION INTVECTORDATA(VEC) RESULT(DAT)
! Function for accessing all the data in the INTVECTOR object.
! The function returns a 2D INTEGER array (matrix) containing all
! values stored in the INTVECTOR object.
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(IN) :: VEC
INTEGER :: DAT(VEC%D, VEC%N)
! Get data all data in the INTVECTOR and return in DAT(:,:).
DAT(:,:) = VEC%DAT(:,1:VEC%N)
RETURN
END FUNCTION INTVECTORDATA

SUBROUTINE INTVECTORFREE(VEC, ISTAT)
! Free all memory associated with the INTVECTOR object.
! ISTAT (optional) is an INTEGER error flag (Success = 0).
IMPLICIT NONE
! Parameters.
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! Local variables.
INTEGER :: LSTAT
! Free the heap memory, record an errors in LSTAT and pass to ISTAT
! (if present).
IF (ALLOCATED(VEC%DAT)) THEN
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF(PRESENT(ISTAT)) ISTAT=LSTAT
END IF
RETURN
END SUBROUTINE INTVECTORFREE

END MODULE AFL
