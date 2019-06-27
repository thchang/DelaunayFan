! FILE : Vector.f90
! AUTHOR: Tyler H. Chang
! LAST UPDATE: June, 2019
! ISO FORTRAN 2003 STANDARD
! 
! ABSTRACT:
!
! Module for FORTRAN vector (dynamic array) operations.
! A Vector is a 2-Dimensional Array where the 1st
! dimension is static and the 2nd dimension is dynamic.
! Options for both 32-bit INTEGER (IntVector) and 64-bit
! REAL (RealVector) types.
! Also defines KIND=R8 for ~64-bit arithmetic on most
! standard machines.
!
! The dynamic memory allocation scheme is as follows.
! If VECTOR%LENGTH == 0, then no memory is allocated.
! If VECTOR%LENGTH overflows VECTOR%MAXLEN, then double
! the size of VECTOR%MAXLEN (unless VECTOR%MAXLEN == 0,
! then initialize to VECTOR%MAXLEN = 8).
! If VECTOR%LENGTH underflows VECTOR%MAXLEN, occupying
! less than half the total storage, reclaim half the
! size allocated to VECTOR.
! 
! USAGE:
!
! The INTVECTOR and REALVECTOR types serve as
! 2-dimensional arrays for INTEGER and REAL(KIND=R8)
! types respectively. 
! The R8 kind is defined herein and gives ~64-bit 
! arithmetic on all known machines.
!
! To initialize a INTVECTOR type, use:
!    TYPE(INTVECTOR) VECTOR
!    VECTOR = INTVECTOR([DIM [, ISTAT]])
! where DIM specifies the size of the 1st (fixed)
! dimension. Similarly, to initilaize an REALVECTOR,
! use:
!    TYPE(REALVECTOR) VECTOR
!    VECTOR = REALVECTOR([DIM [, ISTAT]])
! From hereon, WLOG, both types are referenced as
! 'VECTOR'.
!
! All subroutines related to the Vector class use
! an optional error flag : ISTAT to report successful
! operations.
! The scheme is as follows:
! ISTAT = 0 : SUCCESS
! ISTAT < 0 : ILLEGAL INPUT AT ARG : (-ISTAT)
! ISTAT > 0 : MEMORY ALLOCATION/DEALLOCATION ERROR
!
! The following subroutines are used to manage and
! access the data in a VECTOR object. In all of the
! following, ITEM is either a REAL(KIND=R8) or INTEGER
! valued array, whose length matches the first (static)
! dimension of VECTOR, and IND is an INTEGER index whose
! value is between 1 and VECTOR%LENGTH.
!
! VECTOR%PUSH(ITEM [, ISTAT]) pushes ITEM to the end
! of VECTOR, resizing if necessarry.
!
! VECTOR%POP(ITEM [, ISTAT]) pops the top item off the
! back of VECTOR, resizing if necessarry, and returns
! the value in ITEM.
!
! VECTOR%INSERT(ITEM, IND [, ISTAT]) inserts ITEM at
! location IND.
!
! VECTOR%DELETE(IND [, ISTAT]) deletes the item at
! location IND.
!
! VECTOR%SET(IND, ITEM [, ISTAT]) sets the value of
! the data at location IND to ITEM.
!
! VECTOR%ITEM(IND [, ISTAT]) returns the data at
! location IND as an array.
!
! VECTOR%DATA() returns the entire data container at
! for VECTOR (a 2-dimensional array).
!
! VECTOR%LENGTH() returns the length of the second
! (dynamic) dimension.
!
! VECTOR%DIM() returns the length of the first (static)
! dimension.
!
! VECTOR%FREE() frees all allocated memory associated
! with VECTOR, resetting it to a new object.

MODULE VECTOR_CLASS

PRIVATE
PUBLIC :: INTVECTOR, REALVECTOR, R8

! Hompack type for 64-bit real numbers
INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(13)

! Abstract Vector type
TYPE, ABSTRACT :: VECTOR
   ! private meta data
   INTEGER, PRIVATE :: D, N, MAXLEN
CONTAINS
   PROCEDURE, PASS :: LENGTH => VECTORSIZE
   PROCEDURE, PASS :: DIM => VECTORDIM
END TYPE VECTOR

! Integer Vector type (dynamic array)
TYPE, EXTENDS(VECTOR) :: INTVECTOR
   ! private data
   INTEGER, ALLOCATABLE, PRIVATE :: DAT(:,:)
CONTAINS
   PROCEDURE, PASS :: PUSH => INTVECTORPUSH
   PROCEDURE, PASS :: POP => INTVECTORPOP
   PROCEDURE, PASS :: INSERT => INTVECTORINS
   PROCEDURE, PASS :: DELETE => INTVECTORDEL
   PROCEDURE, PASS :: SET => INTVECTORSET
   PROCEDURE, PASS :: ITEM => INTVECTORITEM
   PROCEDURE, PASS :: DATA => INTVECTORDATA
   PROCEDURE, PASS :: FREE => INTVECTORFREE
END TYPE INTVECTOR

! Real 64-bit Vector type (dynamic array)
TYPE, EXTENDS(VECTOR) :: REALVECTOR
   ! data
   REAL(KIND=R8), ALLOCATABLE, PRIVATE :: DAT(:,:)
CONTAINS
   PROCEDURE, PASS :: PUSH => REALVECTORPUSH
   PROCEDURE, PASS :: POP => REALVECTORPOP
   PROCEDURE, PASS :: INSERT => REALVECTORINS
   PROCEDURE, PASS :: DELETE => REALVECTORDEL
   PROCEDURE, PASS :: SET => REALVECTORSET
   PROCEDURE, PASS :: ITEM => REALVECTORITEM
   PROCEDURE, PASS :: DATA => REALVECTORDATA
   PROCEDURE, PASS :: FREE => REALVECTORFREE
END TYPE REALVECTOR

! Constructor interfaces.
INTERFACE INTVECTOR
   MODULE PROCEDURE NEWINTVECTOR
END INTERFACE INTVECTOR

INTERFACE REALVECTOR
   MODULE PROCEDURE NEWREALVECTOR
END INTERFACE REALVECTOR

CONTAINS

! --- Generic procedures for Vector base type --- !

! check vector length
FUNCTION VECTORSIZE(VEC) RESULT(SIZE)
IMPLICIT NONE
! parameters
CLASS(VECTOR), INTENT(IN) :: VEC
INTEGER :: SIZE
! get size
SIZE = VEC%N
RETURN
END FUNCTION VECTORSIZE

! check vector length
FUNCTION VECTORDIM(VEC) RESULT(DIM)
IMPLICIT NONE
! parameters
CLASS(VECTOR), INTENT(IN) :: VEC
INTEGER :: DIM
! get size
DIM = VEC%D
RETURN
END FUNCTION VECTORDIM

! --- Procedures for IntVector type --- !

! constructor
FUNCTION NEWINTVECTOR(DIM, ISTAT) RESULT(VEC)
IMPLICIT NONE
! parameters
INTEGER, OPTIONAL, INTENT(IN) :: DIM
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
TYPE(INTVECTOR) :: VEC
! set information
VEC%N = 0
VEC%MAXLEN = 0
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for optionals
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

! push function
SUBROUTINE INTVECTORPUSH(VEC, ITEM, ISTAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
INTEGER :: TMP(VEC%D,VEC%N), LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for illegal size
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! no resizing required
IF(VEC%N < VEC%MAXLEN) THEN
   ! increment counter and add
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
! needs to be initialized
ELSE IF(VEC%MAXLEN .EQ. 0) THEN
   ! allocate
   ALLOCATE(VEC%DAT(VEC%D,8), STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 8
   ! increment counter and add
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
! resize data
ELSE
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! reallocate
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
   ! restore current data
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
   ! add new data
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
END IF
RETURN
END SUBROUTINE INTVECTORPUSH

! pop function
SUBROUTINE INTVECTORPOP(VEC, ITEM, ISTAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(OUT) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local variables
INTEGER :: TMP(VEC%D, VEC%N-1), LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for illegal size
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! check for empty vector
IF(VEC%N .EQ. 0) THEN
   ITEM = 0
   RETURN
END IF
! get item
ITEM(:) = VEC%DAT(:,VEC%N)
VEC%N = VEC%N - 1
! reclaim partial memory
IF (VEC%N < VEC%MAXLEN / 2 .AND. &
  VEC%N > 6) THEN
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT=LSTAT
      RETURN
   END IF
   ! reallocate
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT=LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! copy back
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
! reclaim all memory
ELSE IF (VEC%N .EQ. 0) THEN
   ! FREE
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT=LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE INTVECTORPOP

! insert function
SUBROUTINE INTVECTORINS(VEC, ITEM, IND, ISTAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: ITEM(:), IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
INTEGER :: TMP(VEC%D,VEC%N), LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for illegal size
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! check for bad index
IF(IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -2
   RETURN
END IF
! no resizing required
IF(VEC%N < VEC%MAXLEN) THEN
   ! increment counter and insert
   VEC%DAT(:,IND+1:VEC%N+1) = &
      VEC%DAT(:,IND:VEC%N)
   VEC%DAT(:,IND) = ITEM(:)
   VEC%N = VEC%N + 1
! resize data
ELSE
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! reallocate
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
   ! restore current data
   VEC%DAT(:,1:IND-1) = TMP(:,1:IND-1)
   VEC%DAT(:,IND+1:VEC%N+1) = TMP(:,IND:VEC%N)
   ! insert new data
   VEC%N = VEC%N + 1
   VEC%DAT(:,IND) = ITEM(:)
END IF
RETURN
END SUBROUTINE INTVECTORINS

! delete function
SUBROUTINE INTVECTORDEL(VEC, IND, ISTAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
INTEGER :: TMP(VEC%D,VEC%N-1), LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for bad input
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! shift other elements down
VEC%DAT(:,IND:VEC%N-1) = VEC%DAT(:,IND+1:VEC%N)
VEC%N = VEC%N - 1
! reclaim partial memory
IF (VEC%N < VEC%MAXLEN / 2 .AND. &
  VEC%N > 6) THEN
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ! reallocate
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! copy back
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
! reclaim all memory
ELSE IF (VEC%N .EQ. 0) THEN
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE INTVECTORDEL

! set the value of VEC at IND to ITEM
SUBROUTINE INTVECTORSET(VEC, IND, ITEM, ISTAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, INTENT(IN) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for bad input
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
IF (VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -2
   RETURN
END IF
! set VEC%DAT(:,IND) to ITEM
VEC%DAT(:,IND) = ITEM(:)
RETURN
END SUBROUTINE INTVECTORSET

! memory access function (single item)
FUNCTION INTVECTORITEM(VEC, IND, ISTAT) RESULT(ITEM)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(IN) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
INTEGER :: ITEM(VEC%D)
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for bad input
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! get item
ITEM = VEC%DAT(:,IND)
RETURN
END FUNCTION INTVECTORITEM

! memory access function (all data)
FUNCTION INTVECTORDATA(VEC) RESULT(DAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(IN) :: VEC
INTEGER :: DAT(VEC%D, VEC%N)
! get data
DAT(:,:) = VEC%DAT(:,1:VEC%N)
RETURN
END FUNCTION INTVECTORDATA

! free memory
SUBROUTINE INTVECTORFREE(VEC, ISTAT)
IMPLICIT NONE
! parameters
CLASS(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local variables
INTEGER :: LSTAT
! free data
IF (ALLOCATED(VEC%DAT)) THEN
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF(PRESENT(ISTAT)) ISTAT=LSTAT
END IF
RETURN
END SUBROUTINE INTVECTORFREE

! --- Functions for RealVector type --- !

! constructor
FUNCTION NEWREALVECTOR(DIM, ISTAT) RESULT(VEC)
IMPLICIT NONE
! parameters
INTEGER, OPTIONAL, INTENT(IN) :: DIM
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
TYPE(REALVECTOR) :: VEC
! set information
VEC%N = 0
VEC%MAXLEN = 0
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for optionals
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
END FUNCTION NEWREALVECTOR

! push function
SUBROUTINE REALVECTORPUSH(VEC, ITEM, ISTAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(INOUT) :: VEC
REAL(KIND=R8), INTENT(IN) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
REAL(KIND=R8) :: TMP(VEC%D,VEC%N)
INTEGER :: LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for illegal size
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! no resizing required
IF(VEC%N < VEC%MAXLEN) THEN
   ! increment counter and add
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
! needs to be initialized
ELSE IF(VEC%MAXLEN .EQ. 0) THEN
   ! allocate
   ALLOCATE(VEC%DAT(VEC%D,8), STAT=LSTAT)
   IF(LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 8
   ! increment counter and add
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
! resize data
ELSE
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! reallocate
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
   VEC%MAXLEN=VEC%MAXLEN*2
   ! restore current data
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
   ! add new data
   VEC%N = VEC%N + 1
   VEC%DAT(:,VEC%N) = ITEM(:)
END IF
RETURN
END SUBROUTINE REALVECTORPUSH

! pop function
SUBROUTINE REALVECTORPOP(VEC, ITEM, ISTAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(INOUT) :: VEC
REAL(KIND=R8), INTENT(OUT) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local variables
REAL(KIND=R8) :: TMP(VEC%D, VEC%N-1)
INTEGER :: LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for illegal size
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! check for illegal length
IF(VEC%N .EQ. 0) THEN
   ITEM = 0
   RETURN
END IF
! get item
ITEM(:) = VEC%DAT(:,VEC%N)
VEC%N = VEC%N - 1
! reclaim partial memory
IF (VEC%N < VEC%MAXLEN / 2 .AND. &
  VEC%N > 6) THEN
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ! reallocate
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! copy back
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
! reclaim all memory
ELSE IF (VEC%N .EQ. 0) THEN
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE REALVECTORPOP

! insert function
SUBROUTINE REALVECTORINS(VEC, ITEM, IND, ISTAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(INOUT) :: VEC
REAL(KIND=R8), INTENT(IN) :: ITEM(:)
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
REAL(KIND=R8) :: TMP(VEC%D,VEC%N)
INTEGER :: LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for illegal size
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! check for bad index
IF(IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -2
   RETURN
END IF
! no resizing required
IF(VEC%N < VEC%MAXLEN) THEN
   ! increment counter and insert
   VEC%DAT(:,IND+1:VEC%N+1) = &
      VEC%DAT(:,IND:VEC%N)
   VEC%DAT(:,IND) = ITEM(:)
   VEC%N = VEC%N + 1
! resize data
ELSE
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! reallocate
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ALLOCATE(VEC%DAT(VEC%D,VEC%MAXLEN*2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN=VEC%MAXLEN*2
   ! restore current data
   VEC%DAT(:,1:IND-1) = TMP(:,1:IND-1)
   VEC%DAT(:,IND+1:VEC%N+1) = TMP(:,IND:VEC%N)
   ! insert new data
   VEC%N = VEC%N + 1
   VEC%DAT(:,IND) = ITEM(:)
END IF
RETURN
END SUBROUTINE RealVectorIns

! delete function
SUBROUTINE REALVECTORDEL(VEC, IND, ISTAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
REAL(KIND=R8) :: TMP(VEC%D,VEC%N-1)
INTEGER :: LSTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for bad input
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
! shift other elements down
VEC%DAT(:,IND:VEC%N-1) = VEC%DAT(:,IND+1:VEC%N)
VEC%N = VEC%N - 1
! reclaim partial memory
IF (VEC%N < VEC%MAXLEN / 2 .AND. &
  VEC%N > 6) THEN
   ! save current data
   TMP(:,:) = VEC%DAT(:,1:VEC%N)
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   ! reallocate
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! copy back
   VEC%DAT(:,1:VEC%N) = TMP(:,:)
! reclaim all memory
ELSE IF (VEC%N .EQ. 0) THEN
   ! free
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF (LSTAT .NE. 0) THEN
      IF(PRESENT(ISTAT)) ISTAT = LSTAT
      RETURN
   END IF
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE REALVECTORDEL

! set the value of VEC at IND to ITEM
SUBROUTINE REALVECTORSET(VEC, IND, ITEM, ISTAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
REAL(KIND=R8), INTENT(IN) :: ITEM(:)
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for bad input
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF
IF (VEC%D .NE. SIZE(ITEM)) THEN
   IF(PRESENT(ISTAT)) ISTAT = -2
   RETURN
END IF
! set VEC%DAT(:,IND) to ITEM
VEC%DAT(:,IND) = ITEM(:)
RETURN
END SUBROUTINE REALVECTORSET

! memory access function (single item)
FUNCTION REALVECTORITEM(VEC, IND, ISTAT) RESULT(ITEM)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(IN) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
REAL(KIND=R8) :: ITEM(VEC%D)
! initialize ISTAT
IF(PRESENT(ISTAT)) ISTAT = 0
! check for bad input
IF (IND < 1 .OR. IND > VEC%N) THEN
   IF(PRESENT(ISTAT)) ISTAT = -1
   RETURN
END IF 
! get item
ITEM = VEC%DAT(:,IND)
RETURN
END FUNCTION REALVECTORITEM

! memory access function (all data)
FUNCTION REALVECTORDATA(VEC) RESULT(DAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(IN) :: VEC
REAL(KIND=R8) :: DAT(VEC%D, VEC%N)
! retrieve data
DAT(:,:) = VEC%DAT(:,1:VEC%N)
RETURN
END FUNCTION REALVECTORDATA

! free memory
SUBROUTINE REALVECTORFREE(VEC, ISTAT)
IMPLICIT NONE
! parameters
CLASS(REALVECTOR), INTENT(INOUT) :: VEC
INTEGER, OPTIONAL, INTENT(OUT) :: ISTAT
! local data
INTEGER :: LSTAT
! free data
IF (ALLOCATED(VEC%DAT)) THEN
   DEALLOCATE(VEC%DAT, STAT=LSTAT)
   IF(PRESENT(ISTAT)) ISTAT = LSTAT
END IF
RETURN
END SUBROUTINE REALVECTORFREE

END MODULE VECTOR_CLASS
