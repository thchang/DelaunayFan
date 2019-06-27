MODULE VECTOR
! The VECTOR module defines the INTVECTOR type for dynamically expanding
! multidimensional arrays. The INTVECTOR type has 4 fields, the leading
! dimension D, the current length LENGTH, the allocated length MAXLEN, and
! the data DAT(D,:). The INTVECTOR type manages memory similarly to the
! vector type from C++. INTVECTOR starts as an array of size D X 8 and
! doubles in size whenever more memory is needed, or halves in size when
! less memory is needed. The INTVECTOR interface is modeled after a stack
! data structure and features 4 functions/subroutines for accessing data:
!
! INTVECTORPUSH pushes an array of length D to the front of the INTVECTOR
! stack.
!
! INTVECTORPOP pops an array of length D off the top of the INTVECTOR stack.
! 
! INTVECTORDEL deletes an array of length D from a specified index in the
! middle of the stack.
!
! INTVECTORCONTAINS searches the stack for any occurrence of a specific
! D-dimensional vector and returns the index of that vector if found.
!
! The INTVECTOR module also features a constructor NEWINTVECTOR for
! initializing a new INTVECTOR type with default values, and a
! deconstructor INTVECTORFREE for deallocating all memory associated
! with the INTVECTOR type.
!
! Primary Author: Tyler H. Chang
! Last Update: October, 2017
IMPLICIT NONE

! Definition of the INTVECTOR type.
TYPE INTVECTOR
   ! Meta data for the INTVECTOR type.
   INTEGER :: D, LENGTH, MAXLEN
   ! Dynamic array with leading dimension D.
   INTEGER, ALLOCATABLE :: DAT(:,:)
END TYPE INTVECTOR

CONTAINS ! Interface functions and subroutines.

SUBROUTINE NEWINTVECTOR(VEC, IERR, DIM)
! The INTVECTOR constructor creates a new INTVECTOR type.
!
! Output:
!
! On output, VEC is an initialized INTVECTOR object ready to be used.
!
! Optional arguments:
!
! If present, IERR returns an error code associated with execution.
! Succesful initialization = 0.
!
! If present, DIM contains the dimension of the INTVECTOR on input. If
! omitted, the dimension is assumed to be 1.

! Parameter list.
TYPE(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, OPTIONAL, INTENT(OUT) :: IERR
INTEGER, OPTIONAL, INTENT(IN) :: DIM

! Initialize meta data.
VEC%LENGTH = 0
VEC%MAXLEN = 0
! Check for optional arguments.
IF(PRESENT(DIM)) THEN
   IF (DIM < 1) THEN
      IERR = -3
      RETURN
   END IF
   IERR = 0
   VEC%D = DIM
ELSE
   IERR = 0
   VEC%D = 1
END IF
RETURN
END SUBROUTINE NEWINTVECTOR

SUBROUTINE INTVECTORPUSH(VEC, ITEM, IERR)
! The INTVECTORPUSH subroutine adds a new item to the end of the INTVECTOR
! data stack.
!
! Required arguments:
!
! On input, VEC is an initialized INTVECTOR object with leading dimension D,
! and ITEM is an integer array of length D.
!
! On output, ITEM has been appended to VEC.
!
! IERR returns an error code associated with execution.
! Succesful push operation = 0.

! Parameter list.
TYPE(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: ITEM(:)
INTEGER, INTENT(OUT) :: IERR

! Local array in case of reallocation.
INTEGER :: TMP(VEC%D,VEC%LENGTH)

! Check for an illegal size.
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IERR = -2
   RETURN
END IF
IF(VEC%LENGTH < VEC%MAXLEN) THEN ! No resizing required.
   ! Increment counter and append ITEM to VEC%DAT.
   VEC%LENGTH = VEC%LENGTH + 1
   VEC%DAT(:,VEC%LENGTH) = ITEM(:)
ELSE IF(VEC%MAXLEN .EQ. 0) THEN ! VEC%DAT needs to be initialized.
   ! Allocate VEC%DAT.
   ALLOCATE(VEC%DAT(VEC%D,8), STAT=IERR)
   IF(IERR .NE. 0) RETURN
   VEC%MAXLEN = 8
   ! Increment counter and append ITEM to VEC%DAT.
   VEC%LENGTH = VEC%LENGTH + 1
   VEC%DAT(:,VEC%LENGTH) = ITEM(:)
ELSE ! VEC%DAT needs to be resized.
   ! Save the current data in TMP.
   TMP(:,:) = VEC%DAT(:,1:VEC%LENGTH)
   ! Reallocate VEC%DAT.
   DEALLOCATE(VEC%DAT, STAT=IERR)
   IF(IERR /= 0) RETURN
   ALLOCATE(VEC%DAT(VEC%D,VEC%MAXLEN*2), &
      STAT=IERR)
   IF(IERR /= 0) RETURN
   VEC%MAXLEN=VEC%MAXLEN*2
   ! Restore contents to VEC%DAT.
   VEC%DAT(:,1:VEC%LENGTH) = TMP(:,:)
   ! Append ITEM.
   VEC%LENGTH = VEC%LENGTH + 1
   VEC%DAT(:,VEC%LENGTH) = ITEM(:)
END IF
RETURN
END SUBROUTINE INTVECTORPUSH

SUBROUTINE INTVECTORPOP(VEC, ITEM, IERR)
! The INTVECTORPOP subroutine pops an item off the end of the INTVECTOR
! data stack.
!
! Required arguments:
!
! On input, VEC is an initialized INTVECTOR object with leading dimension D,
! and ITEM is an empty integer array of length D.
!
! On output, VEC%LENGTH has decreased by 1 and the removed item is contained
! in the ITEM array.
!
! IERR returns an error code associated with execution.
! Succesful pop operation = 0.

! Parameter list.
TYPE(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(OUT) :: ITEM(:), IERR

! Local array in case of reallocation.
INTEGER :: tmp(vec%d, vec%length-1)

! Check for an illegal leading dimension.
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IERR = -2
   RETURN
END IF
! Check for an illegal length.
IF(VEC%LENGTH .EQ. 0) THEN
   IERR = 0
   ITEM = 0
   RETURN
END IF
! Get ITEM from VEC%DAT.
ITEM(:) = VEC%DAT(:,VEC%LENGTH)
VEC%LENGTH = VEC%LENGTH - 1
! Reclaim memory when possible.
IF (VEC%LENGTH < VEC%MAXLEN / 2 .AND. VEC%LENGTH > 6) THEN
   ! Make a copy of the current data.
   TMP(:,:) = VEC%DAT(:,1:VEC%LENGTH)
   ! Reallocate the memory.
   DEALLOCATE(VEC%DAT, STAT=IERR)
   IF (IERR .NE. 0) RETURN
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=IERR)
   IF (IERR .NE. 0) RETURN
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! Restore current values.
   VEC%DAT(:,1:VEC%LENGTH) = TMP(:,:)
! Reclaim all memory if VEC%LENGTH was 1.
ELSE IF (VEC%LENGTH .EQ. 0) THEN
   ! Free memory and don't reallocate.
   DEALLOCATE(VEC%DAT, STAT=IERR)
   IF (IERR /= 0) RETURN
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE INTVECTORPOP

SUBROUTINE INTVECTORDEL(VEC, IND, IERR)
! The INTVECTORDEL subroutine deletes an item from anywhere in the INTVECTOR
! data stack.
!
! Required arguments:
!
! On input, VEC is an initialized INTVECTOR object with leading dimension D,
! and 0 < IND < VEC%LENGTH+1 is the index of an item in VEC%DAT to delete.
!
! On output, VEC%DAT(:,IND) has been deleted.
!
! IERR returns an error code associated with execution.
! Succesful delete operation = 0.

! Parameter list.
TYPE(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: IND
INTEGER, INTENT(OUT) :: IERR

! Local array in case of reallocation.
INTEGER :: TMP(VEC%D,VEC%LENGTH-1)

! Check for bad inputs.
IF (IND < 1 .OR. IND > VEC%LENGTH) THEN
   IERR = -2
   RETURN
END IF
IERR = 0
! Shift other elements down.
VEC%DAT(:,IND:VEC%LENGTH-1) = VEC%DAT(:,IND+1:VEC%LENGTH)
VEC%LENGTH = VEC%LENGTH - 1
! Reclaim partial memory.
IF (VEC%LENGTH < VEC%MAXLEN / 2 .AND. VEC%LENGTH > 6) THEN
   ! Make a copy of the current data.
   TMP(:,:) = VEC%DAT(:,1:VEC%LENGTH)
   ! Reallocate VEC%DAT.
   DEALLOCATE(VEC%DAT, STAT=IERR)
   IF (IERR .NE. 0) RETURN
   ALLOCATE(VEC%DAT(VEC%D, VEC%MAXLEN / 2), STAT=IERR)
   IF (IERR .NE. 0) RETURN
   VEC%MAXLEN = VEC%MAXLEN / 2
   ! Restore current data.
   VEC%DAT(:,1:VEC%LENGTH) = TMP(:,:)
! Reclaim all memory.
ELSE IF (VEC%LENGTH .EQ. 0) THEN
   DEALLOCATE(VEC%DAT, STAT=IERR)
   IF (IERR .NE. 0) RETURN
   VEC%MAXLEN = 0
END IF
RETURN
END SUBROUTINE INTVECTORDEL

FUNCTION INTVECTORCONTAINS(VEC, ITEM, IERR) RESULT(IND)
! The INTVECTORCONTAINS subroutine searches the contents of the INTVECTOR
! data stack for the contents of ITEM.
!
! Required arguments:
!
! On input, VEC is an initialized INTVECTOR object with leading dimension D,
! and ITEM(1:D) is the integer array to search for.
!
! If found ITEM is found, INTVECTORCONTAINS returns the index. Otherwise, it
! returns 0.
!
! IERR returns an error code associated with execution.
! Succesful search operation = 0.

! Parameter list.
TYPE(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(IN) :: ITEM(:)
INTEGER, INTENT(OUT) :: IERR
! Return value.
INTEGER :: IND
! Local variables.
INTEGER :: I

! Check for an illegal size.
IF(VEC%D .NE. SIZE(ITEM)) THEN
   IERR = -2
   RETURN
END IF
! Search an unsorted list.
IND = 0
DO I=1,VEC%LENGTH
   IF(ALL(VEC%DAT(1:VEC%D-1,I) == ITEM(:))) THEN
      IND = I
      EXIT
   END IF
END DO
RETURN
END FUNCTION INTVECTORCONTAINS

SUBROUTINE INTVECTORFREE(VEC, IERR)
! The INTVECTOR deconstructor frees all memory associated with VEC.
! The error flag IERR returns 0 when memory deallocation was successful.

! Parameter list.
TYPE(INTVECTOR), INTENT(INOUT) :: VEC
INTEGER, INTENT(OUT) :: IERR
! Free VEC%DAT.
IF (ALLOCATED(VEC%DAT)) THEN
   DEALLOCATE(VEC%DAT, STAT=IERR)
END IF
RETURN
END SUBROUTINE INTVECTORFREE

END MODULE VECTOR
