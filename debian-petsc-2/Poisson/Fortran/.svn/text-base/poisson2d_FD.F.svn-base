      module PETScMatContext
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

      DM :: dmContext
      end module PETScMatContext

      subroutine createMat(dm, A, kind)
      use PETScMatContext, only : dmContext
      implicit none 

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

      PetscErrorCode ierr
      PetscBool kind
      Mat A
      DM dm
      PetscInt mx, my, xm, ym
      PetscInt localsize, totalsize
      external PoissonMatMult2d
      
      if (kind) then
         call DMCreateMatrix(dm, MATAIJ, A, ierr)
         call poissonFD2d_petsc(dm, A, ierr)
      else
         call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER   &
     &               ,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,               &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
     &               PETSC_NULL_INTEGER,ierr)

         call DMDAGetCorners(dm,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,     &
     &               PETSC_NULL_INTEGER,xm,ym,                             &
     &               PETSC_NULL_INTEGER, ierr)
         localsize = xm*ym
         totalsize = mx*my

         dmContext = dm

         call MatCreateShell(PETSC_COMM_WORLD, localsize, localsize,       &
     &        totalsize, totalsize, PETSC_NULL_OBJECT, A, ierr)
         call MatShellSetOperation(A, MATOP_MULT, PoissonMatMult2d,ierr)
      end if

      call MatSetFromOptions(A, ierr)

      end

      subroutine PoissonMatMult2d(A, x, y, ierr)
      use PETScMatContext, only : dmContext
      implicit none
#include <finclude/petscmat.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

      PetscErrorCode ierr
      Mat A
      Vec x, y
      DM dm

      !call MatShellGetContext(A, dm, ierr)
      call poissonFD2dMatVec(x, y, dmContext, ierr)

      end

      subroutine poissonFD2dMatVec(x, y, dm, ierr)
      implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

      PetscErrorCode ierr
      DM dm
      PetscScalar hx2, hy2
      PetscScalar cx, cy, cd
      PetscInt mx, my, xm, ym, xs, ys
      PetscInt i, j
      Vec x, y
      Vec xLocal, yLocal
      PetscScalar ,pointer :: px(:,:), py(:,:)

      call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,  &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,ierr)

      call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,            &
     &     PETSC_NULL_INTEGER,ierr)
      
      hx2 = (mx-1)*(mx-1)
      hy2 = (my-1)*(my-1)
      cx=-hx2
      cy=-hy2 
      cd=2.*(hx2+hy2)

      call DMGetLocalVector(dm, xLocal, ierr)
      call DMGlobalToLocalBegin(dm, x, INSERT_VALUES, xLocal, ierr)
      call DMGlobalToLocalEnd(dm, x, INSERT_VALUES, xLocal, ierr)

      call DMGetLocalVector(dm, yLocal, ierr)
      call VecSet(yLocal, 0., ierr) 
      call VecSet(y, 0., ierr)

      call DMDAVecGetArrayF90(dm, xLocal, px, ierr)
      call DMDAVecGetArrayF90(dm, yLocal, py, ierr)

      do i=xs, xs+xm-1
         do j=ys, ys+ym-1
            if (i.eq.0 .or. j.eq.0 .or. i.eq.mx-1 .or. j.eq.my-1) then
               py(i, j) = px(i, j)
            else
               py(i, j) = cd*px(i, j) + cx*(px(i+1, j) + px(i-1, j))    &
     &              + cy*(px(i, j+1) + px(i, j-1))
            end if
         end do
      end do
  
      call DMDAVecRestoreArrayF90(dm, xLocal, px, ierr)
      call DMDAVecRestoreArrayF90(dm, yLocal, py, ierr)
      call DMRestoreLocalVector(dm, xLocal, ierr)
      call DMRestoreLocalVector(dm, yLocal, ierr)

      call DMLocalToGlobalBegin(dm, yLocal, ADD_VALUES, y, ierr)
      call DMLocalToGlobalEnd(dm, yLocal, ADD_VALUES, y, ierr)
      call PetscLogFlops(7.*mx*my, ierr)
      end


      subroutine poissonFD2d_petsc(dm, A, ierr)
      implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

      PetscErrorCode ierr
      Mat A
      DM dm
      PetscScalar hx2, hy2, coef, coef5(5)
      PetscInt mx, my, xm, ym, xs, ys
      PetscInt i, j
      PetscInt indxBeg, indxEnd
      PetscInt indyBeg, indyEnd
      MatStencil row(4), col(4, 5)

      call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,  &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,ierr)

      call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,            &
     &     PETSC_NULL_INTEGER,ierr)
      
      hx2 = 1./((mx-1)*(mx-1))
      hy2 = 1./((my-1)*(my-1))

      coef = 1.
      coef5 = (/ 2./hx2 + 2./hy2, -1./hx2, -1./hx2, -1./hy2, -1./hy2 /)

      do i=xs, xs+xm-1
         do j=ys, ys+ym-1
            row(MatStencil_i) = i
            row(MatStencil_j) = j
            if (i.eq.0 .or. j.eq.0 .or. i.eq.mx-1 .or. j.eq.my-1) then
               coef = 1.
               call MatSetValuesStencil(A, 1, row, 1, row, coef,           &
     &              INSERT_VALUES, ierr)
            else
               col(MatStencil_i, 1) = i
               col(MatStencil_j, 1) = j
               col(MatStencil_i, 2) = i-1
               col(MatStencil_j, 2) = j
               col(MatStencil_i, 3) = i+1
               col(MatStencil_j, 3) = j
               col(MatStencil_i, 4) = i
               col(MatStencil_j, 4) = j-1
               col(MatStencil_i, 5) = i
               col(MatStencil_j, 5) = j+1
               call MatSetValuesStencil(A, 1, row, 5, col, coef5,           &
     &              INSERT_VALUES, ierr)
            end if
         end do
      end do


      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY, ierr)

      end

