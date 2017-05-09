      subroutine init2d(dm, b, ierr)
      implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>
#include <finclude/petscvec.h90>

      DM dm
      Vec b
      PetscErrorCode ierr
      PetscScalar hx, hy
      PetscScalar,pointer :: pb(:,:)
      DMDALocalInfo info
      PetscInt i, j
      PetscInt mx, my, xm, ym, xs, ys
      
      call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,  &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,ierr)

      call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,            &
     &     PETSC_NULL_INTEGER,ierr)

      call DMDAVecGetArrayF90(dm, b, pb, ierr)
      
      hx = 1./(mx-1);
      hy = 1./(mx-1);
  
      do i=xs, xs+xm-1
         do j=ys, ys+ym-1
            if (i.eq.0 .or. j.eq.0 .or. i.eq.mx-1 .or. j.eq.my-1) then
               pb(i, j) = 0.
            else
               pb(i, j) = 1.
            end if
         end do
      end do

      call DMDAVecRestoreArrayF90(dm, b, pb, ierr)
      end 
