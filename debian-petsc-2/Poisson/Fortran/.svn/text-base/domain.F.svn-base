      subroutine createDomain(dm, nx, ny, ierr)
      implicit none

#include <finclude/petsc.h>

      PetscErrorCode ierr
      PetscInt nx, ny
      DM dm

      call DMDACreate2d(PETSC_COMM_WORLD, 
     &     DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, 
     &     DMDA_STENCIL_STAR,
     &     nx, ny, PETSC_DECIDE, PETSC_DECIDE,
     &     1, 1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, dm, ierr)

      end
