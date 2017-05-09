      program Poisson2D
      implicit none

#include <finclude/petsc.h>
      PetscErrorCode ierr
      PetscInt nx, ny 
      PetscInt stage
      PetscBool flg
      Mat A
      Vec u, b
      KSP solver
      PC pc
      DM dm

      nx = 257
      ny = 257

      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
           CHKERRQ(ierr)    
 
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER, "-nx", nx, 
     &          PETSC_NULL_CHARACTER, ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER, "-ny", ny, 
     &          PETSC_NULL_CHARACTER, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-assemble", flg
     &     , ierr)

      call PetscLogStageRegister("Domain creation", stage, ierr)
      call PetscLogStagePush(stage, ierr)
      call createDomain(dm, nx, ny, ierr)
      call PetscLogStagePop(ierr)

      call PetscLogStageRegister("Matrix creation", stage, ierr)
      call PetscLogStagePush(stage, ierr)
      call createMat(dm, A, flg, ierr)
      call PetscLogStagePop(ierr)

      call PetscLogStageRegister("Second member creation", stage, ierr)
      call PetscLogStagePush(stage, ierr)
      call DMCreateGlobalVector(dm, b, ierr)
      call VecDuplicate(b, u, ierr)
      call init2d(dm, b, ierr)
      call PetscLogStagePop(ierr)
      
      call PetscLogStageRegister("KSP creation", stage, ierr)
      call PetscLogStagePush(stage, ierr)
      call KSPCreate(PETSC_COMM_WORLD, solver, ierr)
      call KSPSetOptionsPrefix(solver, "poisson_", ierr)
      call KSPSetOperators(solver, A, A, DIFFERENT_NONZERO_PATTERN,ierr)
      call KSPSetType(solver, KSPCG, ierr)
      call KSPGetPC(solver, pc, ierr)
      call PCSetType(pc, PCNONE, ierr)
      call KSPSetFromOptions(solver, ierr)
      call PetscLogStagePop(ierr)
  
      call PetscLogStageRegister("Solver", stage, ierr)
      call PetscLogStagePush(stage, ierr)
      call KSPSolve(solver, b, u, ierr)
      call PetscLogStagePop(ierr)

      call VecDestroy(u, ierr);
      call VecDestroy(b, ierr);
      call MatDestroy(A, ierr);
      call DMDestroy(dm, ierr);
      call KSPDestroy(solver, ierr);
      call PetscFinalize(ierr)
      end
