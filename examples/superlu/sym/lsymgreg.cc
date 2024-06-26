/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSymGReg.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in regular mode using the
   ARluSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular mode,
      where A and B are obtained from the finite element discretization
      of the 1-dimensional discrete Laplacian
                                  d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: lower triangular part of matrix A 
                                  stored in CSC format.
      {nnzB, irowB, pcolB, valB}: lower triangular part of matrix B 
                                  stored in CSC format.

   3) Library called by this example:

      The SuperLU package is called by ARluSymGenEig to solve
      some linear systems involving B.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lsmatrxc.h       SymmetricMatrixC, a function that generates
                       matrix A in CSC format.
      lsmatrxd.h       SymmetricMatrixD, a function that generates
                       matrix B in CSC format.
      arlsmat.h        The ARluSymMatrix class definition.
      arlgsym.h        The ARluSymGenEig class definition.
      lsymsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lsmatrxc.h"
#include "lsmatrxd.h"
#include "arlsmat.h"
#include "arlgsym.h"
#include "lsymsol.h"


int main()
{

  // Defining variables;

  int    n;              // Dimension of the problem.
  int    nnzA,   nnzB;   // Number of nonzero elements in A and B.
  int    *irowA, *irowB; // pointer to an array that stores the row
                         // indices of the nonzeros in A and B.
  int    *pcolA, *pcolB; // pointer to an array of pointers to the
                         // beginning of each column of A (B) in valA (valB).
  double *valA,  *valB;  // pointer to an array that stores the nonzero
                         // elements of A and B.

  int nev = 4; // Number of requested eigenvalues.

  // Creating matrices A and B.

  n = 100;
  SymmetricMatrixC(n, nnzA, valA, irowA, pcolA);
  ARluSymMatrix<double> A(n, nnzA, valA, irowA, pcolA);

  SymmetricMatrixD(n, nnzB, valB, irowB, pcolB);
  ARluSymMatrix<double> B(n, nnzB, valB, irowB, pcolB);

  // Defining what we need: the four eigenvectors with largest magnitude.

  ARluSymGenEig<double> dprob(nev, A, B);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

  int nconv = dprob.ConvergedEigenvalues();
  
  return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
} // main.

