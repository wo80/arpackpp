/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNSymGRe.cc.
   Example program that illustrates how to solve a real nonsymmetric
   generalized eigenvalue problem in regular mode using the
   ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular
      mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                        (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.
      {nnzA, irowA, pcolA, valA}: matrix B data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by ARluNonSymGenEig to solve
      some linear systems involving B.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lnmatrxc.h       StiffnessMatrix, a function that generates
                       matrix A in CSC format.
      lnmatrxd.h       MassMatrix, a function tha generates matrix
                       B in CSC format.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlgnsym.h       The ARluNonSymGenEig class definition.
      lnsymsol.h       The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lnmatrxc.h"
#include "lnmatrxd.h"
#include "arlnsmat.h"
#include "arlgnsym.h"
#include "lnsymsol.h"

int main()
{

  // Defining variables;

  int     n;               // Dimension of the problem.
  int     nnza,   nnzb;    // Number of nonzero elements in A and B.
  int     *irowa, *irowb;  // pointers to arrays that store the row
                           // indices of the nonzeros in A and B.
  int     *pcola, *pcolb;  // pointers to arrays of pointers to the
                           // beginning of each column of A and B in
                           // valA and ValB.
  double  rho;             // A parameter used in StiffnessMatrix.
  double  *valA,  *valB;   // pointers to arrays that store the
                           // nonzero elements of A and B.

  int nev = 4; // Number of requested eigenvalues.

  // Creating matrices A and B.

  n  = 100;
  rho = 10.0;
  StiffnessMatrix(n, rho, nnza, valA, irowa, pcola);
  ARluNonSymMatrix<double, double> A(n, nnza, valA, irowa, pcola);

  MassMatrix(n, nnzb, valB, irowb, pcolb);
  ARluNonSymMatrix<double, double> B(n, nnzb, valB, irowb, pcolb);

  // Defining what we need: the four eigenvectors with largest magnitude.

  ARluNonSymGenEig<double> dprob(nev, A, B);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

  int nconv = dprob.ConvergedEigenvalues();
  
  return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
} // main

