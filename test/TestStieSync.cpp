/*
This is the test file for the SE-Sync problem defined in StieSync.cpp and StieSync.h 
*/

#ifndef TESTSTIESYNC_CPP
#define TESTSTIESYNC_CPP

/*Output to console*/
#include <iostream>
/*Generate random number*/
#include "Others/randgen.h"
/*Computational time*/
#include <ctime>

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/StieSync/StieSync.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"

/*Trust-region based solvers*/
#include "Solvers/RTRNewton.h"

/*The global head file*/
#include "Others/def.h"

using namespace ROPTLIB;

/*If the file TESTSTIESYNC is defined in def.h file, then the following
main() function will be called. */
#if defined(TESTSTIESYNC)

void testStieSync();

int main(void)
{
	testStieSync();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

void testStieSync()
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	genrandseed(tt);

	// size of the Stiefel manifold (read from input required)
	integer n = 12, p = 8, d = 8, ND = n*d;	

	//load Q matrix from dataset required
	double *Q = new double[ND * ND];
	for (integer i = 0; i < ND; i++)
	{
		for (integer j = i; j < ND; j++)
		{
			Q[i + j * ND] = genrandnormal();
			Q[j + i * ND] = Q[i + j * ND];
		}
	}
	
	//initialize a point on the manifold
	integer numofmanis = 1;
	integer numofmani1 = n;
	StieVariable StieX(p,d);
	ProductElement ProdX(numofmanis , &StieX, numofmani1);
	//if() path is provided
	//else
	ProdX.RandInManifold();

	//Define the manifold domain
	Stiefel Mani(p,d);
	ProductManifold Domain(numofmanis, &Mani, numofmani1);

	//Define the SE-Sync problem
	StieSync Prob(Q,p,d,n);
	
	// Set the problem domain
	Prob.SetDomain(&Domain);

	// output the parameters of the manifold of domain
	Domain.CheckParams();

	// test RTRNewton
	printf("********************************Check RTRNewton*************************************\n");
	RTRNewton RTRNewtonsolver(&Prob, &ProdX);
	RTRNewtonsolver.Debug = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&ProdX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);

	delete[] Q;

	return;	
}

#endif // end of TESTSTIESYNC
#endif // end of TESTSTIESYNC_CPP
