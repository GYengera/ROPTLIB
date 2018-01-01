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
/*Eigen library, however not used*/
#include <Eigen/Sparse>
/*Read and Write to Files*/
#include <fstream>

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

void testStieSync(char* qpath, int inn, int ind, int inp, char* x0path = NULL);

int main(int argc, char* argv[])
{
	if (argc == 5)
	{
		testStieSync(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
	}
	else if (argc == 6)
	{
		testStieSync(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[5]);
	}
	else
	{
		std::cout << "Incorrect input arguments" << std::endl;
		return 0;
	}
#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;	
}

void testStieSync(char* qpath, int inn, int ind, int inp, char* x0path)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	genrandseed(tt);

	// size of the Stiefel manifold (read from input)
	integer n = inn, p = inp, d = ind, ND = n*d;

	//load Q matrix from dataset
	std::ifstream input(qpath,std::ios::in);
	if(!input)
	{
		std::cout << "ERROR: Q matrix data not defined" << std::endl;
		return;
	}
	
	double *Q = new double[ND*ND]; 
	integer row,col;
	double value;
	std::string line;
	for (row = 0; row < ND; row++)
	{
		for (col = row; col < ND; col++)
		{
			Q[row + col * ND] = 0;
			Q[col + row * ND] = Q[row + col * ND];
		}
	}

	while(std::getline(input,line))
	{
		std::stringstream entry(line);	
		entry >> row >> col >> value;
		Q[(row-1) + (col-1)*ND] = value;
	}
	input.close();

	/* The code to store Q in an Eigen::SparseMatrix object, but not applicable in ROPTLIB yet.
	typedef Eigen::Triplet<double> T;
	std::vector<T> elements;
	integer row, col;
	double value;
	std::string line;

	while(std::getline(input,line))
	{
		std::stringstream entry(line);	
		entry >> row >> col >> value;
		elements.push_back(T(row-1,col-1,value));
	}

	Eigen::SparseMatrix<double> Q(ND, ND);
	Q.setFromTriplets(elements.begin(), elements.end());*/

	
	//initialize a point on the manifold
	integer numofmanis = 1;
	integer numofmani1 = n;
	StieVariable StieX(p,d);
	ProductElement ProdX(numofmanis, &StieX, numofmani1);
	if(x0path == NULL)
	{
		ProdX.RandInManifold();
	}
	else
	{
		std::ifstream input(qpath,std::ios::in);
		if(!input)
		{
			std::cout << "ERROR: X matrix could not be initialized" << std::endl;
			return;
		}

		double *x0 = ProdX.ObtainWriteEntireData();
		while(std::getline(input,line))
		{
			std::stringstream entry(line);	
			entry >> row >> col >> value;
			x0[(row-1) + (col-1)*p] = value;
		}
	}

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
	
	//Store xopt into a file
	const double *optimizer = xopt->ObtainReadData();
	std::ofstream outfile;
	//Destination to store xopt can be modified
	outfile.open("~/xopt.txt");
	if(!outfile) {
    		std::cout << "Cannot open file to save xopt";
	}
	for(row = 0; row < p; ++row){
		for(col = 0; col < ND; ++col){
			outfile << (row+1) << " " << (col+1) << " " << optimizer[row + col * p];	
			outfile << std::endl;
	}
	}
	outfile.close();

	delete[] Q;

	return;	
}
#endif // end of TESTSTIESYNC_CPP
