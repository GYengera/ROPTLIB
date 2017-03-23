/*
This defines the class for the SE-Sync final optimization,i.e. after the rank-restricted semidefinite relaxation.
*/

#ifndef STIESYNC_H
#define STIESYNC_H

#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

namespace ROPTLIB{

	class StieSync : public Problem{
	public:
		StieSync(double *inQ, integer inp, integer ind, integer inn);
		virtual ~StieSync();
		virtual double f(Variable *x) const;

		//virtual void RieGrad() const;
		//virtual void RieHessianEta() const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *Q;
		integer p;
		integer d;
		integer n;
	};
};/*end of ROPTLIB namespace*/
#endif
