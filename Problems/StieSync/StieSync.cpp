#include "Problems/StieSync/StieSync.h"

namespace ROPTLIB{

	StieSync::StieSync(double *inQ, integer inp, integer ind, integer inn)
	{
		Q = inQ;
		p = inp;
		d = ind;
		n = inn;
	};

	StieSync::~StieSync(void)
	{
	};

	double StieSync::f(Variable *x) const
	{
		integer N = n, P = p, D = d, ND = n*d;
		const double *xxM = x->ObtainReadData();		
		// xxM \in St(d,p)^n, belongs to R_{p x nd}		
		Vector *xQ = x->ConstructEmpty();
		SharedSpace *Temp = new SharedSpace(xQ);
		double *temp = xQ->ObtainWriteEntireData();
		double result = 0;
		
		//temp = MxxM * MQ
		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		dgemm_(transn, transn, &P, &ND, &ND, &one, const_cast<double *> (xxM), &P, Q, &ND, &zero, temp, &P);
		
		integer length = p*d*n;
		//output = xxM(:)^T * temp(:) = trace(xxM^T * temp) = trace(xxM^T * xxM * Q) ... Eq.(13) in SE-Sync Conference Paper
		result = ddot_(&length, const_cast<double *> (xxM), &GLOBAL::IONE, temp, &GLOBAL::IONE);
		if (UseGrad)
		{
			x->AddToTempData("xQ", Temp);
		}
		else
		{
			delete Temp;
		}
		return (result/2);
	};

	//Equations for Gradient and Hessian are obtained from Eq.(14) of the SE-Sync Conference Paper

	void StieSync::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Temp = x->ObtainReadTempData("xQ");
		Vector *xQ = Temp->GetSharedElement();
		Domain->ScaleTimesVector(x, 1.0, xQ, egf);
	};

	void StieSync::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const double *etaxTV = etax->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		integer N = n, P = p, D = d, ND = n*d;
		double one = 1, zero = 0;
		// exxiTV <- etaxTV * Q, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &ND, &ND, &one, const_cast<double *> (etaxTV), &P, Q, &ND, &zero, exixTV, &P);
		Domain->ScaleTimesVector(x, 1.0, exix, exix);
	};
}; /*end of ROPTLIB namespace*/
