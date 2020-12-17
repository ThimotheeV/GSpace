/***************************************************************************
© F. Rousset 2005-

francois.rousset@umontpellier.fr
 
This file is part of GSpace. This software is a computer program
whose purpose is to perform population genetic simulations.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>

#include "charfunc.hpp"

using namespace std;

/*-----------ln(gamma(xx))----------------------------*/
double gammln(double xx)
{
	double x, tmp, y, ser;
	static const double cof[14] = {57.1562356658629235, -59.5979603554754912,
								   14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
								   .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
								   -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
								   .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5};

	y = x = xx;
	tmp = x + 5.24218750000000000; //Rational 671/128.
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (int j = 0; j < 14; j++)
		ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005 * ser / x);
}

long double Sichelcharfunc(long double xx, std::array<double, 4> DispPars)
{
	// 1-DIMENSIONAL distrib
	/*DispPars[2]<0 is used as indicator that the distribution is the reciprocal Gamma fn.
	DispPars[1] is no longer xi but (limit) omega xi (=kappa dans charfunc.nb)*/
	if (DispPars[2] < 0)
	{
		if (xx == 0)
		{
			return 1;
		}
		else
		{
			// reciprocal gamma for x<>0
			long double arg = sqrt(2 * DispPars[1] * (1 - cos(xx)));
			long double fn_eval = bessel_k(arg, fabs(DispPars[0]), 2); //returns exp(arg) K(arg; DispPars[0])
			if (isnan(fn_eval) == 1)
			{
				cerr << "fn_eval (1) in Sichelcharfunc is NaN";
				cin.get();
				exit(-1);
			}
			fn_eval *= exp(-arg);
			fn_eval /= pow(arg, DispPars[0]);
			if (isnan(fn_eval) == 1)
			{
				cerr << "fn_eval (2) in Sichelcharfunc is NaN";
                cin.get();
				exit(-1);
			}
			fn_eval /= exp(gammln(-DispPars[0]));
			if (isnan(fn_eval) == 1)
			{
				cerr << "fn_eval (3) in Sichelcharfunc is NaN";
				cin.get();
				exit(-1);
			}
			fn_eval *= pow((long double)(2.), (long double)(DispPars[0] + 1));
			return fn_eval;
		}
	}
	else
	{
		//Sichel...
		long double arg = sqrt(1. + 2 * DispPars[1] * (1 - cos(xx)) / DispPars[2]);
		/* 	bessel_k(DispPars[2]*arg, fabs(DispPars[0]), 1);
        works only  if DispPars[2]*arg and DispPars[2] are both <754...
        more generally one needs bessel_k(...,2)
        BUT note that this returns exp(x) K(x; nu)*/

		long double fn_eval = bessel_k(DispPars[2] * arg, fabs(DispPars[0]), 2);
		if (isnan(fn_eval) == 1)
		{
			cerr << "fn_eval (1) in Sichelcharfunc is NaN";
			cin.get();
			exit(-1);
		}
		fn_eval /= pow(arg, DispPars[0]);
		if (isnan(fn_eval) == 1)
		{
			cerr << "fn_eval (2) in Sichelcharfunc is NaN";
			cin.get();
			exit(-1);
		}
		fn_eval /= bessel_k(DispPars[2], fabs(DispPars[0]), 2);
		if (isnan(fn_eval) == 1)
		{
			cerr << "fn_eval (3) in Sichelcharfunc is NaN";
			cin.get();
			exit(-1);
		}
		fn_eval *= exp(DispPars[2] * (1. - arg));
		if (isnan(fn_eval) == 1)
		{
			cerr << "fn_eval (4) in Sichelcharfunc is NaN";
			cin.get();
			exit(-1);
		}
		return fn_eval;
	}
}
