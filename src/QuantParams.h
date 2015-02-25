/*
 * QuantParams.h
 *
 *  Created on: Feb 20, 2015
 *      Author: Jeremy Schwartzentruber
 */

#ifndef QUANT_PARAMS_H_
#define QUANT_PARAMS_H_

using namespace std;


class QuantParams {
public:
	QuantParams();
	QuantParams(double l , double _b0, double _b1) { lambda = l; b0 = _b0; b1 = _b1; }
	
	double lambda;
	double b0;
	double b1;
};


#endif // QUANT_PARAMS_H_
