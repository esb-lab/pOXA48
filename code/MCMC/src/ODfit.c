// ODfit model used by package pNUK73
//
// Created by RPM on 28/5/14.
// Copyright 2014. All rights reserved.
//
// Implements deterministic batch model with one bacterial type and a single-limiting resource to be used with R 
// To compile from R: system("R CMD SHLIB ODfit.c")
// To load from R (in Unix): dyn.load("ODfit.so")
// To call from R use deSolve package and see details in help on using compiled code.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h> 

static double parms[3];
#define rho parms[0]
#define m parms[1]
#define K parms[2]

/* initializer */
void initmod(void (* odeparms)(int * , double *))
{
    int N=3;
    odeparms(&N, parms);
}

/* derivatives and one output variable 
y[0] is R (Resource)
y[1] is B (Bacterial density)
*/
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if(ip[0]<1) error("nout should be at least 1");
    ydot[0] = -m*y[0]/(K+y[0])*y[1];
    ydot[1] = rho*m*y[0]/(K+y[0])*y[1];
    yout[0]=y[1];
}

/* Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    pd[0]=-(K*m*y[1])/((K+y[0])*(K+y[0]));
    pd[1]=-m*y[0]/(K+y[0]);
    pd[(*nrowpd)]=(rho*K*m*y[1])/(K+y[0])*(K+y[0]);
    pd[(*nrowpd)+1]=rho*m*y[0]/(K+y[0]);
    
}


