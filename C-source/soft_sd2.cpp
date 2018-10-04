/* This MEX-function implements a soft output sphere decoder. 
   Based on the paper: C. Studer, M. Wenk, A. Burg, and H. Bölcskei: "Soft-Output Sphere 
   Decoding: Performance and Implementation Aspects", Asilomar 2006
  
   created: 12. April 2007, Christian Mehlführer, chmehl@nt.tuwien.ac.at
  
   compile in Matlab using: "mex -O soft_sd2.cpp" */

#include "mex.h"
#include <string.h>

#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))

#include "soft_sd_c2.cpp"

static char ident[] = "@(#)soft_sd.: Compiled Matlab function,"
                      " <chmehl@nt.tuwien.ac.at> Apr/12/2007";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/*  int mm,nn,nB,nA, debug_level=0, we_must_mxFree_a_im=0;
    double *y_re, *y_im, *a_re, *a_im, *s_tilde_re, *s_tilde_im, *comp;
    static char c[7];
    mxArray *s_hat;
    double *s_hat_re, *s_hat_im; */
  
  /* input variables */
  double *R_re, *R_im, *s_re, *s_im, *dist_ZF, *symbol_alphabet_re, *symbol_alphabet_im;
  bool *bittable;
  int kk, nR, nT, nSym, *M, *symbols_ZF_i, *Nr_of_symbols, total_bits;
  /* output variables */
  double *LLR;
  
  if (nrhs!=7)
    mexErrMsgTxt("7 input arguments required \n\n"
				 "[LLR] = soft_sd(R,s,dist_ZF,symbols_ZF,symbol_alphabet,bittable) ... soft Sphere Decoder\n\n"
				 "  R ... upper triangular matrix obtained from the QR decomposition of the channel H (complex)\n"
				 "  s ... received symbol vector, s=Q^H*y (nR x nSym) (complex)\n"
                 "  dist_ZF ... Distance of the zero forcing solution (real)\n"
                 "  symbols_ZF_i ... indices to symbols of the ZF solution (nT x nSym) (real integer)\n"
                 "  M ... number of bits in the corresponding layer (1 x nR) (real)\n"
				 "  symbol_alphabet ... for the demapping (2^M_max x nT) (complex)\n"
				 "  bittable ... matrix containing the bits according to the symbol_alphabet (M x 2^M) (logical)\n"
				 "  LLR  ... max-log-MAP approximation of the LLR values (M*nR) (real)\n\n");
  if (nlhs>1) 
    mexErrMsgTxt("One output lefthand argument required \n");

  /* check input variables */
  if ( ! mxIsComplex(prhs[0]) )
    mexErrMsgTxt("1st argument 'R' must be a complex-valued (nR x nT) matrix");
  if ( ! mxIsComplex(prhs[1]) )
    mexErrMsgTxt("2nd argument 's' must be a complex-valued (nR x nSym) matrix");
  if (   mxIsComplex(prhs[2]) )
    mexErrMsgTxt("3rd argument 'dist_ZF' must be a real-valued (1 x nSym) matrix");
  if (   mxIsComplex(prhs[3]) )
    mexErrMsgTxt("4th argument 'symbols_ZF_i' must be a real-valued (nT x nSym) integer matrix");
  if (   mxIsComplex(prhs[4]) )
    mexErrMsgTxt("5th argument 'M' must be a real-valued (1 x nT) integer matrix");
  if ( ! mxIsComplex(prhs[5]) )
    mexErrMsgTxt("6th argument 'symbol_alphabet' must be a complex-valued (2^M_max x nT) matrix");
  if ( ! mxIsLogical(prhs[6]) )
    mexErrMsgTxt("7th argument 'bittable' must be a logical (M x 2^M) matrix");
  
  
    nR = mxGetM(prhs[0]);    /*  number of receive antennas */
    nT = mxGetN(prhs[0]);    /*  number of transmit antennas */
    nSym = mxGetN(prhs[1]);  /*  Block size (number of transmitted symbol vectors) */
    
    // fetch input variables
    R_re = (double *)(mxGetPr(prhs[0]));                    /* fetch pointer to real part of R */
	R_im = (double *)(mxGetPi(prhs[0]));                    /* fetch pointer to imag part of R */
    s_re = (double *)(mxGetPr(prhs[1]));                    /* fetch pointer to real part of s */
	s_im = (double *)(mxGetPi(prhs[1]));                    /* fetch pointer to imag part of s */
    dist_ZF = (double *)(mxGetPr(prhs[2]));                 /* fetch ZF distance */
	symbols_ZF_i = (int *)(mxGetPr(prhs[3]));               /* fetch pointer to imag part of ZF solution indices */
	M = (int *)(mxGetPr(prhs[4]));                          /* fetch pointer to number of bits vector */
    symbol_alphabet_re = (double *)(mxGetPr(prhs[5]));      /* fetch pointer to real part of symbol alphabet */
	symbol_alphabet_im = (double *)(mxGetPi(prhs[5]));      /* fetch pointer to imag part of symbol alphabet */
    bittable = (bool *)(mxGetPr(prhs[6]));                  /* fetch pointer to real part of bit mapping table */

    /* allocate memory for output variables */
    total_bits = 0;
    for(kk=0; kk<nT; kk++)
        total_bits += M[kk]; 
    
    plhs[0] = mxCreateDoubleMatrix(total_bits, nSym, mxREAL );    /* output variable with Sphere Decoder solution */

    if(plhs[0] == NULL)
      mexErrMsgTxt("mxCreateDoubleMatrix failed(1)\n");

    LLR = mxGetPr(plhs[0]);                                 /* fetch pointer for output variable */

    if(LLR == NULL)
      mexErrMsgTxt("mxCreateDoubleMatrix failed(2)\n");


/*    LLR[0] = 27.3; */
/*    LLR[1] = -22.3; */
    
/*    mexPrintf("symbol_alphabet_re[0]=%f, symbol_alphabet_re[1]=%f, symbol_alphabet_re[2]=%f, symbol_alphabet_re[3]=%f, symbol_alphabet_re[4]=%f, symbol_alphabet_re[5]=%f, symbol_alphabet_re[6]=%f, symbol_alphabet_re[7]=%f,\n", */ 
/*			symbol_alphabet_re[0], symbol_alphabet_re[1],symbol_alphabet_re[2],symbol_alphabet_re[3],symbol_alphabet_re[4],symbol_alphabet_re[5],symbol_alphabet_re[6],symbol_alphabet_re[7]); */

    
    soft_sd_c2(&R_re[0], &R_im[0], 
               &s_re[0], &s_im[0], 
               dist_ZF, 
               &symbols_ZF_i[0], 
               &symbol_alphabet_re[0], &symbol_alphabet_im[0], 
               &bittable[0], 
               &LLR[0], 
               nSym, &M[0], nT, nR, total_bits);

  
} /* mexFunction */
