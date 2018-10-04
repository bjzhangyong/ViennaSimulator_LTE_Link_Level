// This MEX-function implements a soft output sphere decoder with ML approximation
/************************************************************************/
/* SOFT_SD_IMPROVED_KBEST_CPP											*/
/* by Martin Mayer, February 2011, mmayer@nt.tuwien.ac.at				*/
/* implements improved k-best sphere decoder with soft decisions        */
/*																		*/
/* MEX FUNCTION FILE that parses input data	from Matlab					*/
/************************************************************************/

#include "mex.h"
#include "LTE_rx_soft_sd_improved_kbest_func.cpp"

static char ident[] = "@(#)soft_sd.: Compiled Matlab function,"
                      " <mmayer@nt.tuwien.ac.at> Dec/2010";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // input variables
  double 	*R_re, *R_im, *s_re, *s_im, *symbol_alphabet_re, *symbol_alphabet_im;
  bool 		*bittable;
  int 		kk, nR, nT, nSym, *M, *Nr_of_symbols, total_bits, *k;
  // output variables
  double 	*LLR;
  
  if (nrhs!=6)
    mexErrMsgTxt("6 input arguments required \n\n"
				 "[LLR] = LTE_softsphere_kbest(R,s,dist_ZF,symbol_alphabet,bittable,k) ... soft Sphere Decoder\n\n"
				 "  R ... upper triangular matrix obtained from the QR decomposition of the channel H (complex)\n"
				 "  s ... received symbol vector, s=Q^H*y (nR x nSym) (complex)\n"
                 "  M ... number of bits in the corresponding layer (1 x nR) (real)\n"
				 "  symbol_alphabet ... for the demapping (2^M_max x nT) (complex)\n"
				 "  bittable ... matrix containing the bits according to the symbol_alphabet (M x 2^M) (logical)\n"
				 "  k ... 'k-best nodes'\n"
				 "  LLR  ... max-log-MAP approximation of the LLR values (M*nR) (real)\n\n");
  if (nlhs>1) 
    mexErrMsgTxt("One output lefthand argument required \n");

  // check input variables
  if ( ! mxIsComplex(prhs[0]) )
    mexErrMsgTxt("1st argument 'R' must be a complex-valued (nR x nT) matrix");
  if ( ! mxIsComplex(prhs[1]) )
    mexErrMsgTxt("2nd argument 's' must be a complex-valued (nR x nSym) matrix");
  if (   mxIsComplex(prhs[2]) )
    mexErrMsgTxt("3rd argument 'M' must be a real-valued (1 x nT) integer matrix");
  if ( ! mxIsComplex(prhs[3]) )
    mexErrMsgTxt("4th argument 'symbol_alphabet' must be a complex-valued (2^M_max x nT) matrix");
  if ( ! mxIsLogical(prhs[4]) )
    mexErrMsgTxt("5th argument 'bittable' must be a logical (M x 2^M) matrix");
  if (   mxIsComplex(prhs[5]) )
    mexErrMsgTxt("6th argument 'k' must be a real integer value"); 

 
    nR = mxGetM(prhs[0]);    //  number of receive antennas
    nT = mxGetN(prhs[0]);    //  number of transmit antennas
    nSym = mxGetN(prhs[1]);  //  Block size (number of transmitted symbol vectors)
    
    // fetch input variables
    R_re = (double *)(mxGetPr(prhs[0]));                    // fetch pointer to real part of R
	R_im = (double *)(mxGetPi(prhs[0]));                    // fetch pointer to imag part of R
    s_re = (double *)(mxGetPr(prhs[1]));                    // fetch pointer to real part of s
	s_im = (double *)(mxGetPi(prhs[1]));                    // fetch pointer to imag part of s
	M = (int *)(mxGetPr(prhs[2]));                         	// fetch pointer to number of bits vector
    symbol_alphabet_re = (double *)(mxGetPr(prhs[3]));      // fetch pointer to real part of symbol alphabet
	symbol_alphabet_im = (double *)(mxGetPi(prhs[3]));      // fetch pointer to imag part of symbol alphabet
    bittable = (bool *)(mxGetPr(prhs[4]));                 	// fetch pointer to real part of bit mapping table
	k = (int *)(mxGetPr(prhs[5]));							// fetch pointer to number of survivor nodes k
		
	// check k value
	if( *k < 1 )
		mexErrMsgTxt("k value for improved k-best sphere decoder too small! k has to be >= 1!");

    // allocate memory for output variables
    total_bits = 0;
    for(kk=0; kk<nT; kk++)
        total_bits += M[kk]; 
    
    plhs[0] = mxCreateNumericMatrix(total_bits,nSym,mxDOUBLE_CLASS,mxREAL) ;   // output variable with Sphere Decoder solution
    if(plhs[0] == NULL)	mexErrMsgTxt("mxCreateDoubleMatrix failed(1)\n");

    LLR = (double*)mxGetPr(plhs[0]);    // fetch pointer for output variable
    if(LLR == NULL)		mexErrMsgTxt("mxCreateDoubleMatrix failed(2)\n");
	
	// execute improved k-best soft sphere decoder
    ssd_kbest_improv_c(
				&R_re[0], &R_im[0], 
				&s_re[0], &s_im[0], 
				&symbol_alphabet_re[0], &symbol_alphabet_im[0], 
				&LLR[0], 	// OUTPUT ARRAY, LLR FOR EACH BIT
				&M[0],
				&bittable[0],
				total_bits,
				nSym,
				nT, nR, *k);
				
} /* mexFunction */
