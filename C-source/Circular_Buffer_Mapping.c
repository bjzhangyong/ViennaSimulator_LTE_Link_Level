/* file: Circular_Buffer_Mapping.c

   Given an interleaver length K_pi, this function calculates the mapping
   that works as the circular buffer used in the LTE Turbo Code Rate Matcher.
   MEXed to LTE_common_turbo_rate_matching_circular_buffer_mapping

   The calling syntax is:

      [output] = Circular_Buffer_Mapping(K_pi)

      output = interleaver mapping

   Josep Colom Ikuno, INTHFT
   jcolom@nt.tuwien.ac.at
   www.nt.tuwien.ac.at

*/

#include <mex.h>
#include <stdlib.h>

/* Input Arguments */
#define INPUT       prhs[0]

/* Output Arguments */
#define OUTPUT      plhs[0]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    double *output_p;
    double *input_p;
    mwSize output_length;
    mwSize k_pi;
    mwSize k;
    
    /* Check if the input is of the correct type */
    if ( mxIsDouble(INPUT) != 1)
        mexErrMsgTxt("Input must be double.");
    
    /* Check for proper number of arguments */
    if ((nrhs < 1 )||(nlhs  > 1)) {
        mexErrMsgTxt("Usage: [output] = Circular_Buffer_Mapping(K_pi)");
    } else {
        /* first input is the data word */
        input_p = mxGetPr(INPUT);
        k_pi = (mwSize)(*input_p);
        output_length = k_pi*3;
        
        /* create the output vector */
        OUTPUT = mxCreateDoubleMatrix(1, output_length, mxREAL);
        output_p = mxGetPr(OUTPUT);
        
        /* Populate the output */
        for(k=0;k<k_pi;k++) {
            output_p[k] = (double)k;
            output_p[k_pi+2*k]   = k_pi+k;
            output_p[k_pi+2*k+1] = 2*k_pi+k;
        }
    }
    return;
}
