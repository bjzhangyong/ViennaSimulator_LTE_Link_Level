/* file: Hard_decision.c

   Performs the hard decision at the last part of the Turbo decoder. If LLR >=0,
   then it's 1, if not it's 0.
   MEXed to LTE_rx_hard_decision

   The calling syntax is:

      [output] = Hard_decision(input)

      output = the hard decision

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
    double        *input_p;
    unsigned char *output_p;
    mwSize        i, DataLength;
    
    /* Check if the input is of the correct type */
    if ( mxIsDouble(INPUT) != 1)
        mexErrMsgTxt("Input must be double.");
    
            /* first input is the data word */
    input_p = mxGetPr(INPUT);
    DataLength = mxGetN(INPUT); /* number of data bits */
        
    /* create the output vector */
    OUTPUT = mxCreateLogicalMatrix(1,DataLength);
    output_p = mxGetPr(OUTPUT);
    
    /* Populate the output */
    for(i=0;i<DataLength;i++) {
        if(input_p[i]>=0) {
            output_p[i] =  1;
        } else {
            output_p[i] = 0;
        }
    }
    
    return;
}
