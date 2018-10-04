/* file: byte2bit.c

   Performs a conversion between an uint8 vector representing a sequence 
   of bytes to a vector of logicals representing the bits in the sequence.
   
   MEXed to LTE_common_soft_bit_interleaver

   The calling syntax is:

      [output] = byte2bit(input)

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

/* Bit masks */
#define BIT_1   0x01
#define BIT_2   0x02
#define BIT_3   0x04
#define BIT_4   0x08
#define BIT_5   0x10
#define BIT_6   0x20
#define BIT_7   0x40
#define BIT_8   0x80

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    unsigned char   *input_p;
    unsigned char	*output_p;
    mwSize           DataLength;
    mwSize           output_length;
    mwSize           i;
    mwSize           output_index;
    
    /* Check if the input is of the correct type */
    if ( mxIsUint8(INPUT) != 1)
        mexErrMsgTxt("Input must uint8.");
    
    /* Check for proper number of arguments */
    if ((nrhs < 1 )||(nlhs  > 1)) {
        mexErrMsgTxt("Usage: [output] = [output] = byte2bit(input)");
    } else {
        input_p = mxGetPr(INPUT);
        DataLength = mxGetN(INPUT); /* number of data bytes */
        output_length = DataLength*8;
        
        /* create the output vector */
        OUTPUT = mxCreateLogicalMatrix(1, output_length);
        output_p = mxGetPr(OUTPUT);
        
        /* Populate the output */
        for(i=0;i<DataLength;i++) {
            output_index = 8*i;
            output_p[output_index]   = input_p[i]&BIT_1;
            output_p[output_index+1] = (input_p[i]&BIT_2)>>1;
            output_p[output_index+2] = (input_p[i]&BIT_3)>>2;
            output_p[output_index+3] = (input_p[i]&BIT_4)>>3;
            output_p[output_index+4] = (input_p[i]&BIT_5)>>4;
            output_p[output_index+5] = (input_p[i]&BIT_6)>>5;
            output_p[output_index+6] = (input_p[i]&BIT_7)>>6;
            output_p[output_index+7] = (input_p[i]&BIT_8)>>7;
            #ifdef DEBUG
            mexPrintf("Input number: %d\n",i);
            mexPrintf("input:  %d\n",input_p[i]);
            mexPrintf("output: %d%d%d%d%d%d%d%d\n",output_p[output_index+7],output_p[output_index+6],output_p[output_index+5],output_p[output_index+4],output_p[output_index+3],output_p[output_index+2],output_p[output_index+1],output_p[output_index]);
            #endif
        }
    }
    
    return;
}
