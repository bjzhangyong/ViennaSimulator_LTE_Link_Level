/* file: bit2byte.c

   Performs a conversion between a logical vector representing a sequence 
   of bits to an uint8 vector representing the bits in the sequence.
   MEXed to LTE_common_bit2byte

   The calling syntax is:

      [output] = byte2bit(input)

   Josep Colom Ikuno, INTHFT
   jcolom@nt.tuwien.ac.at
   www.nt.tuwien.ac.at

*/

#include <mex.h>
#include <math.h>
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
    unsigned char   *input_p;
    unsigned char	*output_p;
    mwSize    DataLength;
    mwSize    output_length;
    mwSize    i;
    mwSize    output_index;
    double double_output_data_length;
    double rounded_output_data_length;
    unsigned char current_byte;
    mwSize dims[2];
    
    /* Check if the input is of the correct type */
    if ( (mxIsLogical(INPUT) != 1) && (mxIsClass(INPUT,"uint8") != 1))
        mexErrMsgTxt("Input must Logical or uint8.");
    
    /* Check for proper number of arguments */
    if ((nrhs < 1 )||(nlhs  > 1)) {
        mexErrMsgTxt("Usage: [output] = [output] = bit2byte(input)");
    } else {
        input_p = mxGetPr(INPUT);
        DataLength = mxGetN(INPUT); /* number of data bytes */
        double_output_data_length = (double)DataLength/8;
        rounded_output_data_length = ceil(double_output_data_length);
        output_length = (size_t)rounded_output_data_length;
        #ifdef DEBUG
        /* mexPrintf("sizeof(size_t)=%d,sizeof(mwSize)=%d,sizeof(int)=%d\n",sizeof(size_t),sizeof(mwSize),sizeof(int)); */
        mexPrintf("Input length: %f\n",double_output_data_length);
        mexPrintf("Output length (double): %f\n",rounded_output_data_length);
        mexPrintf("Output length (ulong):  %d\n",output_length);
        #endif
        if(double_output_data_length!=rounded_output_data_length) {
            mexErrMsgTxt("Input length must be a multiple of eight");
        }
        
        
        /* create the output vector */
        dims[0] = 1;
        dims[1] = (mwSize)output_length;
        OUTPUT = mxCreateNumericArray(2,dims,mxUINT8_CLASS,mxREAL);
        output_p = mxGetPr(OUTPUT);
        
        /* Populate the output */
        for(i=0;i<output_length;i++) {
            current_byte = input_p[i*8]; /* LSB */
            current_byte = current_byte + (input_p[i*8+1]<<1);
            current_byte = current_byte + (input_p[i*8+2]<<2);
            current_byte = current_byte + (input_p[i*8+3]<<3);
            current_byte = current_byte + (input_p[i*8+4]<<4);
            current_byte = current_byte + (input_p[i*8+5]<<5);
            current_byte = current_byte + (input_p[i*8+6]<<6);
            current_byte = current_byte + (input_p[i*8+7]<<7);
            output_p[i] = current_byte;
            #ifdef DEBUG
            mexPrintf("Input number: %d\n",i);
            mexPrintf("input: %d%d%d%d%d%d%d%d\n",input_p[i*8+7],input_p[i*8+6],input_p[i*8+5],input_p[i*8+4],input_p[i*8+3],input_p[i*8+2],input_p[i*8+1],input_p[i*8]);
            mexPrintf("output:  %d\n",output_p[i]);
            #endif
        }
    }
    
    return;
}
