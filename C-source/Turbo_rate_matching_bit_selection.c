/* file: Turbo_rate_matching_bit_selection.c

   Performs the puncturing for the turbo code rate matcher, as specified in
   TS36.212, subclause 5.1.4.1.2. Works for bits (logicals)
   The interleave/deinterleave flag signals whether you want to interleave 
   or deinterleave the input. The possible values are:
     1: interleave
     0: deinterleave

 MEXed to LTE_common_turbo_rate_matching_bit_selection_and_pruning

   The calling syntax is:

      [output] = Bit_Interleaver(input, mapping_vector, interleave/deinterleave,output_length)
      Note that when de-interleaving, the size of the output codeword must be inputted.

      output = interleaved (punctured) input

   Josep Colom Ikuno, INTHFT
   jcolom@nt.tuwien.ac.at
   www.nt.tuwien.ac.at

*/

#include <mex.h>
#include <stdlib.h>

/* Input Arguments */
#define INPUT          prhs[0]
#define MAPPING        prhs[1]
#define INTERLEAVE     prhs[2]
#define OUTPUT_LENGTH  prhs[3]

/* Output Arguments */
#define OUTPUT      plhs[0]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    unsigned char   *input_p_char;
    unsigned char	*output_p_char;
    double   *input_p_double;
    double   *output_p_double;
    double   *mapping;
    mwSize    DataLength;
    mwSize i;
    mwSize index;
    double *interleave;
    double *output_length;
    
    /* Check for proper number of arguments */
    if ((nrhs < 3 )||(nlhs  > 1)) {
        mexErrMsgTxt("Usage: [output] = Turbo_rate_matching_bit_selection(input, mapping_vector, interleave/deinterleave,(output_length->only when de-rate-matching))");
    } else {
        interleave = mxGetPr(INTERLEAVE);
        DataLength = mxGetN(MAPPING); /* number of  output data bits */
        mapping = mxGetPr(MAPPING);
        
        /* Code that interleaves */
        if(*interleave==1) {
            if ( mxIsClass(INPUT,"uint8") != 1)
                mexErrMsgTxt("Input must be uint8 when interleaving.");
            if ((nrhs != 3 )||(nlhs  > 1))
                mexErrMsgTxt("Usage: [output] = Turbo_rate_matching_bit_selection(input, mapping_vector, interleave/deinterleave)");
            input_p_char = mxGetPr(INPUT);
            OUTPUT = mxCreateLogicalMatrix(1, DataLength);
            output_p_char = mxGetPr(OUTPUT);
            /* Initialise the output */
            for(i=0;i<DataLength;i++) {
                output_p_char[i] = 0;
            }
            /* Populate the output */
            for(i=0;i<DataLength;i++) {
                index = (int)mapping[i];
                output_p_char[i] = input_p_char[index];
            }
        }
        /* de-interleaving */
        else {
            if ( mxIsDouble(INPUT) != 1)
                mexErrMsgTxt("Input must be double when deinterleaving.");
            if ((nrhs != 4 )||(nlhs  > 1))
                mexErrMsgTxt("Usage: [output] = Turbo_rate_matching_bit_selection(input, mapping_vector, interleave/deinterleave,output_length)");
            output_length = mxGetPr(OUTPUT_LENGTH);
            input_p_double = mxGetPr(INPUT);
            /* Create 2-D, double-precision, floating-point mxArray initialized to 0 */
            OUTPUT = mxCreateDoubleMatrix(1,*output_length,mxREAL);
            output_p_double = mxGetPr(OUTPUT);
            /* Initialise the output, not necessary if initialising to 0 */
            /*for(i=0;i<*output_length;i++) {
                output_p_double[i] = -1; // used to signal an empty position
            } */
            /* Populate the output and combine repeated values */
            for(i=0;i<DataLength;i++) {
                index = (int)mapping[i];
                if(index>=(*output_length))
                    mexErrMsgTxt("Mapping points to a position that does not fit in the given output_length");
                /* if(output_p_double[index]==-1) output_p_double[index] = input_p_double[i]; */
                /* else output_p_double[index] = output_p_double[index] + input_p_double[i]; */
                /* Combine received values. Values are initialised to 0. */
                output_p_double[index] = output_p_double[index] + input_p_double[i];
            }
        }
    }
    return;
}
