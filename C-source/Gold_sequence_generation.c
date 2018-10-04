 /*Gold_sequence_generation - generates the gold sequence needed to scramble and descramble the coded bits in each of the codewords for PDSCH.
 DEFINED IN STANDARD 3GPP TS 36.211 V8.3.0 (2008-06) Section 7.2. 
  From V8.3.0, fast forward was added to the pseudo ramdon sequence generator. The value is N_c=1600. 
 gold_seq = Gold_sequence_generation(output_seq_length, NIDcell, NIDue, subframe, codeword)
 Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
 (c) 2009 by INTHFT
 www.nt.tuwien.ac.at

 MEXed to LTE_common_gen_gold_sequence

 input :   output_seq_length   ... [1 x 1]double - length of coded bits from one codeword = length of gold sequence
           NIDcell             ... [1 x 1]double cell identity
           NIDue               ... [1 x 1]double UE identity
           subframe            ... [1 x 1]double number of the subframe transmitted
           codeword            ... [1 x 1]double number of the codeword transmitted
 output:   gold_seq            ... [1 x # coded bits]double - the gold sequence for one codeword

 date of creation: 2009/03/19
*/

#include <stdlib.h>
#include <math.h>
#include <mex.h>

/* Input Arguments */
#define INPUT_length     prhs[0]
#define INPUT_NIDcell    prhs[1]
#define INPUT_NIDue      prhs[2]
#define INPUT_subframe   prhs[3]
#define INPUT_codeword   prhs[4]

/* Output Arguments */
#define OUTPUT_gold_seq  plhs[0]

/* fast forward offset value for pseudo-random sequence generator after v8.3.0 of 36.211	*/
/*#define	NC_OFFSET	0	*/	/* for older release, until v8.2.0	*/
#define	NC_OFFSET	1600		/* for new release, after v8.3.0	*/
#define	MAX_output_seq_length	800000	/* [NOTE] this is a tempolary value as maximum lenght of output sequence	*/

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{

  /*initial states of m-sequence generators */
  /* PN register state of the first m-seqence   */
  unsigned char mseq_state0[31 + MAX_output_seq_length + NC_OFFSET];
  /* PN register state of the second m-seqence  */
  unsigned char mseq_state1[31 + MAX_output_seq_length + NC_OFFSET];
  unsigned char mseq_state1_temp[31];
  
  long seed1_dec; /* initial state expressed as decimal value   */
  
  mwSize output_seq_length;
  double NIDcell;
  double NIDue;
  double subframe;
  double codeword;
  unsigned char *gold_seq;

  int	deg;		/* order of the polynomials poly0 and poly1 */

  mwSize i, ji;     /* loop variables   */
  int help_convert; /* temporary variable for the dec to binary conversion  */
  int index_help_convert;	/* invert parameter for the bibary conversion	*/
  
  /* Check for proper number of arguments */
  if ((nrhs != 5 )||(nlhs  < 1)||(nlhs  > 2)) 
  {
	mexErrMsgTxt("Usage: s = LTE_common_scrambling(b, NIDcell, NIDue, subframe, codeword, mode)");
  } 
  else 
  {
	/* Get input parameters */
	output_seq_length = (mwSize) (*mxGetPr(INPUT_length));
	NIDcell           = *mxGetPr(INPUT_NIDcell);
	NIDue             = *mxGetPr(INPUT_NIDue);
	subframe          = *mxGetPr(INPUT_subframe);
	codeword          = *mxGetPr(INPUT_codeword);
    
    /* create the output vector */
    OUTPUT_gold_seq = mxCreateLogicalMatrix(1, output_seq_length);
    gold_seq = mxGetPr(OUTPUT_gold_seq);

  /*initial states of m-sequence generators */
  /*DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03) Section 7.2. this initial state is also same with v8.8.0    */
	for(i = 0 ; i < 31 + output_seq_length + NC_OFFSET ; i++)
	{
		mseq_state0[i] = 0;
		mseq_state1[i] = 0;
	}
	mseq_state0[0] = 1;
	
    deg = 31;    
    seed1_dec = pow(2, 14) * NIDue + pow(2, 13) * codeword + pow(2, 9) * (subframe - 1) + NIDcell;
    
    /* Conversion of the seed into binary number for the second m-sequence that is calculated as decimal from input pars */
    for(help_convert = 30 ; help_convert >= 0; help_convert--) /* return to oroginal setting in original version */
    {
        mseq_state1[help_convert] = seed1_dec / (1 << help_convert);
        seed1_dec = seed1_dec - mseq_state1[help_convert] * (1 << help_convert);
    }
    
    for (i = 0; i < output_seq_length + NC_OFFSET ; i++) 
    {
		if(i >= NC_OFFSET)
			gold_seq[i - NC_OFFSET] = mseq_state0[i] ^ mseq_state1[i];
		mseq_state0[i + deg] = mseq_state0[i] ^ mseq_state0[i + 3];
		mseq_state1[i + deg] = (mseq_state1[i] ^ mseq_state1[i + 1]) ^ (mseq_state1[i + 2] ^ mseq_state1[i + 3]);
    }
  }
  return;
}
