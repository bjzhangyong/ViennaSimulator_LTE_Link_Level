/* file: turbo_rate_matcher_bit_selection_and_pruning_mapping.c
 LTE Turbo rate matching bit selection and pruning, as of TS 36.212, Section 5.1.4.2.
 (c) Josep Colom Ikuno, INTHFT
 jcolom@nt.tuwien.ac.at
 www.nt.tuwien.ac.at
 
 MEXed to LTE_common_turbo_rate_matcher_bit_selection_and_pruning_mapping

 input :   G        ... total number of bits available for the
                        transmission of one transport block
           N_ir     ... soft buffer size for the transport block
           Q_m      ... 2 when using QPSK, 4 for 16QAM and 6 for 64QAM
           N_l      ... 1 for blocks mapped onto one transmission layer
                        and 2 for blocks mapped onto two or four 
                        transmission layers.
           rv_idx   ... redundancy version number for this transmission
                        (rv_idx=0,1,2 or 3).
           K_pi     ... Size of the sub-block interleaver (multiple of
                        32). As in TS36.212, Section 5.1.4
           R_tc     ... Number of rows of the matrix in the sub-block
                        interleaver.
           C        ... number of code blocks computed in TS36.121,
                        sublclause 5.1.2.
           r        ... codeword we are processing (0, 1, 2,...) r-th

 output:   mapping  ... the mapping that tells where bits go
           k_0      ... the value of k_0. Can be useful for debugging purposes
*/

#include <mex.h>
#include <stdlib.h>
#include <math.h>

/* Input Arguments */
#define INPUT_G          prhs[0]
#define INPUT_N_IR       prhs[1]
#define INPUT_Q_M        prhs[2]
#define INPUT_N_l        prhs[3]
#define INPUT_RV_IDX     prhs[4]
#define INPUT_K_PI       prhs[5]
#define INPUT_R_TC       prhs[6]
#define INPUT_C          prhs[7]
#define INPUT_F          prhs[8]
#define INPUT_R          prhs[9]
#define INPUT_N_D        prhs[10]
#define INPUT_MODE       prhs[11]
#define INPUT_MODE_DATA  prhs[12]

/* Output Arguments */
#define OUTPUT_1         plhs[0]
#define OUTPUT_2         plhs[1]
#define OUTPUT_3         plhs[2]
#define OUTPUT_4         plhs[3]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{

    double   G;
	double   N_ir;
	double   Q_m;
	double   N_l;
	double   rv_idx;
	double   K_pi;
	double   R_tc;
	double   C;
    double   F;
	double   r;
	double   K_w;
	double   N_cb;
	double   G_prima;
	double   gamma;
	double   k_0;
    double   N_d;
    double   mode;
	
	mwSize   E;
	mwSize   k,j,l,index,null_pos_index;
    mwSize   length_null_array;
	
    double        *mapping;
    unsigned char *tx_bits;
    double        *null_positions;
    double        *null_positions_input;
	double        *output_k0;
    double        *real_number_of_nulls_used;
    bool          null_positions_is_empty;
    
    /* Check for proper number of arguments */
    if ((nrhs != 13 )||(nlhs  < 1)||(nlhs  > 4)) {
        mexErrMsgTxt("Usage: [mapping k_0 null_bits_pos] = LTE_common_turbo_rate_matcher_bit_selection_and_pruning_mapping(G,N_ir,Q_m,N_l,rv_idx,K_pi,R_tc,C,r,N_d,mode,mode_data)");
    } else {
		/* Get input parameters */
        G         = *mxGetPr(INPUT_G);          /* Target number of bits to generate */
		N_ir      = *mxGetPr(INPUT_N_IR);       /* Soft bit buffer size */
		Q_m       = *mxGetPr(INPUT_Q_M);        /* Modulation order */
		N_l       = *mxGetPr(INPUT_N_l);        /* Number of layers */
		rv_idx    = *mxGetPr(INPUT_RV_IDX);     /* Receiver redundancy version */
		K_pi      = *mxGetPr(INPUT_K_PI);       /* Subblock interleaver size */
		R_tc      = *mxGetPr(INPUT_R_TC);       /* Number of rows defined in TS 36.212 subclause 5.1.4.1.1 */
		C         = *mxGetPr(INPUT_C);          /* Total number of Codeblocks C */
        F         = *mxGetPr(INPUT_F);          /* Total number of filler bits F */
		r         = *mxGetPr(INPUT_R);          /* Codeblock index */
        N_d       = *mxGetPr(INPUT_N_D);        /* Number of <NULL> bits per subblock */
		mode      = *mxGetPr(INPUT_MODE);       /* 1 for TX, 2 for RX */
        
        if(mode==1) {
            tx_bits =  mxGetPr(INPUT_MODE_DATA);        /* The data to send */
        } else {
            null_positions_input    = mxGetPr(INPUT_MODE_DATA);  /* Position of the <NULL> bits in the data to send */
            null_positions_is_empty = mxIsEmpty(INPUT_MODE_DATA);
        }
        
        
		/* Total size of the circular buffer */
		K_w = 3*K_pi;
		
		/* Check if the buffer at receive side is big enough */
		N_cb = floor(N_ir/C) <= K_w ? floor(N_ir/C):K_w;
		/* G_prima -> number of symbols sent per layer */
		G_prima = G / (N_l*Q_m);
		/* gamma -> the rest after doing the aforementioned division */
		gamma = fmod(((float)G_prima),((float)C));
		
		/* Set the size of each rate-matched codeword */
		if(r <= (C-gamma-1)) {
			/* Number of bits to transmit */
			E = N_l*Q_m*floor(G_prima/C);
		}
		else {
			/* Number of bits to transmit */
			E = N_l*Q_m*ceil(G_prima/C);
		}
		
		/* Set k_0 */
		k_0 = R_tc * (2*ceil(N_cb/(8*R_tc))*rv_idx+2);

		k=0;
		j=0;
		
		/* Create ouput */
        OUTPUT_1  = mxCreateDoubleMatrix(1, E, mxREAL);
        mapping = mxGetPr(OUTPUT_1);
        OUTPUT_2 = mxCreateDoubleMatrix(1, 1, mxREAL);
        output_k0 = mxGetPr(OUTPUT_2);

        OUTPUT_4 = mxCreateDoubleMatrix(1, 1, mxREAL);
        real_number_of_nulls_used = mxGetPr(OUTPUT_4);
        
        if(mode==1) {
            OUTPUT_3 = mxCreateDoubleMatrix(1,(mwSize)((N_d*3+F*2) * ceil((double)(E/K_pi))), mxREAL); /* The maximum size of the vector containing the <NULL> positions */
            null_positions = mxGetPr(OUTPUT_3);
        }
        
        if(mode==2) {
            length_null_array = (mwSize)mxGetN(INPUT_MODE_DATA);
            /*mexPrintf("%d ",length_null_array);*/
            OUTPUT_3 = mxCreateDoubleMatrix(1,length_null_array, mxREAL);
            null_positions = mxGetPr(OUTPUT_3);
            for(l=0;l<length_null_array;l++) {
                null_positions[l] = null_positions_input[l];
            }
        }
        
		output_k0[0] = k_0;
		
        null_pos_index = 0;
		while(k<E) {
			index = ((mwSize)fmod(((float)k_0)+((float)j),((float)N_cb)));
			
            /* TX */
            if(mode==1) {
                if(tx_bits[index]!=3) {
                    mapping[k] = index;
                    k++;
                } else {
                    /* mexPrintf("%d ",index);mexPrintf("%d ",null_pos_index); */
                    null_positions[null_pos_index] = index;
                    null_pos_index++;
                }
                j++;
            }
            
            /* RX */
            else {
                /* null_pos_index = (mwSize)fmod((float)null_pos_index,(float)(length_null_array)); */
                /* mexPrintf("Index:%d ",index); */
                /* mexPrintf("Null index:%d/%d",null_pos_index,length_null_array); */
                if(null_positions_is_empty || (index!=null_positions[null_pos_index])) {
                    mapping[k] = index;
                    k++;
                    
                } else {
                    null_pos_index++;
                }
                /* mexPrintf("\n"); */
                j++;
            }
		}
    }
    real_number_of_nulls_used[0] = null_pos_index;
    return;
}
