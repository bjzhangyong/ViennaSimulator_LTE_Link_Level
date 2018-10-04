/************************************************************************/
/* SOFT_SD_IMPROVED_KBEST_FUNC_CPP										*/
/* by Martin Mayer, February 2011, mmayer@nt.tuwien.ac.at				*/
/* implements improved k-best sphere decoder with soft decisions        */
/*																		*/
/* SOURCE FILE with main functions										*/
/************************************************************************/

#include "LTE_rx_soft_sd_improved_kbest_func.h"

/*=============================Improved k-best ML Vector Calculation=======================================*\
\*=========================================================================================================*/
// sphere decoder algorithm to get (approximated) maximum likelihood leaf node
void calc_ML_node_kbi(	double *R_re, double *R_im,		// upper triangular [nT x nT] matrix
					double *y_re, double *y_im, 		// received symbol column vector
					double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
					int *num_symbols, int *num_nodes,
					int max_symbols, int max_layer,
					node_t **node, node_t *ML_node 	// tree of nodes
					)
{
	// local variables:
	node_t  *node_ptr, *parent_ptr;
	double	radius;					// search radius
	int		parent, *num_nodesC;	// parents stores parent indices for next layer, num_nodesC stores unpruned node number of layer
	int		i_layer, i_node, kk, num_unpruned, current_symbol; // help variables for loops etc.
	
	// additional variables for improved k-best
	bool next_layer = false;		// indicates if search advances to deeper layer
	int	start_layer = max_layer-1;	// has to be max_layer-1! (start search layer)
	int node_offset = 0;			// node offset in search layer
	int node_offset_start = 0;		// start node offset in search layer
	int *node_offset_ML = 0;			// stores ML node offset of search layer
			node_offset_ML = (int*)mxCalloc(max_layer, sizeof(int));
			
	num_nodesC = (int*)mxCalloc(max_layer, sizeof(int)); 	// num_nodesC changes during execution (num_nodes has to stay the same, therefore use copy of it)
	num_nodesC[max_layer-1] = num_nodes[max_layer-1]; 		// only necessary for top layer, rest determined during algorithm

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// execute one smallest PED depth first search (zero forcing solution) to get initial value for radius
	radius = get_ZF_distance(	&R_re[0], &R_im[0], &y_re[0], &y_im[0], &symbol_alphabet_re[0], &symbol_alphabet_im[0],
								&num_symbols[0], max_symbols, max_layer, node, false );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
	parent = 0;
								
	while( start_layer > 0 )	// begin repeated breadth first search
	{		
		for( i_layer = start_layer; i_layer >= 0; i_layer-- ) // iterate through layers, top to bottom
		{
			// === STEP 1 - PED CALCULATION OF CURRENT LAYER ===
			// calculate PEDs of unpruned nodes in current layer
			num_unpruned = 0;
			
			for( i_node = 0; i_node < num_nodesC[i_layer]; i_node++ ) // iterate through symbols/nodes
			{							
				current_symbol = i_node%num_symbols[i_layer];

				node_ptr   = &node[i_layer][i_node];			// points to current node
				
				if(i_layer < max_layer-1)
				{
					if( (current_symbol==0) && (i_node > 0) )	parent++;
					
					if( parent >= num_nodes[i_layer+1] )
						parent = 0;
					
					parent_ptr = &node[i_layer+1][parent];	
				}
				else
					parent_ptr = NULL;
					
				// calculate partial Euclidean distance of current node
				calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
							current_symbol, max_symbols, i_layer, max_layer, false );
		
				// count unpruned nodes of current layer
				if( node_ptr->PED <= radius ) 
					num_unpruned++;
					
				//printf("layer = %3d %3d\n", i_layer, start_layer);
			}// node iteration
			
			// === STEP 2 - SORT NODES WITH ASCENDING PEDs ===
			node_quick_sort(node[i_layer], 0, num_nodesC[i_layer]-1);	
			// now we've got an array with length num_unpruned of sorted nodes
			
			// === STEP 3 - CHECK IF RESTRECTION IS EXCEEDED IN NEXT LAYER / ASSIGN PARENTS ===
			if( i_layer > 0 ) // only necessary in layers bigger than bottom_layer + 1
			{
				if( i_layer == start_layer )
				{					
					node_offset_start = node_offset;	// we begin at last node offset
					
					// determine number of nodes in next layer
					num_nodesC[i_layer-1] = 0;
					while( node_offset < num_unpruned )
					{
						num_nodesC[i_layer-1] += num_symbols[i_layer-1];	// add children to every unpruned node
						node_offset++;										// another node (and its children) was added
						
						if( num_nodesC[i_layer-1] == num_nodes[i_layer-1] )	// check if max allowed nodes not exceeded in next layer
							break;	// stop if max							
					}
						
					if( node_offset_start >= num_unpruned )	// all nodes of layer checked
					{
						next_layer  = true;		// advance to next layer
						break;
					}
					else
					{
						parent		= node_offset_start;
					}
				}
				else	// not in start layer, just use the regular k-best nodes
				{
					parent = 0;
				
					num_nodesC[i_layer-1] = num_unpruned * num_symbols[i_layer-1];
					if( num_nodesC[i_layer-1] > num_nodes[i_layer-1] )
						num_nodesC[i_layer-1] = num_nodes[i_layer-1];
				}
			}			
		}// layer iteration
				
		// === STEP 5 - CHECK PED FOR NEW RADIUS / new ML node ===

		if( next_layer == true )
		{
			node_offset = 0;	// start at beginning of "node array" (to check all nodes)
			parent 		= node_offset_ML[start_layer];	// begin one layer beneath start layer with "ML parent"
			
			start_layer--;		// start in deeper layer
			
			num_nodesC[start_layer] = num_nodes[start_layer];	
			next_layer = false;
		}	
		else
		{
			if( node[0][0].PED <= radius )
			{
				node_offset_ML[start_layer] = node_offset_start;	// new ML node offset in current start layer
				
				ML_node->PED = node[0][0].PED;				// new ML node
				for( kk = 0; kk < max_layer; kk++ )
					ML_node->PSV[kk] = node[0][0].PSV[kk];	// copy PSV to ML node
				
				radius = ML_node->PED;						// new search radius
			}
		}		
	}	
	
	/*
	for( i_layer = max_layer-1; i_layer>0; i_layer--)
		printf(" %d |", node_offset_ML[i_layer]);
	printf("\n");
		*/
	// free allocated memory
	mxFree(node_offset_ML);
	mxFree(num_nodesC);
}

/*===============================Calculate Counter Hypotheses==============================================*\
\*=========================================================================================================*/
// USES NORMAL K-BEST ALGORITHM FOR COUTNER HYPOTHESES!
void calc_CH_distances_kbi2(	double *R_re, double *R_im,	// upper triangular [nT x nT] matrix
						double *y_re, double *y_im, // received symbol column vector
						double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
						int *num_symbols, int *num_nodes,
						int max_symbols, int max_layer,
						node_t *ML_node, int *num_bits, bool *bittable, int total_bits, double *LLR, node_t **node )
{	
	// local variables:
	node_t  *node_ptr, *parent_ptr;	// pointers to actual node and its parent
	double	radius;					// search radius
	int		parent, *num_nodesC;	// parents stores parent indices for next layer, num_nodesC stores unpruned node number of layer
	int		i_layer, i_node, i_bit, j, b, num_unpruned, current_symbol; // help variables for loops etc.

	// distances for LLR calculation
	double	lambda_ML = ML_node->PED;
	double  *lambda_CH;
			lambda_CH = (double*)mxCalloc( total_bits, sizeof(double) ); // counter hypotheses distances
		
	num_nodesC = (int*)mxCalloc(max_layer, sizeof(int)); 	// num_nodesC changes during execution (num_nodes has to stay the same, therefore use copy of it)
	num_nodesC[max_layer-1] = num_nodes[max_layer-1]; 		// only necessary for top layer, rest determined during algorithm
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	i_bit = 0; // index of counter hypothesis total count
	// iterate through all bits for counter hypotheses (j is layer, b is bit in layer)
	// CH distance and LLR is calculated for each bit
	for( j = 0; j < max_layer; j++ )
	{
		for( b = 0; b < num_bits[j]; b++ )
		{
			//xxxxxxxxxxxxxxxxxxx calculate counter hypothesis of actual bit xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// get CH[i_bit] ZF radius
			radius = get_ZF_distance_ch(&R_re[0], &R_im[0], &y_re[0], &y_im[0], &symbol_alphabet_re[0], &symbol_alphabet_im[0],
										&num_symbols[0], max_symbols, max_layer, node, j, ML_node, &num_bits[0], &bittable[0], total_bits, i_bit, false );						
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
			for( i_layer = max_layer-1; i_layer >= 0; i_layer-- ) // iterate through layers, top to bottom
			{ 
				// === STEP 1 - PED CALCULATION OF CURRENT LAYER ===
				// calculate PEDs of unpruned nodes in current layer
				num_unpruned = 0;
				parent = 0;
				
				for( i_node = 0; i_node < num_nodesC[i_layer]; i_node++ ) // iterate through symbols/nodes
				{				
					current_symbol = i_node%num_symbols[i_layer];
					
					node_ptr   = &node[i_layer][i_node];			// points to current node
					
					if(i_layer < max_layer-1)
					{
						if( (current_symbol==0) && i_node > 0 )	parent++;
						parent_ptr = &node[i_layer+1][parent];	
					}
					else
						parent_ptr = NULL;
							
					// check if current symbol bit equals ML bit
					if ( ( i_layer == j ) && ( bittable[ i_bit + total_bits*current_symbol ] == bittable[ i_bit + total_bits*ML_node->PSV[i_layer] ] ) )
					{ // node is pruned
						//node_ptr->PED = inf; // this will prune ML subset
						calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
									current_symbol, max_symbols, i_layer, max_layer, true );
					}
					else
					{	
						// calculate partial Euclidean distance of current node
						calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
									current_symbol, max_symbols, i_layer, max_layer, false );
						if(node_ptr->PED <= radius) // count unpruned nodes
							num_unpruned++;
					}
				}// node iteration
				// === STEP 2 - SORT NODES WITH ASCENDING PEDs ===
				node_quick_sort(node[i_layer], 0, num_nodesC[i_layer]-1);
				// now we've got an array with length num_unpruned of sorted nodes
				
				// === STEP 3 - CHECK IF RESTRECTION IS EXCEEDED IN NEXT LAYER / ASSIGN PARENTS ===
				if( i_layer > 0 ) // only necessary in layers bigger than bottom_layer + 1
				{	
					num_nodesC[i_layer-1] = num_unpruned * num_symbols[i_layer-1];
					if( num_nodesC[i_layer-1] > num_nodes[i_layer-1] )
						num_nodesC[i_layer-1] = num_nodes[i_layer-1];	
				}			
			}// layer iteration
				
			lambda_CH[i_bit] = node[0][0].PED;	// new CH node found
				
			// calculate LLR of current bit
			if( !bittable[ i_bit + total_bits*ML_node->PSV[j] ] )
				LLR[i_bit] = ML_node->PED - lambda_CH[i_bit];
			else
				LLR[i_bit] = lambda_CH[i_bit] - ML_node->PED;
			
			i_bit++;	// next LLR bit...
		
			//xxxxxxxxxxxxxxx end calculate counter hypothesis of actual bit xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		}
	}
	
	// free allocated memory
	mxFree(num_nodesC);
	mxFree(lambda_CH);
}

/*===============================Calculate Counter Hypotheses==============================================*\
\*=========================================================================================================*/
// USES IMPROVED K-BEST ALGORITHM FOR COUTNER HYPOTHESES
void calc_CH_distances_kbi(	double *R_re, double *R_im,	// upper triangular [nT x nT] matrix
						double *y_re, double *y_im, // received symbol column vector
						double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
						int *num_symbols, int *num_nodes,
						int max_symbols, int max_layer,
						node_t *ML_node, int *num_bits, bool *bittable, int total_bits, double *LLR, node_t **node )
{	
	// local variables:
	node_t  *node_ptr, *parent_ptr;	// pointers to actual node and its parent
	double	radius;					// search radius
	int		parent, *num_nodesC;	// parents stores parent indices for next layer, num_nodesC stores unpruned node number of layer
	int		i_layer, i_node, i_bit, j, b, num_unpruned, current_symbol; // help variables for loops etc.

	// distances for LLR calculation
	double	lambda_ML = ML_node->PED;
	double  *lambda_CH;
			lambda_CH = (double*)mxCalloc( total_bits, sizeof(double) ); // counter hypotheses distances
	
	// additional variables for improved k-best
	bool next_layer = false;		// indicates if search advances to deeper layer
	int	start_layer = max_layer-1;	// has to be max_layer-1! (start search layer)
	int node_offset = 0;			// node offset in search layer
	int *node_offset_CH;			// stores ML node offset of each search layer
		 node_offset_CH = (int*)mxCalloc(max_layer, sizeof(int));
	
	num_nodesC = (int*)mxCalloc(max_layer, sizeof(int)); 	// num_nodesC changes during execution (num_nodes has to stay the same, therefore use copy of it)
	num_nodesC[max_layer-1] = num_nodes[max_layer-1]; 		// only necessary for top layer, rest determined during algorithm
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	i_bit = 0; // index of counter hypothesis total count
	// iterate through all bits for counter hypotheses (j is layer, b is bit in layer)
	// CH distance and LLR is calculated for each bit
	for( j = 0; j < max_layer; j++ )
	{
		for( b = 0; b < num_bits[j]; b++ )
		{
			//xxxxxxxxxxxxxxxxxxx calculate counter hypothesis of actual bit xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// get CH ZF radius
			radius = get_ZF_distance_ch(&R_re[0], &R_im[0], &y_re[0], &y_im[0], &symbol_alphabet_re[0], &symbol_alphabet_im[0],
										&num_symbols[0], max_symbols, max_layer, node, j, ML_node, &num_bits[0], &bittable[0], total_bits, i_bit, false );						
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
			// initialize additional variables for improved k-best for every new CH
			next_layer  = false;		// indicates if search advances to deeper layer
			start_layer = max_layer-1;	// has to be max_layer-1!
			node_offset = 0;			// node offset in search layer
			
			while( start_layer > 0 )
			{
			
				for( i_layer = start_layer; i_layer >= 0; i_layer-- ) // iterate through layers, top to bottom
				{ 
					// === STEP 1 - PED CALCULATION OF CURRENT LAYER ===
					// calculate PEDs of unpruned nodes in current layer
					num_unpruned = 0;
					if( i_layer == start_layer-1 )
						parent = node_offset;	// parent index starts at node offset
					else
						parent = 0;
					
					for( i_node = 0; i_node < num_nodesC[i_layer]; i_node++ ) // iterate through symbols/nodes
					{					
						current_symbol = i_node % num_symbols[i_layer];
						
						node_ptr   = &node[i_layer  ][i_node];			// points to current node
						//parent_ptr = &node[i_layer+1][parents[i_node]]; // points to parent symbol (don't use parent_ptr in max_layer-1!, will be access violation)
						
						if(i_layer < max_layer-1)
						{
							if( (current_symbol==0) && i_node > 0 )	parent++;
							parent_ptr = &node[i_layer+1][parent];	
						}
						else
							parent_ptr = NULL;
						
						// check if current symbol bit equals ML bit
						if ( ( i_layer == j ) && ( bittable[ i_bit + total_bits*current_symbol ] == bittable[ i_bit + total_bits*ML_node->PSV[i_layer] ] ) )
						{ // node is pruned
						//	node_ptr->PED = inf; // this will prune ML subset
							calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
										current_symbol, max_symbols, i_layer, max_layer, true );
						}
						else
						{	
							// calculate partial Euclidean distance of current node
							calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
										current_symbol, max_symbols, i_layer, max_layer, false );

							if(node_ptr->PED <= radius) // count unpruned nodes
								num_unpruned++;
						}
					}// node iteration

					// === STEP 2 - SORT NODES WITH ASCENDING PEDs ===
					node_quick_sort(node[i_layer], 0, num_nodesC[i_layer]-1);
					// now we've got an array with length num_unpruned of sorted nodes
					
					// === STEP 3 - CHECK IF RESTRECTION IS EXCEEDED IN NEXT LAYER / ASSIGN PARENTS ===
					if( i_layer > 0 ) // only necessary in layers bigger than bottom_layer + 1
					{
						if( i_layer == start_layer )
						{					
							// determine remaining nodes
							/* example:
							- number of symbols in next layer: 16 (each node has 16 children)
							- number of unpruned nodes in current layer: num_unpruned = 39
							- k-best restriction: max_nodes = 256
							- number of nodes in next layer if all children are calculated: 39*16 = 624
							- because of our max_nodes restriction, we divide the 624 child nodes into ceil(624/256) = 3 arrays
							- we start with the first array (children 0 to 255) and maybe find a new (smaller) radius in the end
							- in the next run (start_layer is reached again), this radius gives us a new num_unpruned = 28
							- so there is a total of 28*16 = 448 children left, we checked the first 0 to 255 children already
							- in the actual run, we check children 255 to 448 and are finished in the current layer
							*/				
							
							// iterate over "children" and assign parents for next layer calculation
							num_nodesC[i_layer-1] = 0;
							for( i_node = node_offset; i_node < num_unpruned; i_node++ ) 
							{
								num_nodesC[i_layer-1] += num_symbols[i_layer-1];	// add children to every unpruned node
								if( num_nodesC[i_layer-1] > num_nodes[i_layer-1] )				// check if max_nodes not exceeded in next layer
								{
									num_nodesC[i_layer-1] -= num_symbols[i_layer-1];
									break;	// stop if max
								}						
							}// symbol iteration

							node_offset = i_node;	// save node offset for next search in this layer
					
							if( node_offset >= num_unpruned-1 )	// all nodes of layer checked
							{
								next_layer  = true;							// advance to next layer
								node_offset = node_offset_CH[start_layer];	// use node-array which contained the smallest distance
							}	
						}
						else	// not in start layer, just use the regular k-best nodes
						{
							num_nodesC[i_layer-1] = num_unpruned * num_symbols[i_layer-1];
							if( num_nodesC[i_layer-1] > num_nodes[i_layer-1] )
								num_nodesC[i_layer-1] = num_nodes[i_layer-1];
						}
					}			
				}// layer iteration
						
				// === STEP 5 - CHECK PED FOR NEW RADIUS / new ML node ===
				if( node[0][0].PED <= radius )
				{
					node_offset_CH[start_layer] = node_offset;	// new ML node offset in current start layer
					
					lambda_CH[i_bit] = node[0][0].PED;	// new CH node found
					radius  		 = node[0][0].PED;
				}
							
				if( next_layer == true )
				{
					start_layer--;		// start in deeper layer	
					node_offset = 0;
					next_layer = false;
				}
			}
							
			// calculate LLR of current bit
			if( !bittable[ i_bit + total_bits*ML_node->PSV[j] ] )
				LLR[i_bit] = ML_node->PED - lambda_CH[i_bit];
			else
				LLR[i_bit] = lambda_CH[i_bit] - ML_node->PED;
			
			i_bit++;	// next LLR bit...
		
			//xxxxxxxxxxxxxxx end calculate counter hypothesis of actual bit xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		}
	}
	
	// free allocated memory
	mxFree(num_nodesC);
	mxFree(node_offset_CH);
	mxFree(lambda_CH);
}

/*=============================Soft Sphere Decoder=========================================================*\
\*=========================================================================================================*/
// Soft Sphere Decoder
void ssd_kbest_improv_c(
			double *R_re, double *R_im,	// upper triangular [nT x nT] matrix
			double *y_re, double *y_im, 	// received symbol column vector
			double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
			double *LLR,				// maximum likelihood symbol vector (output of sesd_c)
			int   *M,				// number of bits in corresponding layer [nT x 1]
			bool  *bittable,		// bittable
			int total_bits,			// total number of bits
			int nSym,				// number of transmitted symbol vectors
			int nT,					// number of transmit antennas in MIMO system
			int nR,					// number of receive antennas in MIMO system
			int k					// number of "k-best" nodes
			)
{
	// variable declarations:
	node_t	**node, *ML_node; // tree of nodes, ML-node
	
	int	max_layer = nT;		// "height" of tree
	int	max_symbols;		// max number of symbols of every antenna (width of symbol alphabet matrix)
	int	*num_symbols;		// number of symbols per layer
	int	*num_nodes;			// number of nodes   per layer
	int	i_layer, i_node, kk, ss;// help variables for loops etc.
		
	// max likelihood node, stores ML node for every received symbol vector
 	ML_node 		= (node_t*)mxMalloc(sizeof(node_t));
	ML_node->PSV 	= (int*)mxCalloc(max_layer, sizeof(int));		

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// calculate tree information

	// array that indicates number of symbols in each layer (each layer=antenna may have different symbol alphabets)
	num_nodes   = (int*)mxCalloc( max_layer, sizeof(int)); // nT = max_layer layers
	num_symbols = (int*)mxCalloc( max_layer, sizeof(int)); // nT = max_layer layers
	
	max_symbols = 0;

	for( i_layer = max_layer-1; i_layer >= 0; i_layer-- ) // iterate through layers
	{
		// determine number of symbols in current layer
		num_symbols[i_layer] = 2;		// at least 2 symbols (1 bit = 0,1... 2 symbols)

		for( kk = 1; kk < M[i_layer]; kk++ )
			num_symbols[i_layer] *= 2;	// every bit doubles number of symbols in current layer

		// determine max number of symbols of all layers
		if( max_symbols < num_symbols[i_layer] )
			max_symbols = num_symbols[i_layer];

		// determine number of necessary nodes in each layer to ensure k-best
		if( i_layer == max_layer-1 )
			num_nodes[i_layer] = num_symbols[i_layer];
		else
		{
			if( k > num_nodes[i_layer+1] )
				num_nodes[i_layer] = num_symbols[i_layer] * num_nodes[i_layer+1];	// k is bigger than nodes in layer
			else
				num_nodes[i_layer] = num_symbols[i_layer] * k;						// k is smaller than node number, use only k best nodes
		}
		// print structure
		//printf("    layer %ld | %ld symbols | %ld nodes\n", i_layer, num_symbols[i_layer], num_nodes[i_layer]);
	}// layer iteration
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// allocate memory optimized "tree" structure
	node = (node_t**)mxMalloc( sizeof(node_t) * max_layer );
	
	for( i_layer =0; i_layer < max_layer; i_layer++ )
	{
		node[i_layer] = (node_t*)mxMalloc( sizeof(node_t)*num_nodes[i_layer] );
		for( i_node = 0; i_node < num_nodes[i_layer]; i_node++ )
		{
			node[i_layer][i_node].PSV = (int*)mxCalloc( max_layer, sizeof(int) ); // allocate max PSV length for every PSV, easier to handle
			node[i_layer][i_node].PED = 0;
		}
	}
	
	// iterate over received symbol vectors in current channel realization
	for(ss=0; ss<nSym; ss++)
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// get max likelihood symbol vector
		calc_ML_node_kbi(	&R_re[0], &R_im[0], &y_re[ss*max_layer], &y_im[ss*max_layer], &symbol_alphabet_re[0], &symbol_alphabet_im[0], 
						&num_symbols[0], &num_nodes[0],	max_symbols, max_layer, node, ML_node );	
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
		// calculate counter hypotheses and resulting LLRs with k best tree search
		calc_CH_distances_kbi2(	&R_re[0], &R_im[0], &y_re[ss*max_layer], &y_im[ss*max_layer], &symbol_alphabet_re[0], &symbol_alphabet_im[0], 
							&num_symbols[0], &num_nodes[0], max_symbols, max_layer, 
							ML_node, M, &bittable[0], total_bits, &LLR[ss*total_bits], node );
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}
	
	// free allocated memory
	for( i_layer =0; i_layer < max_layer; i_layer++ )
	{
		for( i_node = 0; i_node < num_nodes[i_layer]; i_node++ )
			mxFree(node[i_layer][i_node].PSV);	// free PSV
		mxFree(node[i_layer]);
	}
	mxFree(ML_node->PSV);
	mxFree(ML_node);
	mxFree(node);
	mxFree(num_nodes);
	mxFree(num_symbols);
}