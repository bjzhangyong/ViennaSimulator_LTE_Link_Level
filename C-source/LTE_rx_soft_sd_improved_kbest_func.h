/************************************************************************/
/* SOFT_SD_IMPROVED_KBEST_FUNC_H										*/
/* by Martin Mayer, February 2011, mmayer@nt.tuwien.ac.at				*/
/* implements improved k-best sphere decoder with soft decisions        */
/*																		*/
/* HEADER FILE with additional helper functions							*/
/************************************************************************/

#ifndef SSD_FUNCTIONS_H_
#define SSD_FUNCTIONS_H_

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

static double inf = 99999.999; // "infinity", don't use the FLOAT MAX because if PED calculation adds distance, program will crash!

// node for "tree" structure
typedef struct node_s{
	double PED;			// Partial Euclidean Distance of current node
	int  *PSV;			// Partial Symbol Vector, stores corresponding symbol indices of every layer
}node_t;

/*=============================Node Bubble Sort============================================================*\
\*=========================================================================================================*/
void node_bubble_sort(node_t* nodes, int length)
{
	int i,j;
	node_t tmp_node;

	for( j = length-1; j > 0; j-- )
	{
		for( i = 0; i < j; i++ )
		{
			if( nodes[i].PED > nodes[i+1].PED )
			{
				tmp_node = nodes[i];
				nodes[i] = nodes[i+1];
				nodes[i+1] = tmp_node;
			}
		}
	}
}

/*=============================Node Quick Sort============================================================*\
\*=========================================================================================================*/
void node_quick_sort(node_t* nodes, int left, int right){ //int arr[], int left, int right) {

      int i = left, j = right;
      node_t tmp_node;
      node_t pivot = nodes[ (int)((left + right) / 2) ];


      /* partition */
      while (i <= j) {
            while (nodes[i].PED < pivot.PED)
                  i++;
            while (nodes[j].PED > pivot.PED)
                  j--;
            if (i <= j) {
                  tmp_node = nodes[i];
                  nodes[i] = nodes[j];
                  nodes[j] = tmp_node;
                  i++;
                  j--;
            }
      };

      /* recursion */
      if (left < j)
            node_quick_sort(nodes, left, j);

      if (i < right)
            node_quick_sort(nodes, i, right);
}

/*=============================Random Number Generator=====================================================*\
\*=========================================================================================================*/
int get_rand(int lower_bound, int upper_bound)
{
	int tmp;
	if( lower_bound > upper_bound )
	{	// user input wrong, swap lower and upper bound...
		tmp = lower_bound;
		lower_bound = upper_bound;
		upper_bound = tmp;
	}
	return (int)( rand()%(upper_bound-lower_bound+1)+lower_bound );
}

/*=============================Partial Euclidean Distance Calculation======================================*\
\*=========================================================================================================*/
// calculate partial Euclidean distance of node (depends on parent node and previous partial symbol vectors!)
void calc_PED(node_t *node, node_t *parent,// pointers to current node and its parent
			 double *R_re, double *R_im,	 // upper triangular [nT x nT] matrix
			 double *y_re, double *y_im,  // received symbol column vector
			 double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
			 int current_symbol, int max_symbols, // current symbol of current layer, max number of symbols
			 int current_layer,  int max_layer, bool prune)
{
	int	jj, rr, symbol_index;				// count variables for loops
	double	distance_re = 0, distance_im = 0, symbol_re, symbol_im;	// for calculation of DI and PED

	// arrange partial symbol vector indices
	node->PSV[current_layer] = current_symbol;
	for(jj = current_layer+1; jj < max_layer; jj++ ) // PSV (partial symbol vector) has more than one element, get PSV from parent
		node->PSV[jj] = parent->PSV[jj];

	if( prune )	// this node shall be pruned
	{
		node->PED = inf;
	}
	else
	{
		// calculate distance increment of actual layer and current symbol
		for( jj = max_layer-1; jj >= current_layer; jj-- )
		{
			rr = current_layer+jj*max_layer;
			
			symbol_index =  max_symbols*jj + node->PSV[jj];
			symbol_re = symbol_alphabet_re[ symbol_index ];
			symbol_im = symbol_alphabet_im[ symbol_index ];
			
			distance_re += R_re[rr]*symbol_re - R_im[rr]*symbol_im;
			distance_im += R_re[rr]*symbol_im + R_im[rr]*symbol_re;		
		}
		
		distance_re = y_re[current_layer]-distance_re;
		distance_im = y_im[current_layer]-distance_im;
		
		// current PED:
		if( current_layer == max_layer-1 )
			node->PED = distance_re*distance_re + distance_im*distance_im;	// top layer, no parents have been assigned
		else
		{
			node->PED = parent->PED + distance_re*distance_re + distance_im*distance_im; // deeper layer with parents	
			if(node->PED > inf)
				node->PED = inf;
		}
	}
}

/*=============================Initial Search Radius Calculation 2=========================================*\
\*=========================================================================================================*/
// performs one smallest PED first search until leaf node reached, leaf PED is first radius
// = zero forcing solution
double get_ZF_distance(	
								double *R_re, double *R_im,	// upper triangular [nT x nT] matrix
								double *y_re, double *y_im, // received symbol column vector
								double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
								int *num_symbols, int max_symbols, int max_layer, node_t **node, bool memory_optimized )
{
	// local variables:
	node_t	*node_ptr, *parent_ptr;	// pointer to actual node and its parent
	int	i_layer, i_symbol;		// count variables for loops
	
	// depth first search (branch only to smallest PED until leaf reached)
	for( i_layer = max_layer-1; i_layer >= 0; i_layer--) // iterate through layers, begin at top
	{
		for( i_symbol = 0; i_symbol < num_symbols[i_layer]; i_symbol++ ) // iterate through symbols
		{
			if(memory_optimized)
			{
				node_ptr 	= &node[i_layer%2][i_symbol];		// point to current node
				if(i_layer < max_layer-1)
					parent_ptr 	= &node[(i_layer+1)%2][0]; 			// point to parent (0 because smallest PED has index 0, bubble sort...)			
				else
					parent_ptr = NULL;
			}
			else
			{
				node_ptr 	= &node[i_layer][i_symbol];		// point to current node
				if(i_layer < max_layer-1)
					parent_ptr 	= &node[i_layer+1][0]; 			// point to parent (0 because smallest PED has index 0, bubble sort...)			
				else
					parent_ptr = NULL;
			}
		

			// calculate partial Euclidean distance of current node
			calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
						i_symbol, max_symbols, i_layer, max_layer, false );
		}// symbol iteration
		
		// order nodes (symbols) with ascending PEDs
		if( memory_optimized )
			node_bubble_sort( node[i_layer%2], num_symbols[i_layer] );	
		else
			node_bubble_sort( node[i_layer], num_symbols[i_layer] );

	}// layer iteration		

	return node[0][0].PED; // radius equals min_PED of first leaf node reached...
}

/*=============================Initial Search Radius Calculation for counter hypotheses=========================================*\
\*=========================================================================================================*/
// performs one smallest PED first search until leaf node reached, leaf PED is first radius
// = zero forcing solution
double get_ZF_distance_ch(	
								double *R_re, double *R_im,	// upper triangular [nT x nT] matrix
								double *y_re, double *y_im, // received symbol column vector
								double *symbol_alphabet_re, double *symbol_alphabet_im, // symbol alphabet [max_nr_symbols x nT]
								int *num_symbols, int max_symbols, int max_layer, node_t **node,
								int j, node_t *ML_node, int *num_bits, bool *bittable, int total_bits, int current_bit, bool memory_optimized)
{
	// local variables:
	node_t	*node_ptr, *parent_ptr;	// pointer to actual node and its parent
	int	i_layer, i_symbol;		// count variables for loops
	
	// depth first search (branch only to smallest PED until leaf reached)
	for( i_layer = max_layer-1; i_layer >= 0; i_layer--) // iterate through layers, begin at top
	{
		for( i_symbol = 0; i_symbol < num_symbols[i_layer]; i_symbol++ ) // iterate through symbols
		{
			if(memory_optimized)
			{
				node_ptr 	= &node[i_layer%2][i_symbol];		// point to current node
				if(i_layer < max_layer-1)
					parent_ptr 	= &node[(i_layer+1)%2][0]; 			// point to parent (0 because smallest PED has index 0, bubble sort...)			
				else
					parent_ptr = NULL;
			}
			else
			{
				node_ptr 	= &node[i_layer][i_symbol];		// point to current node
				if(i_layer < max_layer-1)
					parent_ptr 	= &node[i_layer+1][0]; 			// point to parent (0 because smallest PED has index 0, bubble sort...)			
				else
					parent_ptr = NULL;
			}
				
			if ( ( i_layer == j ) && ( bittable[ current_bit + total_bits*i_symbol ] == bittable[ current_bit + total_bits*ML_node->PSV[i_layer] ] ) )
			{ 	// node is pruned
					//node_ptr->PED = inf; // this will prune ML subset
					calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
								i_symbol, max_symbols, i_layer, max_layer, true );
			}
			else
			{
				// calculate partial Euclidean distance of current node
					calc_PED(	node_ptr, parent_ptr, R_re, R_im, y_re, y_im, symbol_alphabet_re, symbol_alphabet_im,
								i_symbol, max_symbols, i_layer, max_layer, false );
			}
		}// symbol iteration
		
		// order nodes (symbols) with ascending PEDs
		if( memory_optimized )
			node_bubble_sort( node[i_layer%2], num_symbols[i_layer] );	
		else
			node_bubble_sort( node[i_layer], num_symbols[i_layer] );
	}// layer iteration		

	return node[0][0].PED; // radius equals min_PED of first leaf node reached...
}

#endif // SSD_FUNCTIONS_H_
