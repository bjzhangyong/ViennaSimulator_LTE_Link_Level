/**
 * \file crc24b.c
 * Functions and types for CRC checks.
 * MEXed to LTE_common_crc24b
 *
 * Generated on Wed Jun 04 13:03:52 2008,
 * by pycrc v0.6.5, http://www.tty1.net/pycrc/
 * MEX-ed by Josep Colom Ikuno, TU Wien, INTHFT
 * jcolom@nt.tuwien.ac.at
 * www.nt.tuwien.ac.at
 * using the configuration:
 *    Width        = 24
 *    Poly         = 0x800063
 *    XorIn        = 0x000000
 *    ReflectIn    = False
 *    XorOut       = 0x000000
 *    ReflectOut   = False
 *    Algorithm    = table-driven
 *
 * gCRC24B= [D24 + D23 + D6 + D5 + D + 1] for a CRC length L = 24
 * poly = dec2hex(binvec2dec([1 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]))
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


 *****************************************************************************/
#include "crc24b.h"
#include "mex.h"
#include <stdlib.h>

/**
 * Static table used for the table_driven implementation.
 *****************************************************************************/
static const crc_t crc_table[256] = {
    0x000000, 0x800063, 0x8000a5, 0x0000c6, 0x800129, 0x00014a, 0x00018c, 0x8001ef,
    0x800231, 0x000252, 0x000294, 0x8002f7, 0x000318, 0x80037b, 0x8003bd, 0x0003de,
    0x800401, 0x000462, 0x0004a4, 0x8004c7, 0x000528, 0x80054b, 0x80058d, 0x0005ee,
    0x000630, 0x800653, 0x800695, 0x0006f6, 0x800719, 0x00077a, 0x0007bc, 0x8007df,
    0x800861, 0x000802, 0x0008c4, 0x8008a7, 0x000948, 0x80092b, 0x8009ed, 0x00098e,
    0x000a50, 0x800a33, 0x800af5, 0x000a96, 0x800b79, 0x000b1a, 0x000bdc, 0x800bbf,
    0x000c60, 0x800c03, 0x800cc5, 0x000ca6, 0x800d49, 0x000d2a, 0x000dec, 0x800d8f,
    0x800e51, 0x000e32, 0x000ef4, 0x800e97, 0x000f78, 0x800f1b, 0x800fdd, 0x000fbe,
    0x8010a1, 0x0010c2, 0x001004, 0x801067, 0x001188, 0x8011eb, 0x80112d, 0x00114e,
    0x001290, 0x8012f3, 0x801235, 0x001256, 0x8013b9, 0x0013da, 0x00131c, 0x80137f,
    0x0014a0, 0x8014c3, 0x801405, 0x001466, 0x801589, 0x0015ea, 0x00152c, 0x80154f,
    0x801691, 0x0016f2, 0x001634, 0x801657, 0x0017b8, 0x8017db, 0x80171d, 0x00177e,
    0x0018c0, 0x8018a3, 0x801865, 0x001806, 0x8019e9, 0x00198a, 0x00194c, 0x80192f,
    0x801af1, 0x001a92, 0x001a54, 0x801a37, 0x001bd8, 0x801bbb, 0x801b7d, 0x001b1e,
    0x801cc1, 0x001ca2, 0x001c64, 0x801c07, 0x001de8, 0x801d8b, 0x801d4d, 0x001d2e,
    0x001ef0, 0x801e93, 0x801e55, 0x001e36, 0x801fd9, 0x001fba, 0x001f7c, 0x801f1f,
    0x802121, 0x002142, 0x002184, 0x8021e7, 0x002008, 0x80206b, 0x8020ad, 0x0020ce,
    0x002310, 0x802373, 0x8023b5, 0x0023d6, 0x802239, 0x00225a, 0x00229c, 0x8022ff,
    0x002520, 0x802543, 0x802585, 0x0025e6, 0x802409, 0x00246a, 0x0024ac, 0x8024cf,
    0x802711, 0x002772, 0x0027b4, 0x8027d7, 0x002638, 0x80265b, 0x80269d, 0x0026fe,
    0x002940, 0x802923, 0x8029e5, 0x002986, 0x802869, 0x00280a, 0x0028cc, 0x8028af,
    0x802b71, 0x002b12, 0x002bd4, 0x802bb7, 0x002a58, 0x802a3b, 0x802afd, 0x002a9e,
    0x802d41, 0x002d22, 0x002de4, 0x802d87, 0x002c68, 0x802c0b, 0x802ccd, 0x002cae,
    0x002f70, 0x802f13, 0x802fd5, 0x002fb6, 0x802e59, 0x002e3a, 0x002efc, 0x802e9f,
    0x003180, 0x8031e3, 0x803125, 0x003146, 0x8030a9, 0x0030ca, 0x00300c, 0x80306f,
    0x8033b1, 0x0033d2, 0x003314, 0x803377, 0x003298, 0x8032fb, 0x80323d, 0x00325e,
    0x803581, 0x0035e2, 0x003524, 0x803547, 0x0034a8, 0x8034cb, 0x80340d, 0x00346e,
    0x0037b0, 0x8037d3, 0x803715, 0x003776, 0x803699, 0x0036fa, 0x00363c, 0x80365f,
    0x8039e1, 0x003982, 0x003944, 0x803927, 0x0038c8, 0x8038ab, 0x80386d, 0x00380e,
    0x003bd0, 0x803bb3, 0x803b75, 0x003b16, 0x803af9, 0x003a9a, 0x003a5c, 0x803a3f,
    0x003de0, 0x803d83, 0x803d45, 0x003d26, 0x803cc9, 0x003caa, 0x003c6c, 0x803c0f,
    0x803fd1, 0x003fb2, 0x003f74, 0x803f17, 0x003ef8, 0x803e9b, 0x803e5d, 0x003e3e
};



/**
 * Update the crc value with new data.
 *
 * \param crc      The current crc value.
 * \param data     Pointer to a buffer of \a data_len bytes.
 * \param data_len Number of bytes in the \a data buffer.
 * \return         The updated crc value.
 *****************************************************************************/
crc_t crc_update(crc_t crc, const unsigned char *data, size_t data_len)
{
    unsigned int tbl_idx;

    while (data_len--) {
        tbl_idx = ((crc >> 16) ^ *data) & 0xff;
        crc = (crc_table[tbl_idx] ^ (crc << 8)) & 0xffffff;

        data++;
    }
    return crc & 0xffffff;
}




/*#include <stdio.h> */
/*#define DEBUG */

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    crc_t crc;
    unsigned char *msg;
    mwSize input_length;
    mwSize i;
    unsigned char *output_bytes;
    const mwSize dims[]={1,3};
    
    /* check for proper number of arguments */
    if(nrhs!=1)
        mexErrMsgTxt("One input required.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    
    /* input must be a string */
    if ( mxIsUint8(prhs[0]) != 1)
        mexErrMsgTxt("Input must be uint8.");
    
    /* Get the input size */
    input_length = mxGetN(prhs[0]); /* The number of columns */
    
    /* Get the input data */
    msg = (unsigned char **)mxGetPr(prhs[0]);
    
    /* call the C subroutine */
    crc = crc_init();
    crc = crc_update(crc, msg, input_length);
    crc = crc_finalize(crc);
    
    /* Return output */
    plhs[0] = mxCreateNumericArray(2,dims,mxUINT8_CLASS,mxREAL);
    output_bytes    = (unsigned char *) mxGetData(plhs[0]);
    output_bytes[0] = (unsigned char)(crc  & 0x000000ff);
    output_bytes[1] = (unsigned char)((crc & 0x0000ff00) >> 8);
    output_bytes[2] = (unsigned char)((crc & 0x00ff0000) >> 16);
    
    #ifdef DEBUG
    /* Print some debug output */
    mexPrintf("LTE CRC24B\n");
    mexPrintf("Message: 0x");
    for(i=0;i<input_length;i++) {
        mexPrintf("%2x",msg[i]);
    }
    mexPrintf("\n");
    mexPrintf("Length:  %d\n",input_length);
    mexPrintf("Poly:    %#x\n",crc_table[1]);
    mexPrintf("CRC:     %#x\n",crc);
    mexPrintf("CRC:     '%#x' '%x' '%x'\n",output_bytes[0],output_bytes[1],output_bytes[2]);
    mexPrintf("CRC:     '%d' '%d' '%d'\n",output_bytes[0],output_bytes[1],output_bytes[2]);
    #endif
}
