
#include "kernel_qr_inverse.hpp"
#include "matrix_multiply.hpp"
#include "x_matrix_utils.hpp"
#include "xf_solver_L1.hpp"

using namespace xf::solver;

extern "C" int func_matrixMultiply(hls::stream<MATRIX_IN_T>& matrixAStrm,
                                   hls::stream<MATRIX_OUT_T>& matrixInverseAStrm) {

	int is_singular = 0;

	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixInverseAStrm);
    return is_singular;
}

extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& matrixAStrm,
                                   hls::stream<MATRIX_OUT_T>& matrixInverseAStrm) {
	int is_singular = 0;
    xf::solver::qrInverse<ROWSCOLSA, MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixInverseAStrm, is_singular);
    return is_singular;
}

extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& matrixAStrm,
                                   hls::stream<MATRIX_OUT_T>& matrixInverseAStrm) {
	int is_singular = 0;
    xf::solver::qrInverse<ROWSCOLSA, MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixInverseAStrm, is_singular);
    return is_singular;
}


void
func_matrixAddition(hls::stream<MATRIX_IN_T>& matrixAStrm,
					hls::stream<MATRIX_IN_T>& matrixBStrm,
		            hls::stream<MATRIX_OUT_T>& matrixAddABStrm)
{
	MATRIX_IN_T A[ROWSCOLSA][ROWSCOLSA];
	MATRIX_IN_T B[ROWSCOLSA][ROWSCOLSA];
	MATRIX_OUT_T Out1[ROWSCOLSA][ROWSCOLSA];

    #pragma HLS inline
    //#pragma HLS stream variable=Constant_out1 depth=ROWS_1
	matrixAddition_assign_Ri:
 	for (int i = 0; i < ROWSCOLSA; i++) {
 	    for (int j = 0; j < ROWSCOLSA; j++) {
 	    		//#pragma HLS PIPELINE
            A[i][j] = matrixAStrm.read();
            B[i][j] = matrixBStrm.read();
 	    	Out1[i][j] = A[i][j]+B[i][j];
 	    	matrixAddABStrm.write(Out1[i][j]);
 	    }
    }

}




// 2 Matrices

//extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& matrixAStrm,
//									hls::stream<MATRIX_IN_T>& matrixBStrm,
//                                   hls::stream<MATRIX_OUT_T>& matrixMultABStrm) {
//
//	matrixMultiply<Transpose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
//	MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixBStrm,matrixMultABStrm);
//    return 0;
//}
