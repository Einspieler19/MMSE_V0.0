
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

extern "C" int func_qr_inverse(hls::stream<MATRIX_IN_T>& matrixAStrm,
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



// 只有mult

extern "C" int kernel_qr_inverse_0a(hls::stream<MATRIX_IN_T>& matrixAStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH) {

//    func_matrixMultiply(matrixAStrm, matrixMMSEH);
	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixMMSEH);
//    func_matrixAddition(matrixA_Mul_ATStrm,
//    					matrixIStrm,
//						matrixMMSEH);

	int qr_inverse_return = 0;
//    qr_inverse_return = func_qr_inverse(matrixA_Mul_ATStrm, matrixMMSEH);
    return qr_inverse_return;
}



//  只有inverse

extern "C" int kernel_qr_inverse_00(hls::stream<MATRIX_IN_T>& matrixAStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH) {


	int qr_inverse_return = 0;
    qr_inverse_return = func_qr_inverse(matrixAStrm, matrixMMSEH);
    return qr_inverse_return;
}

// inverse 和 add

extern "C" int kernel_qr_inverse_01(hls::stream<MATRIX_IN_T>& matrixAStrm,
									hls::stream<MATRIX_IN_T>& matrixIStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH) {

//	hls::stream<MATRIX_OUT_T> matrixA_Mul_ATStrm;
	hls::stream<MATRIX_OUT_T> matrixToInverse;


//    func_matrixMultiply(matrixAStrm, matrixA_Mul_ATStrm);

    func_matrixAddition(matrixAStrm,
    					matrixIStrm,
						matrixToInverse);

	int qr_inverse_return = 0;
    qr_inverse_return = func_qr_inverse(matrixToInverse, matrixMMSEH);
    return qr_inverse_return;
}



// inverse 和 mult

extern "C" int kernel_qr_inverse_0b(hls::stream<MATRIX_IN_T>& matrixAStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH) {

#pragma HLS DATAFLOW

	hls::stream<MATRIX_OUT_T> matrixA_Mul_ATStrm;
#pragma HLS STREAM variable = matrixA_Mul_ATStrm depth = 9


	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixA_Mul_ATStrm);
//    func_matrixMultiply(matrixAStrm, matrixA_Mul_ATStrm);

//    func_matrixAddition(matrixA_Mul_ATStrm,
//    					matrixIStrm,
//						matrixMMSEH);

	int is_singular = 0;
    xf::solver::qrInverse<ROWSCOLSA, MATRIX_IN_T, MATRIX_OUT_T>(matrixA_Mul_ATStrm, matrixMMSEH, is_singular);
    return is_singular;
}


// add 和 mult

extern "C" int kernel_qr_inverse_2(hls::stream<MATRIX_IN_T>& matrixAStrm,
									hls::stream<MATRIX_IN_T>& matrixIStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH) {

	hls::stream<MATRIX_OUT_T> matrixA_Mul_ATStrm;


    func_matrixMultiply(matrixAStrm, matrixA_Mul_ATStrm);

    func_matrixAddition(matrixA_Mul_ATStrm,
    					matrixIStrm,
						matrixMMSEH);

	int qr_inverse_return = 0;
//    qr_inverse_return = func_qr_inverse(matrixToInverse, matrixMMSEH);
    return qr_inverse_return;
}


//  都有

extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& matrixAStrm,
									hls::stream<MATRIX_IN_T>& matrixIStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH) {

	hls::stream<MATRIX_OUT_T> matrixA_Mul_ATStrm;
#pragma HLS STREAM variable = matrixA_Mul_ATStrm depth = 9
	hls::stream<MATRIX_OUT_T> matrixToInverse;
#pragma HLS STREAM variable = matrixToInverse depth = 9

//    func_matrixMultiply(matrixAStrm, matrixA_Mul_ATStrm);
	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixA_Mul_ATStrm);

    func_matrixAddition(matrixA_Mul_ATStrm,
    					matrixIStrm,
						matrixToInverse);

	int is_singular = 0;
//    qr_inverse_return = func_qr_inverse(matrixToInverse, matrixMMSEH);
    xf::solver::qrInverse<ROWSCOLSA, MATRIX_IN_T, MATRIX_OUT_T>(matrixToInverse, matrixMMSEH, is_singular);
    return is_singular;

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
