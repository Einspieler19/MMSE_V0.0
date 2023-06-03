
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



void
func_matrixAddition2(hls::stream<MATRIX_IN_T>& matrixAStrm,
		            hls::stream<MATRIX_OUT_T>& matrixAddABStrm,
					float var_Noise)
{
	MATRIX_IN_T A[ROWSCOLSA][ROWSCOLSA];
	MATRIX_OUT_T Out1[ROWSCOLSA][ROWSCOLSA];
	MATRIX_IN_T var_Noise_complex;
	var_Noise_complex.real(var_Noise);
	var_Noise_complex.imag(0.0);

    #pragma HLS inline
    //#pragma HLS stream variable=Constant_out1 depth=ROWS_1
	matrixAddition_assign_Ri:
 	for (int i = 0; i < ROWSCOLSA; i++) {
 	    for (int j = 0; j < ROWSCOLSA; j++) {
 	    		//#pragma HLS PIPELINE
            A[i][j] = matrixAStrm.read();
            if (i == j) {
//            	Out1[i][j] = A[i][j] + (var_Noise,0.0f);
            	Out1[i][j] = A[i][j] + var_Noise_complex;

            } else {
            	Out1[i][j] = A[i][j];
            }
//        	std::cout<<A[i][j]<<" "<<i << " " <<j<< std::endl;
//        	std::cout<<Out1[i][j]<<std::endl;
 	    	matrixAddABStrm.write(Out1[i][j]);

    }
}
}


extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& matrixAStrm,
								   hls::stream<MATRIX_OUT_T>& matrixMMSEH,
								   float var_Noise) {

	hls::stream<MATRIX_OUT_T> matrixA_Mul_ATStrm;
#pragma HLS STREAM variable = matrixA_Mul_ATStrm depth = 9


//	hls::stream<MATRIX_IN_T> matrixIStrm;
//#pragma HLS STREAM variable = matrixIStrm depth = 9

	hls::stream<MATRIX_OUT_T> matrixToInverse;
#pragma HLS STREAM variable = matrixToInverse depth = 9



//    func_matrixMultiply(matrixAStrm, matrixA_Mul_ATStrm);
	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(matrixAStrm, matrixA_Mul_ATStrm);



	// 单位矩阵
//
//    MATRIX_OUT_T I[ROWSCOLSA][ROWSCOLSA]; // The identity matrix to compare against
//	// Create I to compare against later
//	    for (int r = 0; r < ROWSCOLSA; r++) {
//	        for (int c = 0; c < ROWSCOLSA; c++) {
//	            if (r == c) {
//	                I[r][c].real(var_Noise);
//	                I[r][c].imag(0.0);
//	            } else {
//	                I[r][c] = 0.0;
//	            }
//	        }
//	    }
//
//        for (int r = 0; r < ROWSCOLSA; r++) {
//            for (int c = 0; c < ROWSCOLSA; c++) {
//                matrixIStrm.write(I[r][c]);
//            }
//        }
//	func_matrixAddition(matrixA_Mul_ATStrm,
//    					matrixIStrm,
//						matrixToInverse);


	func_matrixAddition2(matrixA_Mul_ATStrm,
							matrixToInverse,
							var_Noise);

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
