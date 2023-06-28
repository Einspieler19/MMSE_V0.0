
#include "kernel_qr_inverse.hpp"
#include "matrix_multiply.hpp"
#include "x_matrix_utils.hpp"
#include "xf_solver_L1.hpp"

using namespace xf::solver;

extern "C" int func_matrixMultiply(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
                                   hls::stream<MATRIX_OUT_T>& matrixInverseAStrm) {

	int is_singular = 0;

	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(in_channelMatrix_Strm, matrixInverseAStrm);
    return is_singular;
}

extern "C" int func_qr_inverse(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
                                   hls::stream<MATRIX_OUT_T>& matrixInverseAStrm) {
	int is_singular = 0;
    xf::solver::qrInverse<ROWSCOLSA, MATRIX_IN_T, MATRIX_OUT_T>(in_channelMatrix_Strm, matrixInverseAStrm, is_singular);
    return is_singular;
}


void
func_matrixAddition(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
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
            A[i][j] = in_channelMatrix_Strm.read();
            B[i][j] = matrixBStrm.read();
 	    	Out1[i][j] = A[i][j]+B[i][j];
 	    	matrixAddABStrm.write(Out1[i][j]);
 	    }
    }
}



void
func_matrixAddition2(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
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
            A[i][j] = in_channelMatrix_Strm.read();
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


extern "C" int compute_Weights(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
								   hls::stream<MATRIX_T>& in_weights_Strm,
								   float var_Noise) {
#pragma HLS DATAFLOW
	hls::stream<MATRIX_OUT_T> matrixA_Mul_ATStrm;
#pragma HLS STREAM variable = matrixA_Mul_ATStrm depth = ROWSCOLSA*ROWSCOLSA

	hls::stream<MATRIX_IN_T> matrixIStrm;
#pragma HLS STREAM variable = matrixIStrm depth = ROWSCOLSA*ROWSCOLSA

	hls::stream<MATRIX_OUT_T> matrixToInverse;
#pragma HLS STREAM variable = matrixToInverse depth = ROWSCOLSA*ROWSCOLSA

	hls::stream<MATRIX_OUT_T> inv_Strm;
#pragma HLS STREAM variable = inv_Strm depth = ROWSCOLSA*ROWSCOLSA

	hls::stream<MATRIX_OUT_T> ChannelMatrix_Strm;
#pragma HLS STREAM variable = ChannelMatrix_Strm depth = ROWSCOLSA*ROWSCOLSA
	hls::stream<MATRIX_OUT_T> ChannelMatrix_Strm2;
#pragma HLS STREAM variable = ChannelMatrix_Strm depth = ROWSCOLSA*ROWSCOLSA

	MATRIX_T ChannelMatrix[ROWSCOLSA][ROWSCOLSA];
			    for (int r = 0; r < ROWSCOLSA; r++) {
			        for (int c = 0; c < ROWSCOLSA; c++) {
			        	in_channelMatrix_Strm.read(ChannelMatrix[r][c]);
			        	ChannelMatrix_Strm.write(ChannelMatrix[r][c]);
			        }
			    }


//    func_matrixMultiply(in_channelMatrix_Strm, matrixA_Mul_ATStrm);
	matrixMultiply<NoTranspose, ConjugateTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(ChannelMatrix_Strm, matrixA_Mul_ATStrm);
//  H'*H

	// 单位矩阵

    MATRIX_OUT_T I[ROWSCOLSA][ROWSCOLSA]; // The identity matrix to compare against
	// Create I to compare against later
	    compute_Weights_label0:for (int r = 0; r < ROWSCOLSA; r++) {
#pragma HLS UNROLL
	        for (int c = 0; c < ROWSCOLSA; c++) {
	            if (r == c) {
	                I[r][c].real(var_Noise);
	                I[r][c].imag(0.0);
	            } else {
	                I[r][c] = 0.0;
	            }
	        }
	    }

        for (int r = 0; r < ROWSCOLSA; r++) {
#pragma HLS PIPELINE
            for (int c = 0; c < ROWSCOLSA; c++) {
                matrixIStrm.write(I[r][c]);
            }
        }

        //
	func_matrixAddition(matrixA_Mul_ATStrm,
    					matrixIStrm,
						matrixToInverse);
//  H'*H + sigmaI

//	func_matrixAddition2(matrixA_Mul_ATStrm,
//							matrixToInverse,
//							var_Noise);

	int is_singular = 0;
//    qr_inverse_return = func_qr_inverse(matrixToInverse, in_weights_Strm);
    xf::solver::qrInverse<ROWSCOLSA, MATRIX_IN_T, MATRIX_T>(matrixToInverse, inv_Strm, is_singular);

    for (int r = 0; r < ROWSCOLSA; r++) {
        for (int c = 0; c < ROWSCOLSA; c++) {
        	ChannelMatrix_Strm2.write(ChannelMatrix[r][c]);
        }
    }

//	matrixMultiply<NoTranspose, NoTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
//	MATRIX_IN_T, MATRIX_OUT_T>(inv_Strm, ChannelMatrix_Strm, in_weights_Strm);

	matrixMultiply<ConjugateTranspose, NoTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA,
	MATRIX_IN_T, MATRIX_OUT_T>(ChannelMatrix_Strm2, inv_Strm,  in_weights_Strm);


    return is_singular;

}

void my_Precoder(hls::stream<MATRIX_T>& in_weights_Strm,
									hls::stream<MATRIX_IN_T>& in_BitStrm,
								   hls::stream<MATRIX_OUT_T>& out_precodedMatrixStrm) {
	matrixMultiply<NoTranspose, NoTranspose, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, N_FRAMEperOP, ROWSCOLSA, N_FRAMEperOP,
	MATRIX_IN_T, MATRIX_OUT_T>(in_weights_Strm, in_BitStrm, out_precodedMatrixStrm);
}



extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
									hls::stream<MATRIX_IN_T>& in_BitStrm,
									hls::stream<MATRIX_T>& in_weights_Strm,
									hls::stream<MATRIX_T>& out_weightsUpdate_Strm,
								   hls::stream<MATRIX_OUT_T>& out_precodedMatrixStrm,
								   float var_Noise,
								   int flag) {
#pragma HLS DATAFLOW
#pragma HLS STREAM variable = in_weights_Strm depth = ROWSCOLSA*ROWSCOLSA

	const MATRIX_T Weights[ROWSCOLSA][ROWSCOLSA];          // The inverse result from the DUT


	// update W
	if(flag==1){
		int qr_inverse_return = 0;
			qr_inverse_return = compute_Weights(in_channelMatrix_Strm, out_weightsUpdate_Strm, var_Noise);
//		    for (int r = 0; r < ROWSCOLSA; r++) {
//		        for (int c = 0; c < ROWSCOLSA; c++) {
//		        	in_weights_Strm.read(Weights[r][c]);
//		        	out_weightsUpdate_Strm.write(Weights[r][c]);
//		        }
//		    }
	}
	else{
//		for (int r = 0; r < ROWSCOLSA; r++) {
//	        for (int c = 0; c < ROWSCOLSA; c++) {
//	        	in_weights_Strm.write(Weights[r][c]);
//	        }
//	    }
		my_Precoder(in_weights_Strm, in_BitStrm, out_precodedMatrixStrm);
	}
	return 0;
}

