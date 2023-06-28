#include "test_qr_inverse.hpp"
#include "kernel_qr_inverse.hpp"

#include "src/utils.hpp"
#include "hw/utils/x_matrix_utils.hpp"
#include "src/matrix_test_utils.hpp"
#include "src/type_test_utils.hpp"

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

// ---------------------------------------------------------------------------------------------
// Main test program
// ---------------------------------------------------------------------------------------------


int main() {

	int flag[OP_LENGTH];
	int Index_FRAME = 0;


	long unsigned num_tests = 0;
	// 每类矩阵有多少个Test case?
	// long unsigned num_tests = (ROWSCOLSA >= 16 ? 5 : 20); // Small default for HLS
	unsigned int debug = 1;
	// 是否输出矩阵

	double ratio_threshold = 30.0;
	// ratio多少算是差别太大

	double mat_type = 13; // Specify to only run a single matrix type.
	// 0:全测试
	// n:第n个测试

	unsigned print_precision = 10;


	unsigned allowed_ulp_mismatch = 0;
	//调用matrices equal的参数

	// 变量定义
	// ====================================================================
	int qr_inverse_return = 0; // Return code from hls::qr_inverse
	// 求逆成功了吗


	// Matrix arrays
	MATRIX_IN_T A[ROWSCOLSA][ROWSCOLSA]; // The Channel Matrix.  Cast from A_generated

	MATRIX_OUT_T Weights[ROWSCOLSA][ROWSCOLSA];          // The inverse result from the DUT
	MATRIX_OUT_T Weights_expected[ROWSCOLSA][ROWSCOLSA]; // The inverse result from REF in target format

	MATRIX_OUT_T PrecodedMMSE[ROWSCOLSA][N_FRAMEperOP]; // The inverse result from REF in target format
	// Test variables
	QR_INV_TYPE A_cast[ROWSCOLSA][ROWSCOLSA]; // Cast A back to LAPACK compatible type for analysis
	// A切换到另一种数据类型，用于求A的norm

	MATRIX_IN_T in_Data_all[ROWSCOLSA][FRAME_LENGTH]; // The input data array.
	MATRIX_IN_T in_Data[ROWSCOLSA][N_FRAMEperOP]; // The input data array.

	QR_INV_TYPE I_delta[ROWSCOLSA][ROWSCOLSA];
	// 差矩阵1

	MATRIX_OUT_T Weights_delta[ROWSCOLSA]
										 [ROWSCOLSA]; // The difference between the DUT's Weights and LAPACK's Weights
	// 差矩阵2


    MATRIX_T I[ROWSCOLSA][ROWSCOLSA]; // The identity matrix
	// Create I to compare against later
	    for (int r = 0; r < ROWSCOLSA; r++) {
	        for (int c = 0; c < ROWSCOLSA; c++) {
	            if (r == c) {
	                I[r][c].real(1.0);
	                I[r][c].imag(0.0);
	            } else {
	                I[r][c] = 0.0;
	            }
                Weights[r][c] = I[r][c];
	        }
	    }


	// Non-complex type
	double A_norm;
	double Weights_norm;
	double I_delta_norm;

	unsigned int pass_fail = 0; // Pass=0 Fail =1
	// 程序总结果

	bool matched_expected_Weights;
	// 是否一致

	QR_INV_BASE_TYPE I_DUT_ratio;
	// 有多不一致


	// 检测输入数据是什么数据类型并打印
	// 好像不影响程序 没什么卵用
	printf("Running %lu %s tests per matrix type on %d x %d matrices\n", num_tests,
			x_is_float(Weights_expected[0][0])
			? "single precision"
					: x_is_double(Weights_expected[0][0])
					  ? "double precision"
							  : x_is_fixed(Weights_expected[0][0]) ? "fixed point" : "Unknown type",
									  ROWSCOLSA, ROWSCOLSA);



	// Read input matrix and golden matrix from files
	// 找根路径
	std::string data_path = std::string(DATA_PATH);
	std::string base_path;


	//    Weights
	//====================================================================
	// 是哪个数据类型的，从哪读数据
	if (x_is_complex(A[0][0]) == true) {
		base_path = data_path.append("/complex/");
	} else {
		base_path = data_path.append("/float/");
	}

	std::string file_A =
			base_path + "A_matType_" + std::to_string(int(mat_type)) + "_" + std::to_string(num_tests) + ".txt";
	std::string file_Weights_expected =
			base_path + "TestPoint1_A_matType_" + std::to_string(int(mat_type)) + "_" + std::to_string(num_tests) + ".txt";


	// 文件size, 读取的参数
	int A_size = ROWSCOLSA * ROWSCOLSA;
	int Weights_expected_size = ROWSCOLSA * ROWSCOLSA;

	// 数组指针, 读取的参数
	MATRIX_IN_T* A_ptr = reinterpret_cast<MATRIX_IN_T*>(A);
	MATRIX_OUT_T* Weights_expected_ptr = reinterpret_cast<MATRIX_OUT_T*>(Weights_expected);


	//    Input Bitstream
	//====================================================================
	std::string base_path_in_Data_all;
	std::string base_path_flag;
	base_path_in_Data_all = data_path.append("/");
	base_path_flag = data_path.append("/");

	std::string file_in_Data_all =
			base_path_in_Data_all + "4_10_QAM_Input_BitStream" + ".dat";
	std::string file_flag =
			base_path_flag + "W_update_flag" + ".dat";


	// 文件size, 读取的参数
	int in_Data_all_size = ROWSCOLSA*FRAME_LENGTH ;
	int flag_size = OP_LENGTH;


	// 数组指针, 读取的参数
	int* flag_ptr = flag;
	MATRIX_IN_T* in_Data_all_ptr = reinterpret_cast<MATRIX_IN_T*>(in_Data_all);


	//    Stream
	//====================================================================
	// 进出DUT的流
	hls::stream<MATRIX_IN_T> matrixAStrm;
	hls::stream<MATRIX_IN_T> in_Data_alltrm;
	hls::stream<MATRIX_IN_T> in_WeightsStrm;
	hls::stream<MATRIX_OUT_T> out_WeightsStrm;
	hls::stream<MATRIX_OUT_T> matrixMMSEH;


	readTxt(file_flag, flag_ptr, flag_size);


	// 读取
	readTxt(file_in_Data_all, in_Data_all_ptr, in_Data_all_size);


	//    Main loop
	//====================================================================
	for (unsigned int Index_OP = 0; Index_OP < OP_LENGTH;) {

		if(flag[Index_OP]==1){
			// 读取
			readTxt(file_A, A_ptr, A_size);
			readTxt(file_Weights_expected, Weights_expected_ptr, Weights_expected_size);

			// 写入流
			for (int r = 0; r < ROWSCOLSA; r++) {
				for (int c = 0; c < ROWSCOLSA; c++) {
					matrixAStrm.write(A[r][c]);
				}
			}
			// cast
			for (int r = 0; r < ROWSCOLSA; r++) {
				for (int c = 0; c < ROWSCOLSA; c++) {
					A_cast[r][c] = A[r][c];
				}
			}
		}

		else{

			for (int r = 0; r < ROWSCOLSA; r++) {
							for (int c = 0; c < N_FRAMEperOP; c++) {
								in_Data[r][c] = in_Data_all[r][Index_FRAME];
							}
						}

			// 写入流
			for (int r = 0; r < ROWSCOLSA; r++) {
				for (int c = 0; c < N_FRAMEperOP; c++) {
					in_Data_alltrm.write(in_Data[r][c]);
				}
			}
			for (int r = 0; r < ROWSCOLSA; r++) {
				for (int c = 0; c < ROWSCOLSA; c++) {
					in_WeightsStrm.write(Weights[r][c]);
				}
			}
		}


		// Noise, set with random values later
		float var_Noise = 1;

		// run DUT
		qr_inverse_return = kernel_qr_inverse_0(matrixAStrm, in_Data_alltrm, in_WeightsStrm, out_WeightsStrm, matrixMMSEH, var_Noise, flag[Index_OP]);



		if(flag[Index_OP]==1){


			// Check if inverse has ran successfully

						if (qr_inverse_return == 1) {
							printf("ERROR: Inverse failed\n");
							xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
									Weights, "   ", print_precision, 0);
							printf("TB:Fail\n");
							return (1);
						}


			// 读出流
			for (int r = 0; r < ROWSCOLSA; r++) {
				for (int c = 0; c < ROWSCOLSA; c++) {
					out_WeightsStrm.read(Weights[r][c]);
				}
			}


			// Check for NaNs in result

			if (anyNaN<ROWSCOLSA, ROWSCOLSA>(Weights) == 1) {
				printf("ERROR: Caught NaN in Weights\n");
				xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
						Weights, "   ", print_precision, 0);
				printf("TB:Fail\n");
				return (32);
			}


			// Test results
			// ====================================================================

			// Basic check cell by cell check based on allowed_ulp_mismatch value.
			// 求差方法1
			// 用一个函数计算与目标结果的差
			msub<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, QR_INV_TYPE>(Weights, Weights_expected, I_delta);

			// 求差方法2
			// 用一个函数计算与目标结果是否一致，差值多少
			matched_expected_Weights = are_matrices_equal<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T>(
					(MATRIX_OUT_T*)Weights, (MATRIX_OUT_T*)Weights_expected, allowed_ulp_mismatch,
					(MATRIX_OUT_T*)Weights_delta);

			// 求norm
			A_norm = norm1_dbl<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(A_cast);
			Weights_norm = norm1_dbl<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(Weights_expected);
			I_delta_norm = norm1<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(I_delta);

			// 验证norm的异常值

			if (isinf(A_norm)) {
				// Should never be Inf - if it is, we've probably overflowed
				printf("ERROR: Caught unexpected Inf for A_norm\n");
				printf("TB:Fail\n");
				return (4);
			}

			if (isinf(Weights_norm)) {
				// Should never be Inf - if it is, we've probably overflowed
				printf("ERROR: Caught unexpected Inf for Weights_norm\n");
				printf("TB:Fail\n");
				return (5);
			}


			// norm的比率
			I_DUT_ratio =(double)I_delta_norm / (double)A_norm;

			// Check that the norm values are OK and we are not comparing two bad ratios
			// 检验I_DUT_ratio异常值

			if (isnan(I_DUT_ratio)) {
				// Should only be NaN if A_norm was zero, so check that
				if (A_norm != 0) {
					printf("ERROR: Caught unexpected NaN for I_DUT_ratio\n");
					printf("TB:Fail\n");
					return (6);
				}
			}

			if (isinf(I_DUT_ratio)) {
				// Should never be Inf
				printf("ERROR: Caught unexpected Inf for I_DUT_ratio\n");
				printf("TB:Fail\n");
				return (8);
			}

			if (I_DUT_ratio < 0) {
				// Should never be less than zero - if it is, it's either an error code or something went badly
				// wrong
				printf("ERROR: Caught unexpected negative I_DUT_ratio\n");
				printf("TB:Fail\n");
				return (12);
			}


			// Determine if pass or fail.
			// o Check DUT ratio against test threshold, default taken from LAPACK
			if (I_DUT_ratio > ratio_threshold ) {
				std::cout << "ERROR: I_DUT_ratio(" << I_DUT_ratio << ") > ratio_threshold(" << ratio_threshold
						<< "). " << std::endl;
				pass_fail = 1; // Test run fails
			}
		}

		else{
			// 读出流
			for (int r = 0; r < ROWSCOLSA; r++) {
				for (int c = 0; c < N_FRAMEperOP; c++) {
					matrixMMSEH.read(PrecodedMMSE[r][c]);
				}
			}
		}



	// Print matrices for debug
	// 输出所有矩阵
	if ( debug > 0 || I_DUT_ratio > ratio_threshold) {
		std::cout<< "  flag = "<< flag[Index_OP] <<std::endl;
		std::cout<< "  Index_OP = "<< Index_OP <<std::endl;
		std::cout<< "  Index_FRAME = "<< Index_FRAME <<std::endl;
		if(flag[Index_OP]==1){
			printf("  Channel Matrix=\n");
			xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, xf::solver::NoTranspose>(
							A, "   ", print_precision, 1);
					printf("  Weights=\n");
					xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
							Weights, "   ", print_precision, 1);
					printf("  Weights_expected (var_Noise=1)=\n");
					xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
							Weights_expected, "   ", print_precision, 1);
					printf("  Weights_delta=\n");
					xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, xf::solver::NoTranspose>(
							I_delta, "   ", print_precision, 1);
					printf("  ratio= ");
					std::cout<<I_DUT_ratio<<std::endl;
					printf("  matched? ");
					std::cout<< matched_expected_Weights <<std::endl;
		}
		else{
			printf("  Bits Stream=\n");
			xf::solver::print_matrix<ROWSCOLSA, N_FRAMEperOP, MATRIX_IN_T, xf::solver::NoTranspose>(
					in_Data, "   ", print_precision, 1);
			printf("  Precoded=\n");
			xf::solver::print_matrix<ROWSCOLSA, N_FRAMEperOP, QR_INV_TYPE, xf::solver::NoTranspose>(
					PrecodedMMSE, "   ", print_precision, 1);
			Index_FRAME++;
		}
	}
	Index_OP++;
	}


// 循环读文件
// ====================================================================
//---------------- New code post-test review ----------
//    for (unsigned int imat = 0; imat < NUM_MAT_TYPES; imat++) {
//        // Test which matrix type to run
//        if ((mat_type == 0 || imat + 1 == mat_type)) {
//            for (long unsigned i = 0; i < num_tests; i++) {
//                if ((imat == 10) && i > 0) {
//                    // Skip the too large one
//                    break;
//                }

//====================================================================

//            } // End of test loop
//            printf("\n");
//        } // Type test
//    }     // Matrix type loop



if (pass_fail == 1) {
	std::cout << "TB:Fail" << std::endl;
} else {
	std::cout << "TB:Pass" << std::endl;
}
std::cout << "" << std::endl;
return (pass_fail);
}
