/*
 * Copyright 2021 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _KERNEL_QR_INV_HPP_
#define _KERNEL_QR_INV_HPP_

#include "hls_stream.h"

#include "dut_type.hpp"

//const int ROWSCOLSA = QR_INV_ROWSCOLS;
const int ROWSCOLSA = 8;
const int N_FRAMEperOP = 1;
const int OP_LENGTH = 11+1;
const int FRAME_LENGTH = 10+1;


extern "C" int func_matrixMultiply(hls::stream<MATRIX_IN_T>& A,
                                   hls::stream<MATRIX_OUT_T>& B);
extern "C" int func_qr_inverse(hls::stream<MATRIX_IN_T>& matrixAStrm,
                                   hls::stream<MATRIX_OUT_T>& matrixInverseAStrm);
//void
//func_matrixAddition(hls::stream<MATRIX_IN_T>& matrixAStrm,
//					hls::stream<MATRIX_IN_T>& matrixBStrm,
//		            hls::stream<MATRIX_OUT_T>& matrixAddABStrm);

void
func_matrixAddition2(hls::stream<MATRIX_IN_T>& matrixAStrm,
		            hls::stream<MATRIX_OUT_T>& matrixAddABStrm,
					float var_Noise);


//extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& matrixAStrm,
//								   hls::stream<MATRIX_OUT_T>& matrixMMSEH,
//								   float var_Noise);

void my_Precoder(hls::stream<MATRIX_T>& WeightsStrm,
									hls::stream<MATRIX_IN_T>& InputBitStrm,
								   hls::stream<MATRIX_OUT_T>& precodedMatrixStrm);

extern "C" int compute_Weights(hls::stream<MATRIX_IN_T>& matrixAStrm,
								   hls::stream<MATRIX_T>& WeightsStrm,
								   float var_Noise);


extern "C" int kernel_qr_inverse_0(hls::stream<MATRIX_IN_T>& in_channelMatrix_Strm,
									hls::stream<MATRIX_IN_T>& in_BitStrm,
									hls::stream<MATRIX_T>& in_weights_Strm,
									hls::stream<MATRIX_T>& out_weightsUpdate_Strm,
								   hls::stream<MATRIX_OUT_T>& out_precodedMatrixStrm,
								   float var_Noise,
								   int flag) ;


#endif
