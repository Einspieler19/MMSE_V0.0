#
# Copyright 2019-2021 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

source settings.tcl

set PROJ "qr_inverse_test.prj"
set SOLN "sol1"

if {![info exists CLKP]} {
  set CLKP 300MHz
}

open_project -reset $PROJ

add_files "${XF_PROJ_ROOT}/kernel/kernel_qr_inverse_0.cpp" -cflags "-DQR_INV_ROWSCOLS=3 -DSEL_ARCH=0 -D_DATA_PATH=${XF_PROJ_ROOT}/datas/ -I./ -I${XF_PROJ_ROOT}/host/ -I${XF_PROJ_ROOT}/kernel/ -I${XF_PROJ_ROOT}/L1/tests/ -I${XF_PROJ_ROOT}/include/ -I${XF_PROJ_ROOT}/include/hw -I${XF_PROJ_ROOT}/include2 -I${XF_PROJ_ROOT}include3/"
add_files -tb "${XF_PROJ_ROOT}/host/test_qr_inverse.cpp" -cflags "-DQR_INV_ROWSCOLS=3 -DSEL_ARCH=0 -D_DATA_PATH=${XF_PROJ_ROOT}/datas/ -I./ -I${XF_PROJ_ROOT}/host/ -I${XF_PROJ_ROOT}/kernel/ -I${XF_PROJ_ROOT}/tests/ -I${XF_PROJ_ROOT}/include/ -I${XF_PROJ_ROOT}include/hw -I ./host -I${XF_PROJ_ROOT}/include3/"
set_top kernel_qr_inverse_0

open_solution -reset $SOLN




set_part $XPART
create_clock -period $CLKP

if {$CSIM == 1} {
  csim_design
}

if {$CSYNTH == 1} {
  csynth_design
}

if {$COSIM == 1} {
  cosim_design
}

if {$VIVADO_SYN == 1} {
  export_design -flow syn -rtl verilog
}

if {$VIVADO_IMPL == 1} {
  export_design -flow impl -rtl verilog
}

exit
