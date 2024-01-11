#pragma once

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <climits>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <queue>
#define NOMINMAX
#include "exprtk.hpp"
#include <iomanip>
#include <future>
#include <thread>
#include <cstring>
#include <string>
#include <mutex>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>
#include <omp.h>
#include <mpi.h>
#include <fftw3.h>
//#define BLIS_DISABLE_BLAS_DEFS
#include <cblas.h>
#include <lapacke.h>

//#define LAPACKE_dlamch dlamch_
//#define LAPACKE_zhpevx zhpevx_
//#define LAPACK_COL_MAJOR col_major_

//#define MKL_Complex16 std::complex<double>

//#include "mkl.h"
//#include "mkl_vsl.h"
//#include "mkl_dfti.h"

#ifdef USE_RESTRICT
#else
#define __restrict
#endif

typedef enum { UpdateSent, RequestSent, JobSent, Complete, AmInitializing, AmEigenSolving, AmSimulating, AmIdle, AmDone, AmRoot } MPITag;
inline int MPI_Root_Proc = 0;

#include "ThreadPool.h"
#include "ProgressTracker.h"
#include "LinuxFuncs.h"
#include "PhysCon.h"
#include "AbsorptiveRegions.h"
#include "MathTools.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include "Measurers.h"
#include "Potentials.h"
#include "SimulationManager.h"
#include "MultiSimulationManager.h"
#include "HHGFunctions.h"
#include "ConfigParser.h"
#include "ParallelSimParser.h"
#include "ThreadParser.h"