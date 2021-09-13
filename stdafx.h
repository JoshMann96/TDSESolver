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
#include "matplotlibcpp.h"
#include <iomanip>
#include <future>
#include <thread>
#include <cstring>
#include <mutex>
#include <boost/filesystem.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_vsl.h"
#include "mkl_dfti.h"
#include "omp.h"

#ifdef USE_RESTRICT
#else
#define __restrict
#endif

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
#include "TDSEIterator1D.h"
#include "SimulationManager.h"
//#include "SingleSimulationManager.h"
#include "MultiSimulationManager.h"
#include "HHGFunctions.h"
#include "ConfigParser.h"
#include "ParallelSimParser.h"
#include "ThreadParser.h"