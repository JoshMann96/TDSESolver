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
#define NOMINMAX
#include <iomanip>
#include <future>
#include <cstring>
#include <string>

#include <mpi.h>
#include <omp.h>

#ifdef USE_RESTRICT
#else
#define __restrict
#endif

#define lapack_int long int

#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>

typedef enum { UpdateSent, RequestSent, JobSent, Complete, AmInitializing, AmEigenSolving, AmSimulating, AmIdle, AmDone, AmRoot } MPITag;
const boost::unordered_map<MPITag,const char*> tagToString = boost::assign::map_list_of
		(MPITag::AmIdle, "Idle")
		(MPITag::AmEigenSolving, "Eigen")
		(MPITag::AmInitializing, "Init")
		(MPITag::AmSimulating, "Simu")
		(MPITag::AmDone, "Done");
inline int MPI_Root_Proc = 0;

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/tools/rational.hpp>
#include "MathTools.h"
#include "PhysCon.h"