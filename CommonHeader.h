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
#include <thread>
#include <cstring>
#include <string>

#include <mpi.h>
#include <omp.h>

#ifdef USE_RESTRICT
#else
#define __restrict
#endif

typedef enum { UpdateSent, RequestSent, JobSent, Complete, AmInitializing, AmEigenSolving, AmSimulating, AmIdle, AmDone, AmRoot } MPITag;
inline int MPI_Root_Proc = 0;

#include "MathTools.h"
#include "PhysCon.h"