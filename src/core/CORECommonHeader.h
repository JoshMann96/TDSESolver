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
#include <stdexcept>

#include <omp.h>

#ifdef USE_RESTRICT
#else
#define __restrict
#endif


#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/tools/rational.hpp>
#include "MathTools.h"
#include "PhysCon.h"