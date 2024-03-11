#pragma once
#include "STLNCommonHeader.h"
#include "boost/format.hpp"


std::string toHMS(double sec);

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void printRow(std::ostream& os, T0 const& job, T1 const& proc, T2 const& prog,
                T3 const& timeSpent, T6 const& timeStatSpent, T4 const& timeRemaining, T5 const& currentOperation);

std::string progressBar(int prog);

class ProgressTrackerShared
{
private:
	std::vector<double> progs;
	std::vector<int> ts;
	std::vector<int> idxs;
	std::mutex mtx;
public:
	ProgressTrackerShared();
	~ProgressTrackerShared();
	void update(int idx, double prog, int tRemain);
};

class ProgressTrackerMPI
{
private:
	int nJobs;
	int * jobProg, * jobProc;
	std::chrono::time_point<std::chrono::system_clock>* jobStart, *jobStatStart, *jobDone, latestTime;
	double * secondsRemaining;

	std::chrono::duration<double> dur;

	MPITag* currentOperation;

public:
	ProgressTrackerMPI(int nJobs);
	~ProgressTrackerMPI();
	void update(int job, int prog);
	void updateStatus(int job, MPITag stat);
	void output();
	void jobAssigned(int job, int proc);
};