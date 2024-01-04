#pragma once

std::string toHMS(double sec);

template <typename T0, typename T1, typename T2, typename T3, typename T4>
void printRow(std::ostream& os, T0 const& job, T1 const& proc, T2 const& prog,
                T3 const& timeSpent, T4 const& timeRemaining);

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
	std::chrono::time_point<std::chrono::system_clock>* jobStart, latestTime;
	double * secondsSpent, * secondsRemaining;

	std::chrono::duration<double> dur;

public:
	ProgressTrackerMPI(int nJobs);
	~ProgressTrackerMPI();
	void update(int job, int prog);
	void output();
	void jobAssigned(int job, int proc);
};