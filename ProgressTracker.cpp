#include "stdafx.h"


std::string toHMS(double sec){
	int hrs = (int)(sec/3600.0);
	int mns = (int)((sec - 3600*hrs)/60.0);
	int scs = (int)(sec - (3600*hrs + 60*mns));
	boost::format fmt = boost::format("%04d:%02d:%02d") % hrs % mns % scs;
	return fmt.str();
}

template <typename T0, typename T1, typename T2, typename T3, typename T4>
void printRow(std::ostream& os, T0 const& job, T1 const& proc, T2 const& prog,
                T3 const& timeSpent, T4 const& timeRemaining){
	os << std::setw(5) << job << std::setw(5) << proc << std::setw(30) << prog
       << std::setw(12) << timeSpent << std::setw(12) << timeRemaining << "\n";
}

std::string progressBar(int prog){
	std::stringstream str = std::stringstream();
	if (prog != -1){
		str << "[";
		for(int i = 0; i < 20; i++)
			if((i+1)*5 <= prog)
				str << "=";
			else if(i*5 <= prog && prog != 0)
				str << ">";
			else
				str << " ";
		str << "] " << (boost::format("%3d") % prog).str();
	}
	else{
		str << "Progress    \%";
	}
	return str.str();
}


ProgressTrackerShared::ProgressTrackerShared()
{
}

ProgressTrackerShared::~ProgressTrackerShared()
{
}

void ProgressTrackerShared::update(int idx, double prog, int tRemain) {
	int match = 0;
	uint i = 0;
	mtx.lock();
	for (; i < idxs.size() && !match; i++)
		match = idxs[i] == idx;
	i--;
	if (!match) {
		idxs.push_back(idx);
		progs.push_back(prog);
		ts.push_back(tRemain);
	}
	else {
		progs[i] = prog;
		ts[i] = tRemain;
	}

	std::time_t result = std::time(nullptr);
	std::fstream fil;

	try {
		fil = std::fstream("prog.txt", std::ios::out | std::ios::trunc);
	}
	catch (std::exception&) {
		std::cout << "Failed to open prog.txt" << std::endl;
		return;
	}

	fil << std::asctime(std::localtime(&result)) << "\n";

	for (i = 0; i < idxs.size(); i++) {
		if (progs[i] < 1)
			fil << std::setw(3) << idxs[i] << ": " << std::setw(3) << (int)(progs[i] * 100.0) << "\% | ETA: " << std::setw(7) << ts[i] << "s\n";
		else
			fil << std::setw(3) << idxs[i] << ": " << "DONE\n";
	}
	
	
	try {
		fil.close();
	}
	catch (std::exception&) {
		std::cout << "Failed to close prog.txt" << std::endl;
		return;
	}
	mtx.unlock();
	
}


ProgressTrackerMPI::ProgressTrackerMPI(int nJobs)
	: nJobs(nJobs)
{
	jobProg = new int[nJobs];
	std::fill_n(jobProg, nJobs, 0);

	jobProc = new int[nJobs];
	std::fill_n(jobProc, nJobs, -1);

	jobStart = new std::chrono::time_point<std::chrono::system_clock>[nJobs];

	secondsSpent = new double[nJobs];
	std::fill_n(secondsSpent, nJobs, 0.0);

	secondsRemaining = new double[nJobs];
	std::fill_n(secondsRemaining, nJobs, 0.0);
}

ProgressTrackerMPI::~ProgressTrackerMPI(){

}

//update job progress (prog = 0-100, 100 interpreted as complete)
void ProgressTrackerMPI::update(int job, int prog){
	jobProg[job] = prog;

	latestTime = std::chrono::system_clock::now();

	if(prog == 0)
		jobStart[job] = std::chrono::system_clock::now();
	else{
		dur = latestTime - jobStart[job];
		secondsSpent[job] = dur.count();
		secondsRemaining[job] = (secondsSpent[job]/(double)prog)*(100-prog);
	}
}

void ProgressTrackerMPI::jobAssigned(int job, int proc){
	jobProc[job] = proc;
}

void ProgressTrackerMPI::output(){

	std::fstream fil;

	try {
		fil = std::fstream("prog.txt", std::ios::out | std::ios::trunc);
	}
	catch (std::exception&) {
		std::cout << "Failed to open prog.txt" << std::endl;
		return;
	}

	//system time
	std::time_t result = std::time(nullptr);
	fil << std::asctime(std::localtime(&result)) << "\n";

	printRow(fil, "Job", "Proc", progressBar(-1), "TSpent", "TRemain");
	for (int i = 0; i < nJobs; i++) {
		if (jobProc[i] == - 1)
			printRow(fil, i, "UA", progressBar(0), toHMS(0), "N/A");
		else if (jobProg[i] == 0)
			printRow(fil, i, jobProc[i], progressBar(jobProg[i]), toHMS(0), "N/A");
		else
			printRow(fil, i, jobProc[i], progressBar(jobProg[i]), toHMS(secondsSpent[i]), toHMS(secondsRemaining[i]));
	}

	try {
		fil.close();
	}
	catch (std::exception&) {
		std::cout << "Failed to close prog.txt" << std::endl;
		return;
	}

}