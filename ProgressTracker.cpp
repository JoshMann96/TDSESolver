#include "ProgressTracker.h"


std::string toHMS(double sec){
	int hrs = (int)(sec/3600.0);
	int mns = (int)((sec - 3600*hrs)/60.0);
	int scs = (int)(sec - (3600*hrs + 60*mns));
	boost::format fmt = boost::format("%04d:%02d:%02d") % hrs % mns % scs;
	return fmt.str();
}

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void printRow(std::ostream& os, T0 const& job, T1 const& proc, T2 const& prog,
                T3 const& timeSpent, T6 const& timeStatSpent, T4 const& timeRemaining, T5 const& currentOperation){
	os << std::setw(5) << job << std::setw(5) << proc << std::setw(10) << currentOperation << std::setw(30) << prog
       << std::setw(12) << timeSpent << std::setw(12) << timeStatSpent << std::setw(12) << timeRemaining << "\n";
}

std::string progressBar(int prog){
	std::stringstream str = std::stringstream();
	if (prog >= 0){
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
	else if(prog == -2)
		str << "Progress    \%";
	else if(prog == -1)
		str << "";
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
	std::fill_n(jobProg, nJobs, -1);

	jobProc = new int[nJobs];
	std::fill_n(jobProc, nJobs, -1);

	jobStart = new std::chrono::time_point<std::chrono::system_clock>[nJobs];
	std::fill_n(jobStart, nJobs, std::chrono::system_clock::now());

	jobStatStart = new std::chrono::time_point<std::chrono::system_clock>[nJobs];
	std::fill_n(jobStatStart, nJobs, std::chrono::system_clock::now());

	jobDone = new std::chrono::time_point<std::chrono::system_clock>[nJobs];
	std::fill_n(jobStatStart, nJobs, std::chrono::system_clock::now());

	secondsRemaining = new double[nJobs];
	std::fill_n(secondsRemaining, nJobs, 0.0);

	currentOperation = new MPITag[nJobs];
	std::fill_n(currentOperation, nJobs, MPITag::AmIdle);
}

ProgressTrackerMPI::~ProgressTrackerMPI(){

}

//update job progress (prog = 0-100, 100 interpreted as complete)
void ProgressTrackerMPI::update(int job, int prog){
	jobProg[job] = prog;
	latestTime = std::chrono::system_clock::now();

	if(prog == 0)
		jobStatStart[job] = latestTime;
	else{
		dur = latestTime - jobStatStart[job];
		secondsRemaining[job] = (dur.count()/(double)prog)*(100-prog);
	}
}

void ProgressTrackerMPI::updateStatus(int job, MPITag stat){
	latestTime = std::chrono::system_clock::now();
	dur = latestTime - jobStart[job];
	std::printf("Job %3d spent %10s in %s", job, toHMS(dur.count()).c_str(), tagToString.at(currentOperation[job]));
	std::cout << std::endl;

	currentOperation[job] = stat;
	jobStatStart[job] = latestTime;
	if(stat == MPITag::AmDone)
		jobDone[job] = latestTime;
}

void ProgressTrackerMPI::jobAssigned(int job, int proc){
	jobProc[job] = proc;
	jobProg[job] = -1;
	jobStart[job]= std::chrono::system_clock::now();
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

	latestTime = std::chrono::system_clock::now();

	printRow(fil, "Job", "Proc", progressBar(-2), "Job T", "Status T", "ETA", "Status");
	for (int i = 0; i < nJobs; i++) {
		if (currentOperation[i] == MPITag::AmDone)
			dur = jobDone[i] - jobStart[i];
		else
			dur = latestTime - jobStart[i];
		std::chrono::duration<double> statDur = latestTime - jobStatStart[i];
		if (jobProc[i] == - 1)
			printRow(fil, i, "UA", progressBar(-1),  "N/A",  "N/A", "N/A", "Queue");
		else if (jobProg[i] == 0)
			printRow(fil, i, jobProc[i], progressBar(jobProg[i]), toHMS(dur.count()), toHMS(statDur.count()), "N/A", tagToString.at(currentOperation[i]));
		else
			printRow(fil, i, jobProc[i], progressBar(jobProg[i]), toHMS(dur.count()), toHMS(statDur.count()), toHMS(secondsRemaining[i]), tagToString.at(currentOperation[i]));
	}

	try {
		fil.close();
	}
	catch (std::exception&) {
		std::cout << "Failed to close prog.txt" << std::endl;
		return;
	}

}