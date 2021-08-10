#include "stdafx.h"



ProgressTracker::ProgressTracker()
{
}


ProgressTracker::~ProgressTracker()
{
}

void ProgressTracker::update(int idx, double prog, int tRemain) {
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
