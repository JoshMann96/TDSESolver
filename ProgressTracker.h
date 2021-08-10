#pragma once
class ProgressTracker
{
private:
	std::vector<double> progs;
	std::vector<int> ts;
	std::vector<int> idxs;
	std::mutex mtx;
public:
	ProgressTracker();
	~ProgressTracker();
	void update(int idx, double prog, int tRemain);
};