#pragma once
class ParallelSimParser
{
private:
	std::vector<double> var = std::vector<double>(0);
	std::vector<double*> arrVar = std::vector<double*>(0);
	std::vector<std::string> arrVarNames = std::vector<std::string>(0), varNames = std::vector<std::string>(0);
	std::vector<int> arrVarSizes = std::vector<int>(0);
	std::fstream * fil;
	int maxSims = 1;
	int maxElecThreads = 1;
	std::vector<std::string> * flds = new std::vector<std::string>(0);
	std::string curLine = "";
public:
	ParallelSimParser(std::fstream *fil);
	~ParallelSimParser();
	void readConfig();
	int readCommand(std::string input);
	void processNewVariable(std::string input);
	void processNewArray(std::string input);
	int branchOff();
};

