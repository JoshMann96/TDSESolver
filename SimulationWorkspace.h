#pragma once
#include "CommonHeader.h"
#include "ParsingTools.h"
#include "LinuxFuncs.h"
#include "ProgressTracker.h"
#include "ThreadParser.h"


//reads config file and listens to root process for receiving assignments
class SimulationWorkspace
{
private:
	std::vector<double> var = std::vector<double>(0);
	std::vector<double*> arrVar = std::vector<double*>(0);
	std::vector<std::string> arrVarNames = std::vector<std::string>(0), varNames = std::vector<std::string>(0);
	std::vector<int> arrVarSizes = std::vector<int>(0);
	std::fstream * fil;
	std::vector<std::string> * flds = new std::vector<std::string>(0);
	std::string curLine = "";
public:
	SimulationWorkspace(std::fstream *fil);
	~SimulationWorkspace();
	void executeScript();
	int readCommand(std::string input);
	void processNewVariable(std::string input);
	void processNewArray(std::string input);

	int distributeAndCompute();
};
 
