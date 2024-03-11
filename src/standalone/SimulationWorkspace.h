#pragma once
#include "STLNCommonHeader.h"
#include "ParsingTools.h"
#include "LinuxFuncs.h"
#include "ProgressTracker.h"
#include "ThreadParser.h"


//creates variable/array workspace from simulation input file during definition stage
//then dispatches MPI processes to execute rest of input using each unique combination of array elements
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
	//read input file
	void executeScript();
	//read one line of input (definition stage)
	int readCommand(std::string input);
	//creates new variable in workspace
	void processNewVariable(std::string input);
	//creates new array in workspace
	void processNewArray(std::string input);

	//root proc: manage and distribute simulation assignments
	//calc proc: receive assignment and execute remainder of input file using assigned elements of array variables
	int distributeAndCompute();
};
 
