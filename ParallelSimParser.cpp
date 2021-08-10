#include "stdafx.h"



ParallelSimParser::ParallelSimParser(std::fstream *fil){
	ParallelSimParser::fil = fil;
}


ParallelSimParser::~ParallelSimParser()
{
}

void ParallelSimParser::readConfig() {
	do { std::getline(*fil, curLine); } while (readCommand(curLine) && !(fil->eof()));
}

int ParallelSimParser::readCommand(std::string input) {
	parsingTools::split(input, flds, ' ');
	if (input[0] == '#')
		std::cout << input.substr(1, input.length() - 1) << std::endl;
	else if (std::strstr(flds->at(0).c_str(), "MAX_SIMS"))
		maxSims = std::stoi(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "MAX_ELEC"))
		maxElecThreads = std::stoi(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "DEF_ARR"))
		processNewArray(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "FINISH_DEF")) {
		std::cout << "Reached end of global initialization. Now branching off...\n" << std::endl;
		return branchOff();
	}
	else if (std::strstr(flds->at(0).c_str(), "DEF"))
		processNewVariable(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "MASTER_FOL"))
		OSSpecificFuncs::createFolder(flds->at(1).c_str());
	else {
		std::cout << "- " << input << std::endl;
	}
	return 1;
}

void ParallelSimParser::processNewVariable(std::string input) {
	parsingTools::split(input, flds, '=');
	varNames.push_back(flds->at(0));
	parsingTools::split(flds->at(1), flds, '*');
	double val = 1;
	int x;
	for (int i = 0; i < flds->size(); i++) {
		x = parsingTools::findMatch(&(flds->at(i)), &varNames);
		if (x == -1)
			val *= std::stod(flds->at(i));
		else
			val *= var[x];
	}
	var.push_back(val);
	std::cout << "Created new variable: " << varNames.back() << " = " << var.back() << std::endl;
}

void ParallelSimParser::processNewArray(std::string input) {
	parsingTools::split(input, flds, '=');
	arrVarNames.push_back(flds->at(0));
	parsingTools::split(flds->at(1), flds, '|');
	std::string end;
	if (flds->size() > 1)
		end = flds->at(1);
	parsingTools::split(flds->at(0), flds, ',');
	arrVar.push_back(new double[flds->size()]);
	arrVarSizes.push_back(flds->size());
	int x;
	for (int i = 0; i < flds->size(); i++) {
		x = parsingTools::findMatch(&(flds->at(i)), &varNames);
		if (x == -1)
			arrVar.back()[i] = std::stod(flds->at(i));
		else
			arrVar.back()[i] = var[x];
	}
	if (end.size()>0) {
		parsingTools::split(end, flds, '*');
		double mulVal = 0;
		if (flds->size() > 1) {
			for (int j = 1; j < flds->size(); j++) {
				x = parsingTools::findMatch(&(flds->at(j)), &varNames);
				if (x == -1)
					mulVal = std::stod(flds->at(j));
				else
					mulVal = var[x];
				for (int i = 0; i < arrVarSizes.back(); i++)
					arrVar.back()[i] *= mulVal;
			}
		}
	}
	std::cout << "Created new array: " << arrVarNames.back() << " = [";
	for (int i = 0; i < arrVarSizes.back() - 1; i++)
		std::cout << arrVar.back()[i] << ", ";
	std::cout << arrVar.back()[arrVarSizes.back() - 1] << "]" << std::endl;
}

int ParallelSimParser::branchOff() {
	if (arrVar.size() == 0) {
		std::cout << "Expected at least one array. Terminating simulations." << std::endl;
		return 1;
	}
	else {
		int numSims = 1;
		for (int i = 0; i < arrVarSizes.size(); i++) {
			numSims *= arrVarSizes[i];
		}
		std::cout << "Number of simulations to run: " << numSims << ", " << maxSims << " at a time." << std::endl;

		std::stringstream filContents = std::stringstream();
		filContents << fil->rdbuf();
		fil->close();

		ProgressTracker * prg = new ProgressTracker();

		for (int i = 0; i < arrVar.size(); i++) {
			varNames.push_back(arrVarNames.at(i));
		}

		ThreadPool * pool = new ThreadPool(maxSims);
		std::vector<std::future<int>> res;

		std::mutex mtx;

		for (int i = 0; i < numSims; i++) {
			res.emplace_back(
				pool->enqueue([this, i, &filContents, prg, &mtx] {
				omp_set_num_threads(maxElecThreads);
				std::vector<std::string> myVarNames = std::vector<std::string>(varNames);
				std::vector<double> myVar = std::vector<double>(var);
				int curVarGen = i;
				for (int j = 0; j < arrVarSizes.size(); j++) {
					int idx = curVarGen % arrVarSizes.at(j);
					myVar.push_back(arrVar.at(j)[idx]);
					curVarGen /= arrVarSizes.at(j);
				}
				std::stringstream* myFil = new std::stringstream(filContents.str());

				ThreadParser* mySim = new ThreadParser(myFil, prg, myVarNames, myVar, i, &mtx);
				mySim->readConfig();
				delete mySim;
				return 0;
			}));
		}

		for (auto&& result : res)
			result.get();
		res.clear();

		delete pool;
	}
	return 0;
}