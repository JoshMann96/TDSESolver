#include "SimulationWorkspace.h"



SimulationWorkspace::SimulationWorkspace(std::fstream *fil) : fil(fil){}

SimulationWorkspace::~SimulationWorkspace(){}

void SimulationWorkspace::executeScript() {
	//read and execute each line in the file so long as readCommand 
	do { std::getline(*fil, curLine); } while (readCommand(curLine) && !(fil->eof()));
}

int SimulationWorkspace::readCommand(std::string input) {
	parsingTools::split(input, flds, ' ');
	if (input[0] == '#') //comment
		return 1;
	else if (std::strstr(flds->at(0).c_str(), "DEF_ARR")) //defining array variable
		processNewArray(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "FINISH_DEF")) //finished with definitions, computation processes create thread parsers to evaluate the rest of the cfg
		return distributeAndCompute();
	else if (std::strstr(flds->at(0).c_str(), "DEF")) //defining scalar variable
		processNewVariable(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "MASTER_FOL")) //available but obsolete
		OSSpecificFuncs::createFolder(flds->at(1).c_str());
	else {
		std::cout << "- " << input << std::endl;
		return 0;
	}
	return 1;
}

void SimulationWorkspace::processNewVariable(std::string input) {
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
	//std::cout << "Created new variable: " << varNames.back() << " = " << var.back() << std::endl;
}

void SimulationWorkspace::processNewArray(std::string input) {
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
	//std::cout << "Created new array: " << arrVarNames.back() << " = [";
	//for (int i = 0; i < arrVarSizes.back() - 1; i++)
		//std::cout << arrVar.back()[i] << ", ";
	//std::cout << arrVar.back()[arrVarSizes.back() - 1] << "]" << std::endl;
}

int SimulationWorkspace::distributeAndCompute() {
	MPI_Barrier(MPI_COMM_WORLD);
	if (arrVar.size() == 0) {
		std::cout << "Expected at least one array. Terminating simulations." << std::endl;
		return 1;
	}
	else {
		//calc num of sims
		int numSims = 1;
		for (int i = 0; i < arrVarSizes.size(); i++) {
			numSims *= arrVarSizes[i];
		}

		//load file contents
		std::stringstream filContents = std::stringstream();
		filContents << fil->rdbuf();
		fil->close();

		for (int i = 0; i < arrVar.size(); i++) {
			varNames.push_back(arrVarNames.at(i));
		}

		int size, rank;

    	MPI_Comm_size(MPI_COMM_WORLD, &size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		MPI_Status stat;
		int buf;

		if (rank == MPI_Root_Proc){
			int curJob = 0;
			int nProcDone = 0;
			int* assignedJobs = new int[size];
			std::fill_n(assignedJobs, size, MPI_Root_Proc);

			std::chrono::high_resolution_clock::time_point* jt = new std::chrono::high_resolution_clock::time_point[size];
			std::fill_n(jt, size, std::chrono::high_resolution_clock::now());

			ProgressTrackerMPI * prg = new ProgressTrackerMPI(numSims);

			while(nProcDone < size - 1){
				MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
				switch(stat.MPI_TAG){
					case MPITag::UpdateSent :
						prg->update(assignedJobs[stat.MPI_SOURCE], buf);
						break;
					case MPITag::RequestSent :
						//tell progress tracker that job has been complete
						if(assignedJobs[stat.MPI_SOURCE] != MPI_Root_Proc)
							prg->update(assignedJobs[stat.MPI_SOURCE], 100);
						//assign job if available, inform tracker of job assignment
						if(curJob < numSims){
							std::cout << "Sending Job " << std::setw(4) << curJob << " to Proc " << std::setw(4) << stat.MPI_SOURCE << std::endl;
							MPI_Ssend(&curJob, 1, MPI_INT, stat.MPI_SOURCE, MPITag::JobSent, MPI_COMM_WORLD);
							assignedJobs[stat.MPI_SOURCE] = curJob;
							prg->jobAssigned(curJob, stat.MPI_SOURCE);
                        	curJob++;
						}
						else{
							std::cout << "Proc " << std::setw(4) << stat.MPI_SOURCE << " Fired" << std::endl;
							MPI_Ssend(nullptr, 0, MPI_INT, stat.MPI_SOURCE, MPITag::Complete, MPI_COMM_WORLD);
							nProcDone++;
						}
						break;
					default:
						prg->updateStatus(assignedJobs[stat.MPI_SOURCE], (MPITag)stat.MPI_TAG);
						break;
				}
				prg->output();
			}
		}
		else{
			//omp_set_num_threads(omp_get_max_threads());

			//set wisdom file according to number of threads used
			char* wisdomFile = new char[50];
			std::snprintf(wisdomFile, 50, "fftw_nt_%04d.wisdom", omp_get_max_threads());

			int done = 0, job, curVarGen;
			std::vector<std::string> myVarNames;
			std::vector<double> myVar;
			std::stringstream* myFil;
			ThreadParser* mySim;
			while(!done){
				MPI_Ssend(nullptr, 0, MPI_INT, MPI_Root_Proc, MPITag::RequestSent, MPI_COMM_WORLD); //request data
				MPI_Recv(&job, 1, MPI_INT, MPI_Root_Proc, MPI_ANY_TAG, MPI_COMM_WORLD, &stat); //listen for data
				switch(stat.MPI_TAG){
					case MPITag::JobSent : //job available, run calculation
						MPI_Ssend(nullptr, 0, MPI_INT, MPI_Root_Proc, MPITag::AmInitializing, MPI_COMM_WORLD);
						myVarNames = std::vector<std::string>(varNames);
						myVar = std::vector<double>(var);
						curVarGen = job;
						for (int j = 0; j < arrVarSizes.size(); j++) {
							int idx = curVarGen % arrVarSizes.at(j);
							myVar.push_back(arrVar.at(j)[idx]);
							curVarGen /= arrVarSizes.at(j);
						}
						myFil = new std::stringstream(filContents.str());

						mySim = new ThreadParser(myFil, myVarNames, myVar, job);

						//try to import wisdom specifically for number of available threads
						fftw_init_threads();
						fftw_import_wisdom_from_filename(wisdomFile);

						mySim->readConfig();

						fftw_export_wisdom_to_filename(wisdomFile);

						MPI_Ssend(nullptr, 0, MPI_INT, MPI_Root_Proc, MPITag::AmDone, MPI_COMM_WORLD);
						delete mySim;
						break;
					case MPITag::Complete :
						done = 1;
						break;
				}
			}
		}
	}
	return 0;
}