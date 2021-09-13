#include "stdafx.h"



ThreadParser::ThreadParser(std::stringstream *fil, ProgressTracker *prg, std::vector<std::string> varNames, std::vector<double> var, int simIdx, std::mutex *mtx)
{
	ThreadParser::fil = fil;
	ThreadParser::prg = prg;
	ThreadParser::varNames = varNames;
	ThreadParser::var = var;
	ThreadParser::simIdx = simIdx;
	ThreadParser::mtx = mtx;
}

ThreadParser::~ThreadParser() {
	delete sim;
}

void ThreadParser::readConfig() {
	do { std::getline(*fil, curLine); } while (readCommand() && !(fil->eof()));
}

int ThreadParser::addPotential(std::string input) {
	//To add a new potential function, do the instructions in the header file and return here.
	//Add the next case in the list and do what you expect it to do.
	//If you expect a string as one of the variables, it may be accessed with inputText (std::string).
	//All other parameters are in the vector p, in the same order as they are listed in the header file.
	//Once finished here, go add an example use to config_file_samples.cfg
	int potIdx = parsingTools::findMatch(&input, &potNames);
	std::vector<double> p = getPotBlockParameters(potIdx);
	switch (potIdx) {
	case 1:
		sim->addPotential(new Potentials::JelliumPotential(n, x, p[0], p[1], p[2], vtls::findValue(n, x, p[3])));
		break;
	case 0:
		sim->addPotential(new Potentials::JelliumPotentialBacked(n, x, p[0], p[1], p[2], p[3], p[4], vtls::findValue(n, x, p[5])));
		break;
	case 2:
		sim->addPotential(new Potentials::ShieldedAtomicPotential(n, x, p[0], p[1], p[2], p[3]));
		break;
	case 10:
		sim->addPotential(new Potentials::ElectricFieldProfileToPotential(
			n, new ElectricFieldProfiles::CylindricalToLinearProfile(n, x, p[0], p[1], p[2], p[3], p[4]),
			sim->getDX(), p[8], p[7], p[5], new Envelopes::SmoothedInitialGaussianEnvelope(p[6], p[7], p[9]), vtls::findValue(n, x, p[10])));
		break;
	case 4:
		sim->addPotential(
			new Potentials::ElectricFieldProfileToPotential(
				n, new ElectricFieldProfiles::InMetalFieldProfile(
					n, x, p[0], p[1], p[2], p[3], p[8] + PhysCon::im*p[9], p[10]),
				sim->getDX(), p[5], p[6], p[3], new Envelopes::SmoothedInitialGaussianEnvelope(
					p[4], p[6], p[7]),
				vtls::findValue(n, x, p[11])));
		break;
	case 5:
		sim->addPotential(
			new Potentials::ElectricFieldProfileToPotential(
				n, new ElectricFieldProfiles::ConstantFieldProfile(
					n, x, p[2], p[0], p[1]),
				sim->getDX(), p[6], p[5], p[3], new Envelopes::SmoothedInitialGaussianEnvelope(
					p[4], p[5], p[7]),
				vtls::findValue(n, x, p[8])));
		break;
	case 6:
		sim->addPotential(new Potentials::FilePotential(n, x, p[0], inputText.c_str(), vtls::findValue(n, x, p[2])));
		break;
	case 7:
	{
		/*Potentials::WaveFunctionSelfPotential * pot = new Potentials::WaveFunctionSelfPotential(n, (x[n - 1] - x[0]) / n, p[0], p[1], vtls::findValue(n, x, p[2]));
		sim->addPotential(pot);
		spc->push_back(pot);*/
		std::cout << "WaveFunctionSelfPotential is now a template, not a real potential." << std::endl;
		break;
	}
	case 8:
		sim->addPotential(new Potentials::BiasFieldPotential(n, x, p[5], p[6], p[0], p[1], p[2], p[3], p[4], vtls::findValue(n, x, p[7])));
		break;
	case 9:
		sim->addPotential(
			new Potentials::ElectricFieldProfileToPotential(
				n, new ElectricFieldProfiles::FileFieldProfile(
					n, x, p[0], p[3], p[2], p[4], p[5], inputText.c_str()),
				sim->getDX(), p[9], p[8], p[6], new Envelopes::SmoothedInitialGaussianEnvelope(
					p[7], p[8], p[10]),
				vtls::findValue(n, x, p[11])));
		break;
	case 3:
		sim->addPotential(new Potentials::ElectricFieldProfileToPotential(
			n, new ElectricFieldProfiles::CylindricalToCutoffProfile(n, x, p[0], p[1], p[2], p[3], p[4], p[11]),
			sim->getDX(), p[8], p[7], p[5], new Envelopes::SmoothedInitialGaussianEnvelope(p[6], p[7], p[9]), vtls::findValue(n, x, p[10])));
		break;
	case 11:
		sim->addPotential(new Potentials::FiniteBox(n, x, p[0], p[1], p[2], vtls::findValue(n, x, p[3])));
		break;
	case 13:
	{
		Potentials::SurfaceSpaceCharge* pot = new Potentials::SurfaceSpaceCharge(n, sim->getDX(), p[0], nelec, sim->getPotPointer(), wght, dens, vtls::findValue(n, x, p[1]), vtls::findValue(n, x, p[2]), vtls::findValue(n, x, p[3]));
		sim->addPotential(pot);
		spc->push_back(pot);
		break;
	}
	case 12:
	{
		Potentials::FullCylindricalSpaceCharge* pot = new Potentials::FullCylindricalSpaceCharge(n, x, sim->getDX(), p[0], p[1], nelec, sim->getPotPointer(), wght, dens, vtls::findValue(n, x, p[3]), vtls::findValue(n, x, p[4]), vtls::findValue(n, x, p[2]), vtls::findValue(n, x, p[5]));
		sim->addPotential(pot);
		spc->push_back(pot);
		break;
	}
	case 14:
	{
		Potentials::LinearBulkCylindricalFieldSpaceCharge* pot = new Potentials::LinearBulkCylindricalFieldSpaceCharge(n, x, sim->getDX(), p[0], p[1], nelec, sim->getPotPointer(), wght, dens, vtls::findValue(n, x, p[3]), vtls::findValue(n, x, p[4]), vtls::findValue(n, x, p[2]), vtls::findValue(n, x, p[5]));
		sim->addPotential(pot);
		spc->push_back(pot);
		break;
	}
	case 15:
	{
		Potentials::LinearBulkCylSectionFieldSpaceCharge* pot = new Potentials::LinearBulkCylSectionFieldSpaceCharge(n, x, sim->getDX(), p[0], p[1], p[3], nelec, sim->getPotPointer(), wght, dens, vtls::findValue(n, x, p[2]), vtls::findValue(n, x, p[4]));
		sim->addPotential(pot);
		spc->push_back(pot);
		break;
	}
	case 16:
	{
		Potentials::DielectricBulkCylindricalFieldSpaceCharge* pot = new Potentials::DielectricBulkCylindricalFieldSpaceCharge(n, x, sim->getDX(), sim->getDT(), p[0], p[1], p[5], p[6], nelec, sim->getPotPointer(), wght, dens, vtls::findValue(n, x, p[3]), vtls::findValue(n, x, p[4]), vtls::findValue(n, x, p[2]), vtls::findValue(n, x, p[7]));
		sim->addPotential(pot);
		spc->push_back(pot);
		break;
	}
	case 17:
	{
		Potentials::OhmicRetardingPotential* pot = new Potentials::OhmicRetardingPotential(n, sim->getDX(), p[1], p[2], nelec, sim->getPotPointer(), wght, dens, vtls::findValue(n, x, p[0]), vtls::findValue(n, x, p[3]));
		sim->addPotential(pot);
		spc->push_back(pot);
		break;
	}
	}
	return 1;
}

int ThreadParser::addAbsEdge(std::string input) {
	//To add a new absorptive region, do the instructions in the header file and return here.
	//Add the next case in the list and do what you expect it to do.
	//If you expect a string as one of the variables, it may be accessed with inputText (std::string).
	//All other parameters are in the vector p, in the same order as they are listed in the header file.
	//Once finished here, go add an example use to config_file_samples.cfg
	int absIdx = parsingTools::findMatch(&input, &absNames);
	std::vector<double> p = getAbsBlockParameters(absIdx);
	switch(absIdx) {
	case 0:
		sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, (int)((x[n - 1] - p[0] - x[0]) / (x[n - 1] - x[0])*n), n, -p[1]));
		sim->addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(n, (int)((x[n - 1] - p[0] - x[0]) / (x[n - 1] - x[0]) * n), n, p[1]));
		break;
	case 1:
		sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, (int)(p[0] / (x[n - 1] - x[0])*n), 0, -p[1]));
		sim->addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(n, (int)(p[0] / (x[n - 1] - x[0]) * n), 0, p[1]));
		break;
	}
	return 1;
}

int ThreadParser::setWghtCalc(std::string input) {
	int denIdx = parsingTools::findMatch(&input, &denNames);
	std::vector<double> p = getDenBlockParameters(denIdx);
	switch (denIdx) {
	case 0:
		wght = new WfcToRho::FermiGasDistro(p[0]);
		break;
	case 1:
		wght = new WfcToRho::FromDOS(p[1], p[0], p[2], inputText.c_str());
		break;
	}
	return 1;
}

int ThreadParser::setRhoCalc(std::string input) {
	int rhoIdx = parsingTools::findMatch(&input, &rhoNames);
	std::vector<double> p = getRhoBlockParameters(rhoIdx);
	switch (rhoIdx) {
	case 0:
		dens = new WfcToRho::DirectDensity();
		break;
	case 1:
		dens = new WfcToRho::GaussianSmoothedDensity(p[0]);
		break;
	}
	return 1;
}

int ThreadParser::addMeasurer(std::string input) {
	//To add a new measurer, do the instructions in the header file and return here.
	//Add the next case in the list and do what you expect it to do.
	//If you expect a string as one of the variables, it may be accessed with inputText (std::string).
	//All other parameters are in the vector p, in the same order as they are listed in the header file.
	//Once finished here, go add an example use to config_file_samples.cfg
	int meaIdx = parsingTools::findMatch(&input, &meaNames);
	std::vector<double> p = getMeaBlockParameters(meaIdx);
	switch (meaIdx) {
	case 0:
		sim->addMeasurer(new Measurers::BasicMeasurers(inputText.c_str(), n, sim->getDX(), sim->getDT(), x, fol));
		break;
	case 1:
		sim->addMeasurer(new Measurers::Vfunct(n, p[0], p[1], sim->getMaxT(), x, fol));
		break;
	case 2:
		sim->addMeasurer(new Measurers::VDPsi(vtls::findValue(n, x, p[2]), p[1], inputText.c_str(), fol));
		break;
	case 3:
		sim->addMeasurer(new Measurers::VDProbCurrent(n, sim->getDX(), vtls::findValue(n, x, p[2]), p[1], inputText.c_str(), fol));
		break;
	case 4:
		sim->addMeasurer(new Measurers::TotProb(n, sim->getDX(), fol));
		break;
	case 5:
		sim->addMeasurer(new Measurers::ExpectA(n, sim->getDX(), fol));
		break;
	case 6:
		sim->addMeasurer(new Measurers::ExpectP(n, sim->getDX(), fol));
		break;
	case 7:
		sim->addMeasurer(new Measurers::ExpectX(n, x, sim->getDX(), fol));
		break;
	case 8:
		sim->addMeasurer(new Measurers::ExpectE0(n, sim->getDX(), fol));
		break;
	case 9:
		sim->addMeasurer(new Measurers::Psi2t(n, p[0], p[1], sim->getMaxT(), sim->getDT(), x, fol));
		break;
	case 10:
		sim->addMeasurer(new Measurers::OrigPot(n, fol));
		break;
	case 11:
		sim->addMeasurer(new Measurers::TS(fol));
		break;
	case 12:
		sim->addMeasurer(new Measurers::XS(n, x, fol));
		break;
	case 13:
		sim->addMeasurer(new Measurers::DT(sim->getDT(), fol));
		break;
	case 14:
		sim->addMeasurer(new Measurers::DX(sim->getDX(), fol));
		break;
	case 15:
		sim->addMeasurer(new Measurers::NSteps(fol));
		break;
	case 16:
		sim->addMeasurer(new Measurers::NPts(n, fol));
		break;
	case 17:
		sim->addMeasurer(new Measurers::Header(inputText.c_str(), fol));
		break;
	case 18:
		sim->addMeasurer(new Measurers::DoubleConst(p[1], inputText.c_str(), fol));
		break;
	case 19:
		sim->addMeasurer(new Measurers::WignerQPD(n, p[0], p[1], p[2], p[3], p[4], sim->getMaxT(), p[5], x, fol));
		break;
	case 20:
		sim->addMeasurer(new Measurers::PsiT(n, p[2], p[1], inputText.c_str(), fol));
		break;
	case 21:
		sim->addMeasurer(new Measurers::PotT(n, p[2], p[1], inputText.c_str(), fol));
		break;
	case 22:
		sim->addMeasurer(new Measurers::VDPot(vtls::findValue(n, x, p[2]), p[1], inputText.c_str(), fol));
		break;
	case 23:
		sim->addMeasurer(new Measurers::ExpectE(n, sim->getDX(), fol));
		break;
	case 24:
	{
		int vdpos = vtls::findValue(n, x, p[3]);
		if (vdpos >= n - 1)
			vdpos = n - 2;
		sim->addMeasurer(new Measurers::VDPsi(vdpos, p[1], inputText.c_str(), fol));
		sim->addMeasurer(new Measurers::VDPsi(vdpos + 1, p[2], inputText.c_str(), fol));
		break;
	}
	case 25:
	{
		int vdpos = vtls::findValue(n, x, p[3]);
		if (vdpos >= n - 1)
			vdpos = n - 2;
		sim->addMeasurer(new Measurers::VDPot(vdpos, p[1], inputText.c_str(), fol));
		sim->addMeasurer(new Measurers::VDPot(vdpos + 1, p[2], inputText.c_str(), fol));
		break;
	}
	case 26:
		sim->addMeasurer(new Measurers::VDFluxSpec(vtls::findValue(n, x, p[2]), p[1], nelec, p[4], p[3], sim->getMaxT(), inputText.c_str(), fol));
		break;
	case 27:
		sim->addMeasurer(new Measurers::WfcRhoWeights(nelec, n, sim->getDX(), wght, fol));
		break;
	}
	return 1;
}

int ThreadParser::readCommand() {
	parsingTools::split(curLine, flds, ' ');
	if (curLine[0] == '#') {}
	else if (std::strstr(flds->at(0).c_str(), "SAVE_LOC"))
		return setSaveLoc(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "NEW_GEN_SIM"))
		return generalSimInit();
	else if (std::strstr(flds->at(0).c_str(), "NEW_HHG_BASED_SIM"))
		return hhgSimInit();
	else if (std::strstr(flds->at(0).c_str(), "DEN_CAL"))
		return setWghtCalc(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "RHO_CAL"))
		return setRhoCalc(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "ADD_POT"))
		return addPotential(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "ADD_MEA"))
		return addMeasurer(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "ADD_ABS"))
		return addAbsEdge(flds->at(1));
	else if (std::strstr(flds->at(0).c_str(), "FINISH_INIT")) {
		sim->finishInitialization();
	}
	else if (std::strstr(flds->at(0).c_str(), "WAV_DEP_IM_TIME_PROP")) {
		std::complex<double>* tpsi = new std::complex<double>[n];
		std::fill_n(tpsi, n, 1);
		sim->setPsi(tpsi);
		sim->findGroundStateWithWaveFuncDepPot(parseVal(flds->at(1)), parseVal(flds->at(2)));
	}
	else if (std::strstr(flds->at(0).c_str(), "IM_TIME_PROP")) {
		std::complex<double>* tpsi = new std::complex<double>[n];
		std::fill_n(tpsi, n, 1);
		sim->setPsi(tpsi);
		sim->findGroundState(parseVal(flds->at(1)), parseVal(flds->at(2)));
	}
	else if (std::strstr(flds->at(0).c_str(), "FIND_EIGEN_STATE"))
		sim->findEigenState(parseVal(flds->at(1)), parseVal(flds->at(2)));
	else if (std::strstr(flds->at(0).c_str(), "FIND_METAL_STATE_PSM")) {
		//mtx->lock();
		sim->findMetallicInitialState_PSM(parseVal(flds->at(1)), parseVal(flds->at(2)), 0, 0);
		sim->addMeasurer(new Measurers::ElectronNumber(sim->getNElec(), fol));
		nelec[0] = sim->getNElec();
		//mtx->unlock();
	}
	else if (std::strstr(flds->at(0).c_str(), "FIND_METAL_STATE_HOD")) {
		//mtx->lock();
		sim->findMetallicInitialState_HOD(parseVal(flds->at(1)), parseVal(flds->at(2)), 0, 0, parseVal(flds->at(3)));
		sim->addMeasurer(new Measurers::ElectronNumber(sim->getNElec(), fol));
		nelec[0] = sim->getNElec();
		//mtx->unlock();
	}
	else if (std::strstr(flds->at(0).c_str(), "FIND_METAL_STATE")) {
		sim->findMetallicInitialState(parseVal(flds->at(1)), parseVal(flds->at(2)), 0, 0);
		sim->addMeasurer(new Measurers::ElectronNumber(sim->getNElec(), fol));
		nelec[0] = sim->getNElec();
	}
	else if (std::strstr(flds->at(0).c_str(), "NEGATE_SELF_POT_INITIAL_STATE"))
		for (int i = 0; i < spc->size(); i++)
			spc->at(i)->negateGroundEffects(sim->getPsi());
	else if (std::strstr(flds->at(0).c_str(), "RUN_PARALLEL_LOG_PROGRESS"))
		sim->runParallelLogProgress(prg, simIdx);
	else if (std::strstr(flds->at(0).c_str(), "RUN_LOG_PROGRESS"))
		sim->runLogProgress(prg, simIdx);
	else if (std::strstr(flds->at(0).c_str(), "RUN_SCF_LOG_PROGRESS"))
		sim->runSCFLogProgress(prg, simIdx, parseVal(flds->at(1)));
	else if (std::strstr(flds->at(0).c_str(), "RUN_PERT_LOG_PROGRESS"))
		sim->runPertLogProgress(prg, simIdx, (int)(parseVal(flds->at(1))));
	else if (std::strstr(flds->at(0).c_str(), "RUN_OS_U2TU_LOG_PROGRESS"))
		sim->runOS_U2TULogProgress(prg, simIdx);
	else if (std::strstr(flds->at(0).c_str(), "RUN_OS_UW2TUW_LOG_PROGRESS"))
		sim->runOS_UW2TUWLogProgress(prg, simIdx);
	else if (std::strstr(flds->at(0).c_str(), "EXIT"))
		return 0;
	else if (std::strstr(flds->at(0).c_str(), "MULTI_ELEC"))
		multiElec = 1;
	return 1;
	//Add rest of potentials, measurers, absorptive edges, and finish off!
}

int ThreadParser::generalSimInit() {
	double minX = 0;
	double maxX = 0;
	double maxT = 0;
	double dx = 0;
	double dt = 0;
	do {
		std::getline(*fil, curLine);
		parsingTools::split(curLine, flds, ' ');
		if (std::strstr(flds->at(0).c_str(), "X_MIN"))
			minX = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "X_MAX"))
			maxX = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "T_MAX"))
			maxT = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "DX"))
			dx = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "DT"))
			dt = parseVal(flds->at(1));
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	n = (((int)((maxX - minX) / dx))/2)*2; // get num points, force even
	if(multiElec)
		sim = new MultiSimulationManager(n, dx, dt, maxT);
	else {
		std::cout << "SingleSimulationManager is no longer supported -- Use MultiSimulationManager with 1 simulation." << std::endl;
		throw "FEATURE REMOVED";
	}
		//sim = new SingleSimulationManager(n, dx, dt, maxT);
	x = new double[n];
	vtls::linspace(n, minX, maxX, x);
	return 1;
}

int ThreadParser::hhgSimInit() {
	double minX = 0;
	double maxX = 0;
	double maxT = 0;
	double emax = 0;
	double lam = 0;
	double err = 0;
	do {
		std::getline(*fil, curLine);
		parsingTools::split(curLine, flds, ' ');
		if (std::strstr(flds->at(0).c_str(), "X_MIN"))
			minX = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "X_MAX"))
			maxX = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "T_MAX"))
			maxT = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "E_MAX"))
			emax = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "LAM"))
			lam = parseVal(flds->at(1));
		else if (std::strstr(flds->at(0).c_str(), "ERR"))
			err = parseVal(flds->at(1));
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	double pond = HHGFunctions::getPonderomotiveEnergy(emax, lam);

	double dx = HHGFunctions::getIdealDX(pond*100.0, err);
	double dt = HHGFunctions::getIdealDT(dx);

	n = (int)((maxX - minX) / dx);
	if (multiElec)
		sim = new MultiSimulationManager(n, dx, dt, maxT);
	else {
		std::cout << "SingleSimulationManager is no longer supported -- Use MultiSimulationManager with 1 simulation." << std::endl;
		throw "FEATURE REMOVED";
	}
		//sim = new SingleSimulationManager(n, dx, dt, maxT);
	x = new double[n];
	vtls::linspace(n, minX, maxX, x);
	return 1;
}

int ThreadParser::setSaveLoc(std::string input) {
	int folLen = input.length() + 40;
	fol = new char[folLen];
	std::vector<std::string> * nflds = new std::vector<std::string>(0);
	parsingTools::split(input, nflds, '~');
	//std::stringstream tempstr;
	int fi = 0;
	if (nflds->size() % 2 == 0) {
		std::cout << "Expected extra ~ in save location string (did you forget last /?)." << std::endl;
		throw("Unsafe File Format");
	}
	else {
		for (int i = 0; i < nflds->size() / 2; i++) {
			for (int j = 0; j < nflds->at(i * 2).length(); j++)
				fol[j + fi] = nflds->at(i * 2).at(j);
			fi += nflds->at(i * 2).length();
			//tempstr << nflds->at(i * 2);
			std::string tstr = std::to_string(std::lrint(parseVal(nflds->at(i * 2 + 1))));
			for (int j = 0; j < tstr.length(); j++)
				fol[j + fi] = tstr.at(j);
			fi += tstr.length();
			//tempstr << (int)parseVal(nflds->at(i * 2 + 1));
		}
		for (int j = 0; j < nflds->back().length(); j++)
			fol[j + fi] = nflds->back().at(j);
		fi += nflds->back().length();
		fol[fi] = '\0';
		//tempstr << nflds->back();
	}
	//std::string temp = tempstr.str();
	//fol = new char[temp.length()];
	//std::strcpy(fol, temp.c_str());
	OSSpecificFuncs::createFolder(fol);
	return 1;
}

void ThreadParser::printSuccessPot(int potIdx, std::vector<double> params) {
	std::cout << "Successfully added " << potNames[potIdx] << "Potential with parameters: " << std::endl;
	for (int i = 0; i < potFields[potIdx].size(); i++)
		std::cout << potFields[potIdx][i] << " = " << params[i] << std::endl;
}

void ThreadParser::printSuccessMea(int meaIdx, std::vector<double> params) {
	std::cout << "Successfully added " << meaNames[meaIdx] << "Measurer with parameters: " << std::endl;
	for (int i = 0; i < meaFields[meaIdx].size(); i++)
		std::cout << meaFields[meaIdx][i] << " = " << params[i] << std::endl;
}

double ThreadParser::parseValMul(std::string str) {
	std::vector<std::string> * nflds = new std::vector<std::string>(0);
	parsingTools::split(str, nflds, '*');
	int x;
	double val = 1;
	for (int i = 0; i < nflds->size(); i++) {
		x = parsingTools::findMatch(&(nflds->at(i)), &varNames);
		if (x == -1) {
			try {
				val *= std::stod(nflds->at(i));
			}
			catch (std::invalid_argument) {
				std::cout << "Unrecognized variable: " << nflds->at(i) << std::endl;
				std::__throw_invalid_argument(nflds->at(i).c_str());
			}
		}
		else
			val *= var[x];
	}
	return val;
}

double ThreadParser::parseVal(std::string str) {
	if (str[0] == '\"') {
		inputText = str.substr(1);
		return 0;
	}
	std::vector<std::string> * nflds = new std::vector<std::string>(0);
	parsingTools::split(str, nflds, '+');
	double val = 0;
	for (int i = 0; i < nflds->size(); i++)
		val += parseValMul(nflds->at(i));
	return val;
}

std::vector<double> ThreadParser::getPotBlockParameters(int potNum) {
	std::vector<double> out = std::vector<double>(potFields[potNum].size(), 0);
	int fldNum = 0;
	do {
		std::getline(*fil, curLine);
		if (std::strstr(curLine.c_str(), "END"))
			break;
		parsingTools::split(curLine, flds, ' ');
		fldNum = parsingTools::findMatch(&(flds->at(0)), &potFields[potNum]);
		if (fldNum < 0) {
			std::cout << "Unknown parameter name: " << flds->at(0) << std::endl;
			return out;
		}
		else {
			out[fldNum] = parseVal(flds->at(1));
		}
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	return out;
}

std::vector<double> ThreadParser::getAbsBlockParameters(int absNum) {
	std::vector<double> out = std::vector<double>(absFields[absNum].size(), 0);
	int fldNum = 0;
	do {
		std::getline(*fil, curLine);
		if (std::strstr(curLine.c_str(), "END"))
			break;
		parsingTools::split(curLine, flds, ' ');
		fldNum = parsingTools::findMatch(&(flds->at(0)), &absFields[absNum]);
		if (fldNum < 0) {
			std::cout << "Unknown parameter name: " << flds->at(0) << std::endl;
			return out;
		}
		else {
			out[fldNum] = parseVal(flds->at(1));
		}
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	return out;
}

std::vector<double> ThreadParser::getMeaBlockParameters(int meaNum) {
	std::vector<double> out = std::vector<double>(meaFields[meaNum].size(), 0);
	if (out.size() < 1)
		return out;
	int fldNum = 0;
	do {
		std::getline(*fil, curLine);
		if (std::strstr(curLine.c_str(), "END"))
			break;
		parsingTools::split(curLine, flds, ' ');
		fldNum = parsingTools::findMatch(&(flds->at(0)), &meaFields[meaNum]);
		if (fldNum < 0) {
			std::cout << "Unknown parameter name: " << flds->at(0) << std::endl;
			return out;
		}
		else {
			out[fldNum] = parseVal(flds->at(1));
		}
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	return out;
}

std::vector<double> ThreadParser::getDenBlockParameters(int denNum) {
	std::vector<double> out = std::vector<double>(denFields[denNum].size(), 0);
	if (out.size() < 1)
		return out;
	int fldNum = 0;
	do {
		std::getline(*fil, curLine);
		if (std::strstr(curLine.c_str(), "END"))
			break;
		parsingTools::split(curLine, flds, ' ');
		fldNum = parsingTools::findMatch(&(flds->at(0)), &denFields[denNum]);
		if (fldNum < 0) {
			std::cout << "Unknown parameter name: " << flds->at(0) << std::endl;
			return out;
		}
		else {
			out[fldNum] = parseVal(flds->at(1));
		}
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	return out;
}

std::vector<double> ThreadParser::getRhoBlockParameters(int rhoNum) {
	std::vector<double> out = std::vector<double>(rhoFields[rhoNum].size(), 0);
	if (out.size() < 1)
		return out;
	int fldNum = 0;
	do {
		std::getline(*fil, curLine);
		if (std::strstr(curLine.c_str(), "END"))
			break;
		parsingTools::split(curLine, flds, ' ');
		fldNum = parsingTools::findMatch(&(flds->at(0)), &rhoFields[rhoNum]);
		if (fldNum < 0) {
			std::cout << "Unknown parameter name: " << flds->at(0) << std::endl;
			return out;
		}
		else {
			out[fldNum] = parseVal(flds->at(1));
		}
	} while (!(fil->eof()) && !std::strstr(curLine.c_str(), "END"));
	return out;
}