#pragma once
class ThreadParser
{
private:
	std::stringstream *fil;
	ProgressTracker *prg;
	std::vector<std::string> varNames;
	std::vector<double> var;
	MultiSimulationManager *sim;
	std::string curLine;
	std::vector<std::string> *flds = new std::vector<std::string>(0);
	std::vector<Potentials::WaveFunctionSelfPotential*> *spc = new std::vector<Potentials::WaveFunctionSelfPotential*>(0);
	WfcToRho::Weight* wght;
	WfcToRho::Density* dens;
	double *x;
	char * fol;
	int n;
	int simIdx;
	int multiElec = 0;
	std::string inputText = "";
	int* nelec = new int[1];
	std::mutex* mtx;

	std::vector<std::string> absNames =
	{
		"SmoothedRight",
		"SmoothedLeft"
	};
	//Add an absorptive function name here and then add its parameters in absFields.

	std::vector<std::vector<std::string>> absFields = 
	{
		{"SIZE", "RATE"},
		{"SIZE", "RATE"}
	};
	//Add the fields as an array. No need to specify the type, all are assumed to be double. Only one variable may be stored as a string.
	//Proceed to modify addAbsEdge in the CPP file.

	std::vector<std::string> kinNames =
	{
		"FreeElecPSM",
		"MathExprPSM",
		"EffMassBdyPSM",
		"MathExprBdyPSM"
	};

	std::vector<std::vector<std::string>> kinFields =
	{
		{"M_EFF"},
		{"EXPR"},
		{"M_EFF_L", "M_EFF_R", "TRANS_RATE", "TRANS_POS", "EDGE_RATE", "EXP_ORDER", "FORCE_NORM"},
		{"EXPR", "EXP_ORDER", "FORCE_NORM"}
	};

	std::vector<std::string> rhoNames =
	{
		"Direct",
		"GaussSmooth"
	};

	std::vector<std::vector<std::string>> rhoFields =
	{
		{},
		{"SIGMA"}
	};

	std::vector<std::string> denNames =
	{
		"FermiGas",
		"FromDOSFile"
	};

	std::vector<std::vector<std::string>> denFields =
	{
		{"FERMI_E"},
		{"FERMI_E", "FERMI_LEVEL", "WELL_LENGTH", "FILE"}
	};

	std::vector<std::string> potNames = 
	{
		"BackedJell",
		"Jellium",
		"ShieldedAtomicPotential",
		"CutoffGaussianPulseCylindrical",
		"BasicPenetratingField",
		"GaussianPulseFlat",
		"FilePot",
		"SelfCoulomb",
		"BiasField",
		"FileField", 
		"GaussianPulseCylindrical",
		"FiniteBox",
		"FullCylindricalSpaceCharge",
		"UniformSurfaceSpaceCharge",
		"LinearBulkCylindricalFieldSpaceCharge",
		"LinearBulkCylSectionFieldSpaceCharge",
		"DielectricBulkCylindricalFieldSpaceCharge",
		"OhmicRetardingCurrent",
		"CylindricalImageCharge"
	};
	//Add a potential function name here and then add its parameters in potFields.

	std::vector<std::vector<std::string>> potFields =
	{
		{"CENTER", "FERMI_E", "WORK_F", "BACK_START", "BACK_WIDTH", "REF_POINT"},
		{"CENTER", "FERMI_E", "WORK_F", "REF_POINT"},
		{"CENTER", "LATTICE_SPACING", "Z_PROTONS", "DECAY_CONST"},
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "ENHANCEMENT", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT", "BUFFER_LENGTH" },
		{"MIN_X", "MAX_X", "MAX_E", "LAM", "TAU", "PHASE", "PEAK_T", "BUFFER_TIME", "EPSILON_REL_REAL", "EPSILON_REL_IMAG", "CONDUCTIVITY", "REF_POINT"},
		{"MIN_X", "MAX_X", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"},
		{"OFFSET", "FILE", "REF_POINT"},
		{"STRENGTH", "PERP_DIST", "REF_POINT"},
		{"MIN_X", "MAX_X", "BUFFER_X_MIN", "BUFFER_X_MAX", "E_FIELD", "T_START", "BUFFER_TIME", "REF_POINT"},
		{"OFFSET", "FILE", "BUFFER_X_MIN", "BUFFER_X_MAX", "BUFFER_LENGTH", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"},
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "ENHANCEMENT", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"},
		{"MIN_X", "MAX_X", "V_IN", "REF_POINT"},
		{"FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "REF_POINT"},
		{"FERMI_E", "MIN_X", "MAX_X", "REF_POINT"},
		{"FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "REF_POINT"},
		{"FERMI_E", "RADIUS", "SURF_X", "THETA_0", "REF_POINT"},
		{"FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "WELL_WIDTH", "DAMP_RATE", "REF_POINT"},
		{"SURF_X", "TRANS_LEN", "RESISTIVITY", "REF_POINT"},
		{"WORK_F", "FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "REF_POINT"},
	};
	//Add the fields as an array. No need to specify the type, all are assumed to be double. Only one variable may be stored as a string.
	//Proceed to modify addPotential in the CPP file.

	std::vector<std::string> meaNames =
	{
		"BasicMeasurers",
		"VFuncT",
		"VDPsi",
		"VDProbCurrent",
		"TotProb",
		"ExpectA",
		"ExpectP",
		"ExpectX",
		"ExpectE0",
		"Psi2T",
		"V0",
		"TS",
		"XS",
		"DT",
		"DX",
		"NSteps",
		"NPts",
		"Header",
		"DoubleConst",
		"WignerQPD",
		"PsiT",
		"PotT",
		"VDPot",
		"ExpectE",
		"VDBiPsi",
		"VDBiPot",
		"VDFluxSpec",
		"WfcRhoWeights"
	};
	//Add a measurement name here and then add its parameters in potFields.

	std::vector<std::vector<std::string>> meaFields =
	{
		{"NAME"},
		{"N_X", "N_T"},
		{"NAME", "VD_NUM", "VD_POS"},
		{"NAME", "VD_NUM", "VD_POS"},
		{},
		{},
		{},
		{},
		{},
		{"N_X", "N_T"},
		{},
		{},
		{},
		{},
		{},
		{},
		{},
		{"NAME"},
		{"FIL_NAME", "CONST"},
		{"N_X", "N_P", "N_T", "P_MIN", "P_MAX", "WINDOW_WIDTH"},
		{"NAME", "VD_NUM", "MEA_T"},
		{"NAME", "VD_NUM", "MEA_T"},
		{"NAME", "VD_NUM", "VD_POS"},
		{},
		{"NAME", "VD_NUM_0", "VD_NUM_1", "VD_POS"},
		{"NAME", "VD_NUM_0", "VD_NUM_1", "VD_POS"},
		{"NAME", "VD_NUM", "VD_POS", "E_MAX", "N_SAMP"},
		{}
	};
	//Add the fields as an array. No need to specify the type, all are assumed to be double. Only one variable may be stored as a string.
	//Proceed to modify addMeasurer in the CPP file.

public:
	ThreadParser(std::stringstream *fil, ProgressTracker *prg, std::vector<std::string> varNames, std::vector<double> var, int simIdx, std::mutex* mtx);
	~ThreadParser();
	void readConfig();
	int readCommand();
	int generalSimInit();
	int hhgSimInit();
	int setSaveLoc(std::string input);
	int addPotential(std::string input);
	int addMeasurer(std::string input);
	int addAbsEdge(std::string input);
	int setWghtCalc(std::string input);
	int setRhoCalc(std::string input);
	int setKin(std::string input);
	double parseVal(std::string str);
	double parseValMul(std::string str);

	std::vector<double> getBlockParameters(int num, std::vector<std::vector<std::string>> fields);

	void printSuccessPot(int potIdx, std::vector<double> params);
	void printSuccessMea(int meaIdx, std::vector<double> params);
};

