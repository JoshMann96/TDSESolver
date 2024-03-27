#pragma once
#include "STLNCommonHeader.h"
#include "LinuxFuncs.h"
#include "ParsingTools.h"
#include "WfcRhoTools.h"
#include "AbsorptiveRegions.h"
#include "Potentials.h"
#include "SimulationManager.h"
#include "HHGFunctions.h"

//Reads simulation information after definition stage, executes simulation
class ThreadParser
{
private:
	std::stringstream *fil;
	//ProgressTracker *prg;
	std::vector<std::string> varNames;
	std::vector<double> var;
	SimulationManager *sim;
	std::string curLine;
	std::vector<std::string> *flds = new std::vector<std::string>(0);
	std::vector<Potentials::WaveFunctionSelfPotential*> *spc = new std::vector<Potentials::WaveFunctionSelfPotential*>(0);
	WfcToRho::Weight* wght;
	WfcToRho::Density* dens;
	double *x;
	char * fol;
	int n;
	std::string inputText = "";
	int* nElec = new int[1];
	//std::mutex* mtx;
	int mpiJob;
	
	//MPI update progress callback
	void mpiCallbackFunc(int prog);
	std::function<void(int)> callbackFunc;


	//Names of absorptive boundary types.
	std::vector<std::string> absNames =
	{
		"SmoothedRight",
		"SmoothedLeft"
	};
	//Add an absorptive function name here and then add its parameters in absFields.

	//Input fields for each abs boundary type.
	std::vector<std::vector<std::string>> absFields = 
	{
		{"SIZE", "RATE"},
		{"SIZE", "RATE"}
	};
	//Add the fields as an array. No need to specify the type, all are assumed to be double. Only one variable may be stored as a string.
	//Proceed to modify addAbsEdge in the CPP file.

	//Names of kinetic operators
	std::vector<std::string> kinNames =
	{
		"FreeElecPSM",
		"MathExprPSM",
		"EffMassBdyPSM",
		"MathExprBdyPSM"
	};

	//Input fields for each kinetic operator
	std::vector<std::vector<std::string>> kinFields =
	{
		{"M_EFF"},
		{"EXPR"},
		{"M_EFF_L", "M_EFF_R", "TRANS_RATE", "TRANS_POS", "EDGE_RATE", "EXP_ORDER", "FORCE_NORM"},
		{"EXPR", "EXP_ORDER", "FORCE_NORM"}
	};

	//Names of simulation density (post-process) calculators
	std::vector<std::string> rhoNames =
	{
		"Direct",
		"GaussSmooth"
	};

	//Input fields for each density (post-process) calculators
	std::vector<std::vector<std::string>> rhoFields =
	{
		{},
		{"SIGMA"}
	};

	//Names of simulation density calculators
	std::vector<std::string> denNames =
	{
		"FermiGas",
		"FromDOSFile"
	};

	//Input fields for each density calculator
	std::vector<std::vector<std::string>> denFields =
	{
		{"FERMI_E"},
		{"FERMI_E", "FERMI_LEVEL", "WELL_LENGTH", "FILE"}
	};

	//Names of potential functions
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
		"GaussFileField", 
		"GaussianPulseCylindrical",
		"FiniteBox",
		"FullCylindricalSpaceCharge",
		"UniformSurfaceSpaceCharge",
		"LinearBulkCylindricalFieldSpaceCharge",
		"LinearBulkCylSectionFieldSpaceCharge",
		"DielectricBulkCylindricalFieldSpaceCharge",
		"OhmicRetardingCurrent",
		"CylindricalImageCharge",
		"GaussianPulseExponential",
		"Cos2FileField",
		"Cos2PulseCylindrical",
		"Cos2PulseExponential",
		"CutoffCos2PulseCylindrical"
	};
	//Add a potential function name here and then add its parameters in potFields.

	//Input fields for each potential function
	std::vector<std::vector<std::string>> potFields =
	{
		{"CENTER", "FERMI_E", "WORK_F", "BACK_START", "BACK_WIDTH", "REF_POINT"}, //BackedJell
		{"CENTER", "FERMI_E", "WORK_F", "REF_POINT"}, //Jellium
		{"CENTER", "LATTICE_SPACING", "Z_PROTONS", "DECAY_CONST"}, //ShieldedAtomicPotentiakl
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "ENHANCEMENT", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT", "BUFFER_LENGTH" }, //CutoffGaussianPulseCylindrical
		{"MIN_X", "MAX_X", "MAX_E", "LAM", "TAU", "PHASE", "PEAK_T", "BUFFER_TIME", "EPSILON_REL_REAL", "EPSILON_REL_IMAG", "CONDUCTIVITY", "REF_POINT"}, //BasicPenetratingField
		{"MIN_X", "MAX_X", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"}, //GaussianPulseFlat
		{"OFFSET", "FILE", "REF_POINT"}, //FilePot
		{"STRENGTH", "PERP_DIST", "REF_POINT"}, //SelfCoulomb
		{"MIN_X", "MAX_X", "BUFFER_X_MIN", "BUFFER_X_MAX", "E_FIELD", "T_START", "BUFFER_TIME", "REF_POINT"}, //BiasField
		{"OFFSET", "FILE", "BUFFER_X_MIN", "BUFFER_X_MAX", "BUFFER_LENGTH", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"}, //GaussFileField
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "ENHANCEMENT", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"}, //GaussianPulseCylindrical
		{"MIN_X", "MAX_X", "V_IN", "REF_POINT"}, //FiniteBox
		{"FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "REF_POINT"}, //FullCylindricalSpaceCharge
		{"FERMI_E", "MIN_X", "MAX_X", "REF_POINT"}, //UniformSurfaceSpaceCharge
		{"FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "REF_POINT"}, //LinearBulkCylindricalFieldSpaceCharge
		{"FERMI_E", "RADIUS", "SURF_X", "THETA_0", "REF_POINT"}, //LinearBulkCylSectionFieldSpaceCharge
		{"FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "WELL_WIDTH", "DAMP_RATE", "REF_POINT"}, //DielectricBulkCylindricalFieldSpaceCharge
		{"SURF_X", "TRANS_LEN", "RESISTIVITY", "REF_POINT"}, //OhmicRetardingCurrent
		{"WORK_F", "FERMI_E", "RADIUS", "SURF_X", "MIN_X", "MAX_X", "REF_POINT"}, //CylindricalImageCharge
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "BUFFER_TIME", "REF_POINT"}, //GaussianPulseExponential
		{"OFFSET", "FILE", "BUFFER_X_MIN", "BUFFER_X_MAX", "BUFFER_LENGTH", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "REF_POINT"}, //Cos2FileField
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "ENHANCEMENT", "LAM", "TAU", "PEAK_T", "PHASE", "REF_POINT"}, //Cos2PulseCylindrical
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "LAM", "TAU", "PEAK_T", "PHASE", "REF_POINT"}, //Cos2PulseExponential
		{"MIN_X", "MAX_X", "RADIUS", "MAX_E", "ENHANCEMENT", "LAM", "TAU", "PEAK_T", "PHASE", "REF_POINT", "BUFFER_LENGTH" } //CutoffCos2PulseCylindrical
	};
	//Add the fields as an array. No need to specify the type, all are assumed to be double. Only one variable may be stored as a string.
	//Proceed to modify addPotential in the CPP file.

	//Names of measurers
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

	//Input fields for each measurer
	std::vector<std::vector<std::string>> meaFields =
	{
		{"NAME"}, //BasicMeasurers
		{"N_X", "N_T"}, //VFuncT
		{"NAME", "VD_NUM", "VD_POS"}, //VDPsi
		{"NAME", "VD_NUM", "VD_POS"}, //VDProbCurrent
		{}, //TotProb
		{}, //ExpectA
		{}, //ExpectP
		{}, //ExpectX
		{}, //ExpectE0
		{"N_X", "N_T"}, //Psi2T
		{}, //V0
		{}, //TS
		{}, //XS
		{}, //DT
		{}, //DX
		{}, //NSteps
		{}, //NPts
		{"NAME"}, //Header
		{"FIL_NAME", "CONST"}, //DoubleConst
		{"N_X", "N_P", "N_T", "P_MIN", "P_MAX", "WINDOW_WIDTH"}, //WignerQPD
		{"NAME", "VD_NUM", "MEA_T"}, //PsiT
		{"NAME", "VD_NUM", "MEA_T"}, //PotT
		{"NAME", "VD_NUM", "VD_POS"}, //VDPot
		{}, //ExpectE
		{"NAME", "VD_NUM_0", "VD_NUM_1", "VD_POS"}, //VDBiPsi
		{"NAME", "VD_NUM_0", "VD_NUM_1", "VD_POS"}, //VDBiPot
		{"NAME", "VD_NUM", "VD_POS", "E_MAX", "N_SAMP"}, //VDFluxSpec
		{} //WfcRhoWeights
	};
	//Add the fields as an array. No need to specify the type, all are assumed to be double. Only one variable may be stored as a string.
	//Proceed to modify addMeasurer in the CPP file.

public:
	ThreadParser(std::stringstream *fil, std::vector<std::string> varNames, std::vector<double> var, int mpiJob);
	~ThreadParser();

	//execute rest of file
	void executeScript();
	//read/execute line
	int readCommand();
	//initialize simulation
	int generalSimInit();
	//specifically for HHG (obsolete)
	int hhgSimInit();
	//set data save location
	int setSaveLoc(std::string input);
	//add potential element to simulation
	int addPotential(std::string input);
	//add measurer element to simulation
	int addMeasurer(std::string input);
	//add absorptive element boundary to simulation
	int addAbsEdge(std::string input);
	//sets weight calculator (for density calculation)
	int setWghtCalc(std::string input);
	//sets density calculator
	int setRhoCalc(std::string input);
	//sets kinetic operator
	int setKin(std::string input);
	//parse a value which may include addition (+) and multiplication (*) (by splitting about addition first)
	double parseVal(std::string str);
	//parse a value which may include multiplication (*), for use by parseVal
	double parseValMul(std::string str);

	//gets all field values for an element
	std::vector<double> getBlockParameters(int num, std::vector<std::vector<std::string>> fields);

	//printing functions (not in use)
	void printSuccessPot(int potIdx, std::vector<double> params);
	void printSuccessMea(int meaIdx, std::vector<double> params);
};

