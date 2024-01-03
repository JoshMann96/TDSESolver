#include "stdafx.h"

int main(int argc, char** argv)
{
	if (argc != 2) {
		std::cout << "Needs one input argument: config file address." << std::endl;
		if (argc > 2)
			std::cout << "If address includes spaces, it needs to be surrounded by quotes." << std::endl;
	}
	else {
		MPI_Init(&argc, &argv);
		cfgParse::readCFG(argv[1]);
	}

	return 0;
}