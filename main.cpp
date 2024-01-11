#include "stdafx.h"

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cout << "Takes two arguments:" << std::endl;
		std::cout << "    config file (relative or full address, required)" << std::endl;
		std::cout << "    root test (am I the MPI root process or not?)" << std::endl;
	}
	else{
		MPI_Init(&argc, &argv);
		if(argc >= 3){ //is manual root selection being used?
			if(std::stoi(argv[2])==1){ //am I the root process?
				MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Root_Proc); //get the root process (my rank)
				std::cout << "MPI Root ID: " << MPI_Root_Proc << std::endl;
				std::cout << "Root has " << omp_get_max_threads() << " OMP threads available." << std::endl;
				int size;
				MPI_Comm_size(MPI_COMM_WORLD, &size);
				for(int i = 0; i < size; i++){ //loop over all (other) processes
					if(i != MPI_Root_Proc)
						MPI_Ssend(&MPI_Root_Proc, 1, MPI_INT, i, MPITag::AmRoot, MPI_COMM_WORLD); //send root process ID
				}
			}
			else if (std::stoi(argv[2])==0){ //am I not the root process
				MPI_Status stat;
				MPI_Recv(&MPI_Root_Proc, 1, MPI_INT, MPI_ANY_SOURCE, MPITag::AmRoot, MPI_COMM_WORLD, &stat); //listen for root process
				std::cout << "Non-root has " << omp_get_max_threads() << " OMP threads available." << std::endl;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD); //wait until all processes have received root
		//everyone has figured out who is root. If a second directive was not provided then root process is still 0

		cfgParse::readCFG(argv[1]);
		MPI_Finalize();
	}

	return 0;
}