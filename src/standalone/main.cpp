#include "STLNCommonHeader.h"
#include "SimulationWorkspace.h"


int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cout << "Takes two positional arguments:" << std::endl;
		std::cout << "    config file (required)" << std::endl;
		std::cout << "    root test (optional, 1 for root, 0 for calculation process)" << std::endl;
		return 1;
	}
	MPI_Init(&argc, &argv);

	//manual root process selection
	if(argc >= 3){
		if(std::stoi(argv[2])==1){ //I am the root process?
			MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Root_Proc); //set the root process (my rank)

			std::cout << "MPI Root ID: " << MPI_Root_Proc << std::endl;
			std::cout << "Root has " << omp_get_max_threads() << " OMP threads available." << std::endl;

			//get num processes
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);

			//tell all other processes the root process rank, wait for each one to be received
			for(int i = 0; i < size; i++){ 
				if(i != MPI_Root_Proc)
					MPI_Ssend(&MPI_Root_Proc, 1, MPI_INT, i, MPITag::AmRoot, MPI_COMM_WORLD);
			}
		}
		else if (std::stoi(argv[2])==0){ //I am not the root process
			//listen for a process to declare root
			MPI_Status stat;
			MPI_Recv(&MPI_Root_Proc, 1, MPI_INT, MPI_ANY_SOURCE, MPITag::AmRoot, MPI_COMM_WORLD, &stat);

			std::cout << "Non-root has " << omp_get_max_threads() << " OMP threads available." << std::endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD); //wait until all processes have received root

	//everyone has figured out who is root. If a second directive was not provided then root process is still 0

	//read and execute script
	std::fstream* fil = new std::fstream(argv[1], std::ios::in);
	SimulationWorkspace *sims = new SimulationWorkspace(fil);
	sims->executeScript();
	
	//all done
	MPI_Finalize();

	return 0;
}