#include "stdafx.h"

namespace OSSpecificFuncs {
	void createFolder(const char* fol) {
		
		boost::filesystem::create_directories(fol);

		/*
		struct stat info;
		if (stat(fol, &info) != 0)
			if(mkdir(fol, 0755) == -1)
				std::cerr << "Error :  " << fol << " : " << strerror(errno) << std::endl;
				*/
		/*{
			int l1 = 9, l2 = std::strlen(fol);

			char* cstr = new char[l1 + l2 + 1];
			strncpy(cstr, "mkdir -p ", l1);
			strcpy(&cstr[l1], fol);
			std::cout << cstr << std::endl;

			std::system(cstr);
			delete[] cstr;
		}*/
	}
}