#include "stdafx.h"

namespace cfgParse {
	void readCFG(char *trg) {
		fil = new std::fstream(trg, std::ios::in);
		std::getline(*fil, curLine);
		if(std::strstr(curLine.c_str(), "MULTIPLE_SIMS_PARALLEL")){
			ParallelSimParser *sims = new ParallelSimParser(fil);
			std::cout << "Beginning global initializations..." << std::endl;
			sims->readConfig();
		}
	}
}

namespace parsingTools {
	void split(const std::string str, std::vector<std::string> * fields, char delim) {
		fields->clear();
		int startPos = 0;
		for (int i = 0; i < str.length(); i++) {
			if (str.at(i) == delim) {
				fields->push_back(str.substr(startPos, i-startPos));
				startPos = i + 1;
			}
		}
		if (startPos < str.length())
			fields->push_back(str.substr(startPos, str.length() - startPos));
		if (fields->size() < 1)
			fields->push_back("");
	}

	int findMatch(const std::string *cmp, const std::vector<std::string> *lst) {
		for (int i = 0; i < lst->size(); i++)
			if(std::strstr(cmp->c_str(), lst->at(i).c_str()))
				return i;
		return -1;
	}
}
