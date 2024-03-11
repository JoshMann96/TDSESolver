#pragma once
#include "CORECommonHeader.h"

namespace parsingTools {
	void split(const std::string str, std::vector<std::string> * fields, char delim);
	int findMatch(const std::string *cmp, const std::vector<std::string> *lst);
}