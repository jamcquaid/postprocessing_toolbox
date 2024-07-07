#include "pptb.h"

using real_t = float;

int main(int argc, char** argv)
{
	// Parse input filename
	std::vector<int> args;
	for (int n = 0; n<argc; n++)
	{
		args.push_back(std::atoi(argv[n]));
	}
    std::size_t min_nt   = 0;
	std::size_t max_nt   = 1e9;
	std::size_t max_vars = 100;
    if (args.size() > 1)
    {
		min_nt = args[1];
		if (args.size() > 2) max_nt = args[2];
		if (args.size() > 3) max_vars = args[3];
    }

	// Parse directory
	pptb::utils::file_directory_t fileStruct;

	// Run time average
	pptb::analysis::timeAverage_t<real_t> timeAverage(min_nt, max_nt, max_vars, fileStruct);

	// Write output vtk
	timeAverage.export_data();
	
	return 0;
}
