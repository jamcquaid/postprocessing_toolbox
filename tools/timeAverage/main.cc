#include "pptb.h"
#include "scidf.h"

using real_t = float;

int main()
{
	// input filename
	const std::string ifile = "input.sdf";
	scidf::node_t input;
	scidf::read(ifile, input);
	
    std::size_t min_nt   = input["Config"]["min_nt"];
	std::size_t max_nt   = input["Config"]["max_nt"];

	// Parse directory
	pptb::utils::file_directory_t fileStruct;

	// Run time average
	pptb::analysis::timeAverage_t<real_t> timeAverage(min_nt, max_nt, input["Config"]["out_file"], fileStruct);

	// Write output vtk
	timeAverage.export_data();
	
	return 0;
}
