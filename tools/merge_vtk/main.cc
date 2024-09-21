#include "pptb.h"

using real_t = double;

int main(int argc, char** argv)
{
	// Parse input filename
	std::vector<std::string> args;
	for (int n = 0; n<argc; n++)
	{
		args.push_back(std::string(argv[n]));
	}
	std::string output_filename = "geom_merge.vtk";
    if (argc > 1)
    {
        output_filename  = args[1];
    }

	// Parse directory
	pptb::utils::file_directory_t fileStruct;
	
	// Modify components
	const auto geom_merge = pptb::analysis::transform_vtk<pptb::analysis::merge, real_t>(fileStruct);
	
	// Write VTK file
	pptb::io::write_vtk_data<real_t>(output_filename, geom_merge);
	
	return 0;
}
