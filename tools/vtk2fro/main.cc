#include "pptb.h"

using real_t = double;

int main(int argc, char** argv)
{
	// Parse input filename
	std::vector<std::string> args;
	for (int n = 0; n<argc; n++) args.push_back(std::string(argv[n]));
    std::string input_filename  = "geom.vtk";
	std::string output_filename = "geom.fro";
    if (argc > 1)
    {
        input_filename  = args[1];
		if (argc > 2) output_filename = args[2];
    }

	// Call VTK file importer
	using geom_t = pptb::geom::surf_geom_t<real_t>;
	geom_t geom;
	pptb::io::import_vtk(input_filename, geom);

	// Write fro file
	pptb::io::write_fro<real_t>(output_filename, geom);
	
	// Set 
	return 0;
}
