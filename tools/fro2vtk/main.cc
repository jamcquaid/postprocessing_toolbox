#include "pptb.h"

using real_t = float;

int main(int argc, char** argv)
{
	// Parse input filename
	std::vector<std::string> args;
	for (int n = 0; n<argc; n++) args.push_back(std::string(argv[n]));
    std::string input_filename  = "geom.fro";
	std::string output_filename = "geom.vtk";
    if (args.size() == 2)
    {
        input_filename  = args[1];
		if (args.size() == 3) output_filename = args[2];
    }

	// Call fro file importer
	using geom_t = pptb::io::surf_geom_t<real_t>;
	geom_t geom;
	pptb::io::import_fro(input_filename, geom);

	// Write VTK file
	pptb::io::write_vtk<real_t, geom_t, pptb::io::surface_vtk>(output_filename, geom);
	
	// Set 
	return 0;
}
