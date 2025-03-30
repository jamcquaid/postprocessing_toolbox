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
	std::string input_filename  = "geom_in.vtk";
    if (argc > 1)
    {
        input_filename  = args[1];
    }

	// Read vtk file
	using geom_t = pptb::geom::surf_geom_t<real_t>;
	geom_t geom;
	pptb::io::import_vtk(input_filename, geom);

	// Compute loads
	pptb::analysis::integrate_loads(geom);
	
	return 0;
}
