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
	std::string output_filename = "geom_out.vtk";
    if (argc > 1)
    {
        input_filename  = args[1];
		if (argc > 2) output_filename = args[2];
    }

	// Read vtk file
	using geom_t = pptb::geom::surf_geom_t<real_t>;
	geom_t geom_in;
	pptb::io::import_vtk(input_filename, geom_in);
	
	// Rotate geometry
	pptb::analysis::transform_vtk<pptb::analysis::rotate>(geom_in, 2, 180.0);

	// Component lambda
	std::string varname = "Components";
	auto modify = [&](const auto& xyz, const auto& oldComp)
	{
		int newComp = 3;
		
		return newComp;
	};
	
	// Modify components
	pptb::analysis::modify_scalars(geom_in, varname, modify);
	
	// Write VTK file
	pptb::io::write_vtk_data<real_t>(output_filename, geom_in);
	
	return 0;
}
