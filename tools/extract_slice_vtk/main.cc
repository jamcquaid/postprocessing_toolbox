#include "pptb.h"

using real_t = float;

int main(int argc, char** argv)
{
	// Parse input filename
	std::vector<std::string> args;
	for (int n = 0; n<argc; n++) args.push_back(std::string(argv[n]));
    std::string input_filename  = "surf.vtk";
	std::string output_filename = "slice_export.dat";
    if (args.size() > 1)
    {
        input_filename  = args[1];
		if (args.size()>2) output_filename = args[2];
    }
	
	// Call vtk file importer
	using geom_t = pptb::geom::surf_geom_t<real_t>;
	geom_t geom;
	pptb::io::import_vtk(input_filename, geom);
	
	// Setup slice operator
	const std::size_t num_slices = 7;
	const std::array<std::string, num_slices> slice_plane{"XY","XY","XY","XY","XY","XY","XY"};
	const std::array<real_t, num_slices> slice_pos{0.23926, 0.526372, 0.777595, 0.95704, 1.07667, 1.148448, 1.184337};
	const bool normalize_coord = true;
	const int normalize_dir    = 0;
	pptb::analysis::plane_slice_t<real_t, num_slices> plane_slice(slice_plane, slice_pos, normalize_coord, normalize_dir);

	// Run slice operation on specified geometry structure
	plane_slice.extract_data(geom);

	// Write data to file
	plane_slice.export_data(output_filename);
	
	return 0;
}
