#include "pptb.h"
#include "scidf.h"

using real_t = double;

int main(int argc, char** argv)
{
	// input filename
	const std::string ifile = "input.sdf";
	scidf::node_t input;
	scidf::read(ifile, input);

	std::string input_filename  = input["Config"]["input_filename"];
	std::string output_filename = input["Config"]["output_filename"];
	const int nsmooth           = input["Config"]["nsmooth"];

	// Import VTK with spade
	print("Importing geometry...");
	spade::geom::vtk_geom_t<3, 3, real_t> geom;
	spade::geom::read_geom(input_filename, geom, true);

	// Conversion
	spade::geom::voronoi_geom_t<3, real_t> voronoi;
	voronoi = pptb::analysis::dual_conversion(nsmooth, geom);
	
	// Write output voronoi tesselation
	pptb::io::write_vtk_data(output_filename, voronoi);
	
	return 0;
}
