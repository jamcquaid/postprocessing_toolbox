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

	// Slice settings
	scidf::node_t& node = input["Config"]["Slices"];
	std::vector<std::string> slice_name;
	std::vector<std::string> slice_plane;
	std::vector<real_t> slice_pos;
	const bool normalize_coord = input["Config"]["normalize_coord"];
	const int  normalize_dir   = input["Config"]["normalize_dir"];
	for (auto& p: node.children)
	{
		const std::string name    = p.first;
		scidf::node_t& node_local = p.second;

		// Slice settings
		std::string plane_tmp = node_local["slice_plane"];
		real_t pos_tmp        = node_local["position"];

		// Push into vector
		slice_name.push_back(name);
		slice_plane.push_back(plane_tmp);
		slice_pos.push_back(pos_tmp);
	}
	
	// Call vtk file importer
	using geom_t = pptb::geom::surf_geom_t<real_t>;
	geom_t geom;
	pptb::io::import_vtk(input_filename, geom);
	
	// Setup slice operator
	pptb::analysis::plane_slice_t<real_t> plane_slice(slice_name, slice_plane, slice_pos, normalize_coord, normalize_dir);

	// Run slice operation on specified geometry structure
	plane_slice.extract_data(geom);

	// Write data to file
	plane_slice.export_data(output_filename);
	
	return 0;
}
