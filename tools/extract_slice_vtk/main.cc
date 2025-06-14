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
	std::vector<std::string> slice_file;
	std::vector<pptb::analysis::truncate_t<real_t>> truncate;
	const bool normalize_coord = scidf::default_to<bool>(false, input["Config"]["normalize_coord"]);
	const int  normalize_dir   = scidf::default_to<int>(0, input["Config"]["normalize_dir"]);
	const real_t scaling       = scidf::default_to<real_t>(1.0, input["Config"]["slice_scaling"]);
	for (auto& p: node.children)
	{
		const std::string name    = p.first;
		scidf::node_t& node_local = p.second;

		// Slice settings
		std::string plane_tmp = node_local["slice_plane"];
		real_t pos_tmp = -10.0;
		std::string file_tmp = "";
		if (plane_tmp == "file" || plane_tmp == "3-points")
		{
			std::string test = node_local["file"];
			file_tmp = test;
		}
		else
		{
			pos_tmp = node_local["position"];
		}

		// Truncation settings
		pptb::analysis::truncate_t<real_t> settings{-1,0,0.0};
		auto& set_type = node_local["cutoff_type"];
		if (set_type.assigned_value)
		{
			const std::string type = node_local["cutoff_type"];
			if (type == "max") settings.type = 1;
			else if (type == "min") settings.type = 0;
			settings.dir  = node_local["cutoff_variable"];
			settings.val  = node_local["cutoff_value"];
		}
		truncate.push_back(settings);

		// Push into vector
		slice_name.push_back(name);
		slice_plane.push_back(plane_tmp);
		slice_pos.push_back(pos_tmp);
		slice_file.push_back(file_tmp);
	}
	
	// Call vtk file importer
	using geom_t = pptb::geom::surf_geom_t<real_t>;
	geom_t geom;
	pptb::io::import_vtk(input_filename, geom);
	
	// Setup slice operator
	pptb::analysis::plane_slice_t<real_t> plane_slice(slice_name, slice_plane, slice_pos, slice_file, truncate, scaling, normalize_coord, normalize_dir);

	// Run slice operation on specified geometry structure
	plane_slice.extract_data(geom);

	// Write data to file
	plane_slice.export_data(output_filename);
	
	return 0;
}
