namespace pptb::io
{
	enum vtk_file_type
	{
		surface_vtk
	};

	template <typename rtype, typename geom_t, const vtk_file_type file_type = surface_vtk>
		inline void write_vtk(const std::string& filename, const geom_t& geom)
	{
		if constexpr(file_type == surface_vtk)
		{
			print("Writing surface VTK file...");

			// Open file
			std::ofstream fh(filename);

			// Write header
			fh << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << geom.nodes.size() << " double\n";

			for (std::size_t i = 0; i < geom.nodes.size(); ++i)
			{
				fh << geom.nodes[i][0] << " ";
				fh << geom.nodes[i][1] << " ";
				fh << geom.nodes[i][2] << "\n";
			}
			fh << "POLYGONS " << geom.connect.size() << " " << 4*geom.connect.size() << "\n";
			for (std::size_t i = 0; i < geom.connect.size(); ++i)
			{
				fh << "3 " << geom.connect[i][0]-1 << " " << geom.connect[i][1]-1 << " " << geom.connect[i][2]-1 << "\n";
			}
			fh << "CELL_DATA " << geom.connect.size() << "\n";
			fh << "SCALARS " << "Components" << " double\nLOOKUP_TABLE default\n";
			for (std::size_t i = 0; i < geom.connect.size(); ++i)
			{
				fh << geom.components[i] << "\n";
			}
		}
		else
		{
			print("INVALID VTK FILE TYPE SELECTED!");
		}
	}
}
