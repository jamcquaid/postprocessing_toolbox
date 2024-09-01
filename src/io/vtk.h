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
				fh << "3 " << geom.connect[i][0] << " " << geom.connect[i][1] << " " << geom.connect[i][2] << "\n";
			}
			fh << "CELL_DATA " << geom.connect.size() << "\n";
			fh << "SCALARS " << "Components" << " double\nLOOKUP_TABLE default\n";
			for (std::size_t i = 0; i < geom.connect.size(); ++i)
			{
				fh << geom.components[i] << "\n";
			}
			
			// Close file
			fh.close();
			
		}
		else
		{
			print("INVALID VTK FILE TYPE SELECTED!");
		}
	}

	template <typename rtype, typename geom_t, const vtk_file_type file_type = surface_vtk>
	inline void write_vtk_data(const std::string& filename, const geom_t& geom, const std::size_t& maxVars = 10000)
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
			if (geom.maxNpe>2)
			{
				fh << "POLYGONS " << geom.connect.size() << " " << (geom.maxNpe+1)*geom.connect.size() << "\n";
			}
			else
			{
				fh << "LINES " << geom.connect.size() << " " << (geom.maxNpe+1)*geom.connect.size() << "\n";
			}
			for (std::size_t i = 0; i < geom.connect.size(); ++i)
			{
				fh << geom.maxNpe << " ";
				for (int j = 0; j<geom.maxNpe; ++j)
				{
					fh<<geom.connect[i][j]<<" ";
				}
				fh<<std::endl;
			}
			fh << "CELL_DATA " << geom.connect.size() << "\n";

			for (int ivar = 0; ivar<std::min(geom.varnames.size(), maxVars); ivar++)
			{
				fh << "SCALARS " << geom.varnames[ivar] << " double\nLOOKUP_TABLE default\n";
				for (std::size_t i = 0; i < geom.connect.size(); ++i)
				{
					fh << geom.data[ivar][i] << "\n";
				}
			}
			
			// Close file
			fh.close();
			
		}
		else
		{
			print("INVALID VTK FILE TYPE SELECTED!");
		}
	}

	template <typename geom_t>
	inline void import_vtk(const std::string& filename, geom_t& geom, const bool& verbose = true)
	{
		using rtype = geom_t::value_type;

		print("Importing vtk file --> ",filename);

		// Open file
		utils::ascii_file_t fh(filename);

		// Skip header lines
		fh.next_line();
		fh.next_line();
		fh.next_line();
		fh.next_line();
		fh.next_line();
		std::size_t nVerts;
		std::string str1, str2;
		fh.parse(str1, nVerts, str2);

		// Initialize memory for nodes
		geom.nodes.resize(nVerts);

		// Import nodes
		for (int n = 0; n<nVerts; n++)
		{
			// Read line
			fh.next_line();

			// Node vector
			auto& nodes = geom.nodes[n];

			// Parse input
			fh.parse(nodes[0], nodes[1], nodes[2]);
		}

		// Import number of elements
		std::size_t nElem, n0;
		fh.next_line();
		fh.parse(str1, nElem, n0);

		// Store memory for elements
		geom.connect.resize(nElem);

		// Import connectivity
		for (int n = 0; n<nElem; n++)
		{
			// Read line
			fh.next_line();

			// Node vector
			auto& connect = geom.connect[n];

			// Parse input
			std::size_t tmp;
			fh.parse(tmp, connect[0], connect[1], connect[2]);
		}

		// Skip line
		fh.next_line();
		
		// Find number of variables
		int nvars = 0;
		while (!fh.eof())
		{
			// Skip headers
			fh.next_line();
			if (fh.eof()) break;
			fh.next_line();

			// Count variables
			nvars++;

			// Skip stored data
			for (int n = 0; n<nElem; n++)
			{
				fh.next_line();
			}
		}

		// Close file
		fh.fh.close();

		// Print grid size
		if (verbose)
		{
			print("   nElem    = ",nElem);
			print("   nVerts   = ",nVerts);
		}
		
		// Reopen file
		utils::ascii_file_t file(filename);

		// Skip lines we've already read
		for (int n = 0; n<5+nVerts+1+nElem+1; n++) file.next_line();
		
		// Allocate memory for variable data
		geom.varnames.resize(nvars);
		geom.data.resize(nvars);
		for (int v = 0; v<nvars; v++)
		{
			// Read in variable name
			file.next_line();
			file.parse(str1, geom.varnames[v], str2);
			file.next_line();

			// Print to screen
			if (verbose) print("   Var, varname = ",v,geom.varnames[v]);

			// Allocate memory
			auto& data_loc = geom.data[v];
			data_loc.resize(nElem);

			// Import data
			for (int n = 0; n<nElem; n++)
			{
				// Read data
				file.next_line();
				file.parse(data_loc[n]);
			}
		}

		// Close file
		file.fh.close();

		if (verbose) print("Import complete!");
		
	}
}
