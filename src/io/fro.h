#pragma once
#include <iomanip>

namespace pptb::io
{
	template <typename geom_t>
	inline void import_fro(const std::string& filename, geom_t& geom)
	{
		using rtype = geom_t::value_type;

		print("Importing fro file --> ",filename);

		// Open file
		spade::utils::ascii_file_t fh(filename);

		// Read first line
		fh.next_line();
		std::size_t nElem, nVerts, n0, n1, n2, n3;
		fh.parse(nElem, nVerts, n0, n1, n2, n3);

		// Allocate memory
		geom.nodes.resize(nVerts);
		geom.connect.resize(nElem);
		geom.components.resize(nElem);
		
		// Sweep node list
		for (int n = 0; n<nVerts; n++)
		{
			// Next line
			fh.next_line();

			// Data types
			std::size_t ielem;
			auto& node = geom.nodes[n];

			// Parse this line
			fh.parse(ielem, node[0], node[1], node[2]);
		}

		// Import connectivity and component numbers
		for (int n = 0; n<nElem; n++)
		{
			// Next line
			fh.next_line();

			// Data types
			std::size_t ielem;
			auto& connect = geom.connect[n];
			auto& comp    = geom.components[n];

			// Parse line
			fh.parse(ielem, connect[0], connect[1], connect[2], comp);

			// Switch to 0-based system
			connect[0]--;
			connect[1]--;
			connect[2]--;
		}

		fh.fh.close();
		
		print("Import fro complete! nElem, nVerts = ",nElem, nVerts);
		
	}

	template <typename real_t, typename geom_t>
	inline void write_fro(const std::string& filename, const geom_t& geom)
	{
		print("Writing surface fro file...");

		// Open file
		print(filename);
		std::ofstream fh(filename);
		fh << std::fixed << std::setprecision(10);

		// Write header
		const int n0 = 0, n1 = 0, n2 = 0;
		int maxCompID = 0, iv = 0;
		if (geom.components.size() == geom.connect.size())
		{
			for (int i = 0; i<geom.connect.size(); ++i) maxCompID = std::max(maxCompID, geom.components[i]);
		}
		else
		{
			for (int v = 0; v<geom.varnames.size(); ++v)
			{
				if (geom.varnames[v] == "Components")
				{
					iv = v;
					for (int i = 0; i<geom.connect.size(); ++i)	maxCompID = std::max(maxCompID, int(geom.data[v][i]));
				}
			}
		}
		fh << geom.connect.size() << " " << geom.nodes.size() << " " << n0 << " " << n1 << " " << n2 << " " << maxCompID << "\n";

		// Write out node data
		std::cout << std::setprecision(16);
		for (int i = 0; i<geom.nodes.size(); ++i)
		{
			fh << i+1 << " " << geom.nodes[i][0] << " " << geom.nodes[i][1] << " " << geom.nodes[i][2] << "\n";
		}

		// Write out connectivity
		for (int i = 0; i<geom.connect.size(); ++i)
		{
			fh << i+1 << " " << geom.connect[i][0]+1 << " " << geom.connect[i][1]+1 << " " << geom.connect[i][2]+1 << " " << geom.data[iv][i] << "\n";
		}
		
		// Close file
		fh.close();
	}
}
