#pragma once
#include <iomanip>

namespace pptb::io
{
	template <typename geom_t>
	inline void import_obj(const std::string& filename, geom_t& geom)
	{
		using rtype = geom_t::value_type;

		print("Importing obj file --> ",filename);

		// Open file
		spade::utils::ascii_file_t fh(filename);

		// Read first line
		fh.next_line();
		fh.next_line();
		std::string header    = "v";
		std::string firstchar = "v";

		rtype x,y,z;

		// Count nodes
		int nVerts = -1;
		while (firstchar == header)
		{
			fh.parse(firstchar, x,y,z);
			if (firstchar==header) ++nVerts;
			fh.next_line();
		}

		// Close and reopen
		fh.fh.close();
		spade::utils::ascii_file_t fh2(filename);

		// Skip first two lines
		fh2.next_line();
		fh2.next_line();

		// Allocate memory
		geom.nodes.resize(nVerts);

		for (int i = 0; i<nVerts; ++i)
		{
			fh2.parse(firstchar,geom.nodes[i][0],geom.nodes[i][1],geom.nodes[i][2]);
			fh2.next_line();
		}

		// Look for connectivity
		header    = "f";
		firstchar = "f";
		std::size_t n1,n2,n3;

		fh2.next_line();
		fh2.next_line();

		// Count faces
		int nElem = 0;
		while (!fh2.eof())
		{
			fh2.parse(firstchar,n1,n2,n3);
			if (firstchar==header) ++nElem;
			fh2.next_line();
		}
		
		// Allocate memory
		geom.connect.resize(nElem);
		geom.components.resize(nElem);

		// Close and reopen
		fh2.fh.close();
		spade::utils::ascii_file_t fh3(filename);

		fh3.next_line();
		fh3.next_line();
		for (int i = 0; i<nVerts; ++i) fh3.next_line();
		fh3.next_line();
		fh3.next_line();

		for (int i = 0; i<nElem; ++i)
		{
			fh3.parse(firstchar,n1,n2,n3);
			geom.connect[i][0] = n1-1;
			geom.connect[i][1] = n2-1;
			geom.connect[i][2] = n3-1;
			geom.components[i] = 0;
			fh3.next_line();
		}

		fh3.fh.close();
		
		print("Import fro complete! nElem, nVerts = ",nElem, nVerts);		
	}
}
