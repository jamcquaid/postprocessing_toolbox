
namespace pptb::io
{
	template <typename geom_t>
	inline void import_fro(const std::string& filename, geom_t& geom)
	{
		using rtype = geom_t::value_type;

		print("Importing fro file --> ",filename);

		// Open file
		utils::ascii_file_t fh(filename);

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
		}

		print("Import fro complete! nElem, nVerts = ",nElem, nVerts);
		
	};
}
