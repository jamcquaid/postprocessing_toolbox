namespace pptb::analysis
{
	template <typename rtype>
	struct plane_t
	{
		using value_type = rtype;
		
		// Member variables
		std::array<value_type, 3> nvec; // Plane normal
		std::array<value_type, 3> tvec; // Plane tangential vector 1
		std::array<value_type, 3> mvec; // plane tangential vector 2
		std::array<value_type, 3> pt;
	};
	
	template <typename rtype, const std::size_t num_slices>
	struct plane_slice_t
	{
		using value_type  = rtype;
		using surf_geom_t = geom::surf_geom_t<value_type>;
		using vec_t = std::array<value_type, 3>;
		
		// Member variables
		std::array<surf_geom_t, num_slices> slice_data; // Cross-section geometry and data
		std::array<std::size_t, num_slices> slice_plane; // Slice plane
		std::array<value_type, num_slices> slice_pos; // Slice position in plane normal direction
		std::array<plane_t<value_type>, num_slices> planes; // Slice plane basis
		bool normalize_coord = false;
		int normalize_dir    = 0;
		
		// Some functions
		constexpr static int nslice(){return num_slices;}

		// Get coordinate on plane
		inline vec_t get_plane_pt(const auto& iplane)
		{
			vec_t pt{0.0, 0.0, 0.0};
			pt[slice_plane[iplane]] = slice_pos[iplane];
			return pt;
		}
		
		// Constructor --> slice select based on Cartesian planes
	    plane_slice_t(const std::array<std::string, num_slices>& slice_plane_in, const auto& slice_pos_in, const bool& norm_coord_in = false, const int& norm_dir_in = 0)
			: slice_pos{slice_pos_in}, normalize_coord{norm_coord_in}, normalize_dir{norm_dir_in}
		{
			print("Initializing slice operator...");

			for (int n = 0; n<nslice(); n++)
			{
				if (slice_plane_in[n] == "XY" || slice_plane_in[n] == "YX")
				{
					slice_plane[n] = 0;
					planes[n].nvec = {0.0, 0.0, 1.0};
					planes[n].tvec = {1.0, 0.0, 0.0};
					planes[n].mvec = utils::cross_prod(planes[n].nvec, planes[n].tvec);
				}
				else if (slice_plane_in[n] == "XZ" || slice_plane_in[n] == "ZX")
				{
					slice_plane[n] = 1;
					planes[n].nvec = {0.0, 1.0, 0.0};
					planes[n].tvec = {1.0, 0.0, 0.0};
					planes[n].mvec = utils::cross_prod(planes[n].nvec, planes[n].tvec);
				}
				else if (slice_plane_in[n] == "YZ" || slice_plane_in[n] == "ZY")
				{
					slice_plane[n] = 2;
					planes[n].nvec = {1.0, 0.0, 0.0};
					planes[n].tvec = {0.0, 1.0, 0.0};
					planes[n].mvec = utils::cross_prod(planes[n].nvec, planes[n].tvec);
				}
				else
				{
					slice_plane[n] = 1e9;
					print("Invalid slice plane selected!");
				}

				// Set plane point
				planes[n].pt = get_plane_pt(n);
				
				print("   Slice = ",n," plane,ID = ",slice_plane[n],slice_plane_in[n]);
				print("      nvec = ",planes[n].nvec[0], planes[n].nvec[1], planes[n].nvec[2]);
				print("      tvec = ",planes[n].tvec[0], planes[n].tvec[1], planes[n].tvec[2]);
				print("      mvec = ",planes[n].mvec[0], planes[n].mvec[1], planes[n].mvec[2]);
				print("      pt   = ",planes[n].pt[0], planes[n].pt[1], planes[n].pt[2]);
			}
		}

		// Distance to plane
		inline value_type dist_to_plane(const auto& islice, const vec_t& xyz)
		{
			vec_t vec;
			vec[0] = xyz[0] - planes[islice].pt[0];
			vec[1] = xyz[1] - planes[islice].pt[1];
			vec[2] = xyz[2] - planes[islice].pt[2];

			return utils::dot_prod(vec, planes[islice].nvec);
		}

		// Distance to plane
		inline vec_t project_onto_plane(const auto& islice, const vec_t& xyz)
		{
			vec_t vec, projPt;
			vec[0] = xyz[0] - planes[islice].pt[0];
			vec[1] = xyz[1] - planes[islice].pt[1];
			vec[2] = xyz[2] - planes[islice].pt[2];
			const auto dist = utils::dot_prod(vec, planes[islice].nvec);

			// Projection
			projPt[0] = xyz[0] - dist * planes[islice].nvec[0];
			projPt[1] = xyz[1] - dist * planes[islice].nvec[1];
			projPt[2] = xyz[2] - dist * planes[islice].nvec[2];

			return projPt;
		}

		// Function to extract data from specified geometry
		template <typename geom_t>
		inline void extract_data(const geom_t& surface)
		{
			print("Running slice plane extractor...");
			
			// First, find average dx on surface geometry
			value_type dx_avg = 0.0;
			for (int n = 0; n<surface.connect.size(); n++)
			{
				// Node list for this element
				const auto& connect = surface.connect[n];

				// Node coordinates
				const auto& node0 = surface.nodes[connect[0]];
				const auto& node1 = surface.nodes[connect[1]];
				const auto& node2 = surface.nodes[connect[2]];

				// Grid spacings
				const auto dx0   = utils::array_norm(node0, node1);
				const auto dx1   = utils::array_norm(node0, node2);
				const auto dx2   = utils::array_norm(node1, node2);

				// Average spacing
				dx_avg += (dx0 + dx1 + dx2) / value_type(3.0 * surface.connect.size());
			}

			print("Avg. surface spacing = ",dx_avg);

			// Search radius in plane vicinity --> based on cell diagonal 
			const value_type search_radius = sqrt(3*std::pow(dx_avg,2));
			print("Search radius        = ",search_radius);
			
			// Loop slice planes
			for (int islice = 0; islice<nslice(); islice++)
			{
				print("   Extracting data on slice plane ",islice);

				// Sweep all elements on provided surface
				std::vector<std::size_t> elemRange; // List of elements close to plane
				std::vector<std::size_t> elemIntersect; // List of intersecting elements
				for (int n = 0; n<surface.connect.size(); n++)
				{
					// Node list for this element
					const auto& list = surface.connect[n];
					
					// Sweep nodes
					int count_pos = 0, count_neg = 0;
					for (int inode = 0; inode<3; inode++)
					{
						// Node coordinate
						const auto& node = surface.nodes[list[inode]];
						
						// Distance from point to plane
						const value_type dist = dist_to_plane(islice, node);
						
						if (dist>0) count_pos++;
						if (dist<0) count_neg++;
					}
					
					// Did we intersect the plane?
					if (count_pos<3 && count_neg<3)	elemIntersect.push_back(n);

					// Element centroid
					const auto centroid = surface.centroid(n);

					// Distance from centroid to plane
					const value_type dist = fabs(dist_to_plane(islice, centroid));
					if (dist<search_radius) elemRange.push_back(n);					
				}

				// Select slice plane data storage
				auto& slice = slice_data[islice];

				// Set cross-section memory
				slice.nodes.resize(elemIntersect.size());
				slice.connect.resize(elemIntersect.size());
				slice.varnames.resize(surface.varnames.size());
				slice.varnames = surface.varnames;
				slice.data.resize(surface.varnames.size());
				for (int i = 0; i<slice.varnames.size(); i++)
				{
					auto& data_loc = slice.data[i];
					data_loc.resize(elemIntersect.size());
				}
				
				// First, sweep points intersecting the plane and build the plane cross-section
				value_type minL =   1e9;
				value_type maxL = - 1e9;
				for (int n = 0; n<elemIntersect.size(); n++)
				{
					// Element ID
					const std::size_t fglob = elemIntersect[n];
					
					// Element centroid
					const auto centroid = surface.centroid(fglob);
					
					// Project point onto plane and store on the slice geometry
					slice.nodes[n] = project_onto_plane(islice, centroid);

					// Copy data onto slice plane
					for (int ivar = 0; ivar<slice.varnames.size(); ivar++) slice.data[ivar][n] = surface.data[ivar][fglob];

					// Find maximum coordinate extent in normalization direction
					if (slice.nodes[n][normalize_dir]<minL) minL = slice.nodes[n][normalize_dir];
					if (slice.nodes[n][normalize_dir]>maxL) maxL = slice.nodes[n][normalize_dir];
				}

				// If normalizing, adjust coordinates now
				if (normalize_coord)
				{
					// Tangential directions to normalization
					int tdir1  = normalize_dir + 1;
					if (tdir1>=3) tdir1 -= 3;
					int tdir2 = normalize_dir + 2;
					if (tdir2>=3) tdir2 -=3;
					// Loop over nodes again
					for (int n = 0; n<slice.nodes.size(); n++)
					{
						// Correct coordinates
						slice.nodes[n][normalize_dir] = (slice.nodes[n][normalize_dir] - minL) / (maxL - minL);
						slice.nodes[n][tdir1]         = slice.nodes[n][tdir1] / (maxL - minL);
						slice.nodes[n][tdir2]         = slice.nodes[n][tdir2] / (maxL - minL);
					}
				}
				
			} // End slice loop
		} // End extraction function

		inline void export_data(const std::string& filename)
		{
			print("Exporting data to tecplot file...");

			// Open file
			std::ofstream fh(filename);

			// Slice number
			const std::string slice_name = "SLICE EXTRACTION";
				
			// Write file title
			fh << "TITLE = \"" << slice_name << "\" \n";

			// Build next header line
			std::string var_header = "VARIABLES = \"x\",\"y\",\"z\",";
			if (normalize_coord) var_header = "VARIABLES = \"x/L\",\"y/L\",\"z/L\",";
			for (int n = 0; n<slice_data[0].varnames.size(); n++)
			{
				if (n<slice_data[0].varnames.size()-1)
				{
					var_header += "\""+slice_data[0].varnames[n]+"\",";
				}
				else
				{
					var_header += "\""+slice_data[0].varnames[n]+"\"";
				}
			}
			fh << var_header << "\n";
			
			for (int islice = 0; islice<nslice(); islice++)
			{
				print("   Exporting slice = ",islice);

				// Get slice data
				const auto& slice = slice_data[islice];

				// Last header line
				const std::string slice_no = "SLICE EXTRACTION "+std::to_string(islice);
				fh << "ZONE T = \"" << slice_no << "\", I = " << slice.nodes.size() << "\n";

				// Sweep data points
				for (int n = 0; n<slice.nodes.size(); n++)
				{
					// Write node coordinates to file
					const auto& node = slice.nodes[n];
					for (int inode = 0; inode<node.size(); inode++) fh << node[inode] << " ";

					// Write data to file
					for (int ivar = 0; ivar<slice.varnames.size(); ivar++) fh << slice.data[ivar][n] << " ";

					// To next line
					fh << "\n";
				}				
			} // End slice loop

			// Close file
			fh.close();
			
		} // End export
	};
}
