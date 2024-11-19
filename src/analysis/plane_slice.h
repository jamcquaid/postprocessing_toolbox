#include <algorithm>

namespace pptb::analysis
{
	template<std::floating_point rtype>
	struct plane_t
	{
		using value_type = rtype;
		
		// Member variables
		std::array<value_type, 3> nvec; // Plane normal
		std::array<value_type, 3> tvec; // Plane tangential vector 1
		std::array<value_type, 3> mvec; // plane tangential vector 2
		std::array<value_type, 3> pt;
	};

	template<std::floating_point rtype>
	struct sliceElement_t
	{
		rtype xyz[3];
		rtype data[100];
	};
	
	template <std::floating_point rtype>
	struct plane_slice_t
	{
		using value_type  = rtype;
		using surf_geom_t = geom::surf_geom_t<value_type>;
		using vec_t = std::array<value_type, 3>;
		using container_t = std::vector<vec_t>;
		
		// Member variables
		std::vector<std::string> slice_name;
		std::vector<surf_geom_t> slice_data; // Cross-section geometry and data
		std::vector<std::size_t> slice_plane; // Slice plane
		std::vector<value_type> slice_pos; // Slice position in plane normal direction
		std::vector<std::string> slice_file; // Import slice plane points from file
		std::vector<container_t> slice_file_points; // Imported point coordinates from slice file
		std::vector<plane_t<value_type>> planes; // Slice plane basis
		bool normalize_coord = false;
		int normalize_dir    = 0;
		int nslice;

		// Get coordinate on plane
		inline vec_t get_plane_pt(const auto& iplane)
		{
			vec_t pt{0.0, 0.0, 0.0};
			pt[slice_plane[iplane]] = slice_pos[iplane];
			return pt;
		}

		void write_slice_details(const int& n, const auto& slice_plane_in)
		{
			if (slice_plane[n] < 3)
			{
				print("   Slice = ",n," name = ",slice_name[n], " plane,ID = ",slice_plane[n],slice_plane_in[n]);
				print("      nvec = ",planes[n].nvec[0], planes[n].nvec[1], planes[n].nvec[2]);
				print("      tvec = ",planes[n].tvec[0], planes[n].tvec[1], planes[n].tvec[2]);
				print("      mvec = ",planes[n].mvec[0], planes[n].mvec[1], planes[n].mvec[2]);
				print("      pt   = ",planes[n].pt[0], planes[n].pt[1], planes[n].pt[2]);
			}
			else
			{
				print("   Slice = ",n," name = ",slice_name[n], " plane,ID = ",slice_plane[n],slice_plane_in[n]);
				print("      file, npts = ",slice_file[n],slice_file_points[n].size());
			}
		}

		void import_slice_points(const int& n)
		{
			auto& slice_points = slice_file_points[n];

			// First check for number of points in file
			int count = 0;
			utils::ascii_file_t fh(slice_file[n]);
			while (!fh.eof())
			{
				fh.next_line();
				if (fh.eof()) break;
				++count;
			}

			// Resize container
			slice_points.resize(count);

			// Close file
			fh.fh.close();

			// Reopen file
			utils::ascii_file_t file(slice_file[n]);

			// Import data
			for (int i = 0; i<slice_points.size(); ++i)
			{
				auto& ptArr = slice_points[i];
				file.next_line();
				file.parse(ptArr[0], ptArr[1], ptArr[2]);
			}
		}
		
		// Constructor --> slice select based on Cartesian planes
	    plane_slice_t(const auto& slice_name_in, const auto& slice_plane_in, const auto& slice_pos_in,
					  const auto& slice_file_in, const bool& norm_coord_in = false, const int& norm_dir_in = 0)
		: slice_pos{slice_pos_in}, slice_file{slice_file_in}, normalize_coord{norm_coord_in}, normalize_dir{norm_dir_in}
		{
			print("Initializing slice operator...");

			// Number of plane slices
			nslice = slice_pos.size();

			// Allocate memory
			slice_name.resize(nslice);
			slice_data.resize(nslice);
			slice_plane.resize(nslice);
			slice_pos.resize(nslice);
			slice_file_points.resize(nslice);
			planes.resize(nslice);
			
			for (int n = 0; n<nslice; n++)
			{
				slice_name[n] = slice_name_in[n];
				if (slice_plane_in[n] == "XY" || slice_plane_in[n] == "YX")
				{
					slice_plane[n] = 2;
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
					slice_plane[n] = 0;
					planes[n].nvec = {1.0, 0.0, 0.0};
					planes[n].tvec = {0.0, 1.0, 0.0};
					planes[n].mvec = utils::cross_prod(planes[n].nvec, planes[n].tvec);
				}
				else if (slice_plane_in[n] == "file")
				{
					// Do something here
					slice_plane[n] = 3;
					import_slice_points(n);
				}
				else
				{
					slice_plane[n] = 1e9;
					print("Invalid slice plane selected!");
				}

				// Set plane point (planes specified by Cartesian directions)
				if (slice_plane[n] < 3) planes[n].pt = get_plane_pt(n);

				// Write details to screen
				write_slice_details(n, slice_plane_in);
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

		// Slice provided data set using specified plane vectors
		template <typename geom_t>
		surf_geom_t slice_by_plane(const int& islice, const int& tdir1, const int& tdir2, const geom_t& surface)
		{
			// Output type
			surf_geom_t slice;

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
			}

			// Select slice plane data storage
			slice.maxNpe = 2; // Planar vtk

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
				int ntdir1  = normalize_dir + 1;
				if (ntdir1>=3) ntdir1 -= 3;
				int ntdir2 = normalize_dir + 2;
				if (ntdir2>=3) ntdir2 -= 3;
					
				// Loop over nodes again
				for (int n = 0; n<slice.nodes.size(); n++)
				{
					// Correct coordinates
					slice.nodes[n][normalize_dir] = (slice.nodes[n][normalize_dir] - minL) / (maxL - minL);
					slice.nodes[n][ntdir1]         = slice.nodes[n][ntdir1] / (maxL - minL);
					slice.nodes[n][ntdir2]         = slice.nodes[n][ntdir2] / (maxL - minL);
				}
			}

			return slice;
		}

		// Extract data using preset coordinates
		template <typename geom_t>
		surf_geom_t slice_by_file(const int& islice, const geom_t& surface)
		{
			// Get slice file points
			const auto& file_points = slice_file_points[islice];
			
			// Initialize slice operator
			surf_geom_t slice;

			// Select slice plane data storage
			slice.maxNpe = 2; // Planar vtk

			// Set cross-section memory
			const int npts = slice_file_points[islice].size();
			slice.nodes.resize(npts);
			slice.connect.resize(npts);
			slice.varnames.resize(surface.varnames.size());
			slice.varnames = surface.varnames;
			slice.data.resize(surface.varnames.size());
			for (int i = 0; i<slice.varnames.size(); i++)
			{
				auto& data_loc = slice.data[i];
				data_loc.resize(npts);
			}

			// Going to brute force this for now
			for (int ipt = 0; ipt<npts; ++ipt)
			{
				// Coordinate on slice plane
				const auto& xyz = file_points[ipt];

				// Store in slice struct
				slice.nodes[ipt][0] = xyz[0];
				slice.nodes[ipt][1] = xyz[1];
				slice.nodes[ipt][2] = xyz[2];
				
				int minIdx = 0;
				value_type minDist = 1e9;
				for (int n = 0; n<surface.connect.size(); ++n) // Sweep centroids
				{
					if (surface.data[0][n] > 0)
					{
						// Element centroid
						const auto centroid = surface.centroid(n);

						// Compute distance
						const value_type dist = utils::array_norm(xyz, centroid);

						if (dist < minDist)
						{
							minDist = dist;
							minIdx  = n;
						}
					}
				}

				// Assign data from identified element onto slice plane
				for (int n = 0; n<slice.varnames.size(); n++)
				{
					slice.data[n][ipt] = surface.data[n][minIdx];
				}
			}

			// Get min/max L
			value_type minL =   1e9;
			value_type maxL = - 1e9;
			for (int n = 0; n<npts; n++)
			{
				// Find maximum coordinate extent in normalization direction
				if (slice.nodes[n][normalize_dir]<minL) minL = slice.nodes[n][normalize_dir];
				if (slice.nodes[n][normalize_dir]>maxL) maxL = slice.nodes[n][normalize_dir];
			}

			// If normalizing, adjust coordinates now
			if (normalize_coord)
			{
				// Tangential directions to normalization
				int ntdir1  = normalize_dir + 1;
				if (ntdir1>=3) ntdir1 -= 3;
				int ntdir2 = normalize_dir + 2;
				if (ntdir2>=3) ntdir2 -= 3;
					
				// Loop over nodes again
				for (int n = 0; n<slice.nodes.size(); n++)
				{
					// Correct coordinates
					slice.nodes[n][normalize_dir] = (slice.nodes[n][normalize_dir] - minL) / (maxL - minL);
					slice.nodes[n][ntdir1]         = slice.nodes[n][ntdir1] / (maxL - minL);
					slice.nodes[n][ntdir2]         = slice.nodes[n][ntdir2] / (maxL - minL);
				}
			}
			
			return slice;
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
			for (int islice = 0; islice<nslice; islice++)
			{
				print("   Extracting data on slice plane ",islice);
				
				const int tdir1 = std::abs(0*planes[islice].tvec[0] + 1*planes[islice].tvec[1] + 2*planes[islice].tvec[2]);
				const int tdir2 = std::abs(0*planes[islice].mvec[0] + 1*planes[islice].mvec[1] + 2*planes[islice].mvec[2]);
				
				auto& slice = slice_data[islice];
				if (slice_plane[islice]<3)
				{
					slice = slice_by_plane(islice, tdir1, tdir2, surface);
				}
				else
				{
					// Slice by points in file
					slice = slice_by_file(islice, surface);
				}

				// Sorting only works when aligned with Cartesian directions for now...
//				if (slice_plane[islice] < 3)
//				{
//					// Get min/max L and geom centroid
//					const int npts = slice.nodes.size();
//					std::vector<value_type> geom_centroid;
//					geom_centroid.resize(3, 0.0);
//					for (int n = 0; n<npts; n++)
//					{
//						// Geometry midpoint
//						geom_centroid[0] += slice.nodes[n][0] / value_type(npts);
//						geom_centroid[1] += slice.nodes[n][1] / value_type(npts);
//						geom_centroid[2] += slice.nodes[n][2] / value_type(npts);
//					}
//					
//					// Assemble structure to be sorted
//					std::vector<sliceElement_t<value_type>> sliceElements;
//					sliceElements.resize(slice.nodes.size());
//					for (int n = 0; n<slice.nodes.size(); n++)
//					{
//						sliceElements[n].xyz[0] = slice.nodes[n][0];
//						sliceElements[n].xyz[1] = slice.nodes[n][1];
//						sliceElements[n].xyz[2] = slice.nodes[n][2];
//						for (int ivar = 0; ivar<slice.varnames.size(); ivar++)
//							sliceElements[n].data[ivar] = slice.data[ivar][n];
//					}
//
//					// Comparison function
//					auto compare = [&](const sliceElement_t<value_type>& elemA, const sliceElement_t<value_type>& elemB)
//					{
//						const value_type phiA = atan2(elemA.xyz[tdir2] - geom_centroid[tdir2], elemA.xyz[tdir1] - geom_centroid[tdir1]);
//						const value_type phiB = atan2(elemB.xyz[tdir2] - geom_centroid[tdir2], elemB.xyz[tdir1] - geom_centroid[tdir1]);
//						return phiA > phiB;
//					};
//				
//					// Sort points
//					std::sort(sliceElements.begin(), sliceElements.end(), compare);
//
//					// Push data back into slice geom
//					for (int n = 0; n<npts; n++)
//					{
//						slice.nodes[n][0] = sliceElements[n].xyz[0];
//						slice.nodes[n][1] = sliceElements[n].xyz[1];
//						slice.nodes[n][2] = sliceElements[n].xyz[2];
//						for (int ivar = 0; ivar<slice.varnames.size(); ivar++)
//							slice.data[ivar][n] = sliceElements[n].data[ivar];
//					}
//				}
				
				// Set connectivity
				for (int n = 0; n<slice.nodes.size(); n++)
				{
					auto& list = slice.connect[n];
					list[0] = n;
					list[1] = n+1;
					if (list[1]>slice.nodes.size()-1) list[1] = 0;
				}
				
			} // End slice loop
		} // End extraction function

		inline void export_data(const std::string& filename)
		{
			print("Exporting data to tecplot file...");

			// Open file
			std::ofstream fh(filename);

			// Slice number
			const std::string slice_tmp = "SLICE EXTRACTION";
				
			// Write file title
			fh << "TITLE = \"" << slice_tmp << "\" \n";

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
			
			for (int islice = 0; islice<nslice; islice++)
			{
				print("   Exporting slice = ",islice);

				// Get slice data
				const auto& slice = slice_data[islice];

				// Last header line
				const std::string slice_no = slice_name[islice];
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
