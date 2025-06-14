#pragma once

namespace pptb::analysis
{
	template <std::floating_point rtype>
	struct spanwiseAverage_t
	{
		using value_type = rtype;

		// Geometry structure
		geom::surf_geom_t<value_type> geom;

		const std::string out_file;

	    spanwiseAverage_t(const geom::surf_geom_t<value_type>& geom3D, const std::string& out_file_in, const int& idir) : out_file{out_file_in}
		{
			print("Spanwise averaging...");

			// First, slice the 3D geometry we were provided to get sampling points on 2D plane
			print("Slicing 3D structure...");
			std::vector<std::string> slice_name{"spanwiseAverage"};
			std::vector<std::string> slice_plane;
			std::vector<std::string> slice_file{""};
			if (idir == 0) slice_plane.push_back("YZ");
			if (idir == 1) slice_plane.push_back("XZ");
			if (idir == 2) slice_plane.push_back("XY");
			const std::vector<value_type> slice_pos{0.0};
			const bool normalize_coord    = false;
			const bool normalize_dir      = 0;
			std::vector<truncate_t<value_type>> truncate;
			plane_slice_t<value_type> plane_slice(slice_name, slice_plane, slice_pos, slice_file, truncate, value_type(1.0), normalize_coord, normalize_dir);

			// Slice data
			plane_slice.extract_data(geom3D);

			// Setup spanwise averaged geometry structure
			const auto& slice = plane_slice.slice_data[0];
			const int npts    = slice.connect.size();
			const int nvars   = slice.varnames.size();
			
			// Tangential directions to normalization
			int tdir1  = idir + 1;
			if (tdir1>=3) tdir1 -= 3;
			int tdir2 = idir + 2;
			if (tdir2>=3) tdir2 -= 3;

			// Sweep 3D structure
			print("   Mapping time averaged data onto 2D slice plane...");
			std::vector<int> counter;
			counter.resize(npts, 0);
			std::vector<int> mapping;
			mapping.resize(geom3D.connect.size(), -10);
			for (int n = 0; n<geom3D.connect.size(); ++n)
			{
				// Only valid elements
				if (geom3D.data[0][n]>0)
				{
					// Element centroid
					value_type centroid[3] = {0.0, 0.0, 0.0};
					const auto& list = geom3D.connect[n];
					for (int j = 0; j<geom3D.maxNpe; ++j)
					{
						centroid[0] += geom3D.nodes[list[j]][0] / value_type(geom3D.maxNpe);
						centroid[1] += geom3D.nodes[list[j]][1] / value_type(geom3D.maxNpe);
						centroid[2] += geom3D.nodes[list[j]][2] / value_type(geom3D.maxNpe);
					}

					// Gonna brute force this part because I'm lazy
					value_type dist_min = 1e9;
					int minIdx          = -10;
					for (int i = 0; i<npts; ++i)
					{
						// Find closest point on 2D geom
						value_type dist = sqrt(pow(centroid[tdir1] - slice.nodes[i][tdir1],2) + pow(centroid[tdir2] - slice.nodes[i][tdir2],2));
						
						if (dist<dist_min)
						{
							dist_min = dist;
							minIdx = i;
						}
					}

					// Run counter
					++counter[minIdx];

					// Mapping
					mapping[n] = minIdx;
				}
			}

			int npts2 = 0;
			for (int j = 0; j<npts; ++j)
				if (counter[j]>0)
					++npts2;

			print("   Initializing spanwise average data structure...");
			geom.maxNpe = 2;
			geom.nodes.resize(npts2);
			geom.connect.resize(npts2);
			geom.varnames.resize(nvars);
			geom.varnames = slice.varnames;
			geom.data.resize(nvars);
			for (int i = 0; i<nvars; i++)
			{
				auto& data_loc = geom.data[i];
				data_loc.resize(npts2);
				for (int j = 0; j<npts2; j++) data_loc[j] = 0.0;
			}

			// Set nodes and connectivity
			print("   Setting nodes and connectivity...");
			int count = 0;
			std::vector<int> mapping2;
			mapping2.resize(npts, -10);
			for (int n = 0; n<npts; ++n)
			{
				if (counter[n]>0)
				{
					// Set coordinates
					geom.nodes[count][0] = slice.nodes[n][0];
					geom.nodes[count][1] = slice.nodes[n][1];
					geom.nodes[count][2] = slice.nodes[n][2];

					// Set connectivity
					auto& list = geom.connect[count];
					list[0] = count;
					list[1] = count+1;
					if (list[1]>npts2-1) list[1] = 0;

					// Fix mapping
					mapping2[n] = count;
					
					// Counter
					++count;
				}
			}
			
			// Final averaging --> data
			print("   Averaging data onto plane...");
			for (int ivar = 0; ivar<nvars; ++ivar)
			{
				auto& data_loc = geom.data[ivar];

				// Sweep 3D geom
				for (int n = 0; n<geom3D.connect.size(); ++n)
				{
					// Only valid elements
					if (mapping[n]>=0)
					{
						const int idx = mapping2[mapping[n]];
						data_loc[idx] += geom3D.data[ivar][n] / value_type(counter[mapping[n]]);
					}
				}
			}
			print("Complete!");
		}

		// Write output VTK
		inline void export_data()
		{
			print("Exporting data to tecplot file...");

			// Open file
			std::string delimiter = ".";
			std::string token = out_file.substr(0, out_file.find(delimiter));
			std::string out_file2 = token+".dat";
			std::ofstream fh(out_file2);

			const std::string title = "SPANWISE AVERAGE";
				
			// Write file title
			fh << "TITLE = \"" << title << "\" \n";

			// Build next header line
			std::string var_header = "VARIABLES = \"x\",\"y\",\"z\",";
			for (int n = 0; n<geom.varnames.size(); n++)
			{
				if (n<geom.varnames.size()-1)
				{
					var_header += "\""+geom.varnames[n]+"\",";
				}
				else
				{
					var_header += "\""+geom.varnames[n]+"\"";
				}
			}
			fh << var_header << "\n";
			

			// Last header line
			const std::string zone_name = "SPANWISE AVERAGE";
			fh << "ZONE T = \"" << zone_name << "\", I = " << geom.nodes.size() << "\n";

			// Sweep data points
			for (int n = 0; n<geom.nodes.size(); n++)
			{
				// Write node coordinates to file
				const auto& node = geom.nodes[n];
				for (int inode = 0; inode<node.size(); inode++) fh << node[inode] << " ";
				
				// Write data to file
				for (int ivar = 0; ivar<geom.varnames.size(); ivar++) fh << geom.data[ivar][n] << " ";
				
				// To next line
				fh << "\n";
			}				
			
			// Close file
			fh.close();

			// Also write a VTK file just because
			io::write_vtk_data<value_type>(out_file, geom);
			
		} // End export
	};
}
