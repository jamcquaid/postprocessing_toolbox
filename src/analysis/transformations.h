#pragma once

namespace pptb::analysis
{
	enum transform_type
	{
		rotate,
		translate,
		scale,
		merge
	};

	template <const transform_type ttype, typename real_t>
	requires(ttype == rotate)
	static void transform_vtk(geom::surf_geom_t<real_t>& geom, const int& rotate_dir, const real_t& rotate_angle)
	{
		utils::mat_t<real_t, 3> rot;

		print("Building rotation matrix...");
		
		// Angle in radians
		const real_t pi    = 3.1415926535;
		const real_t theta = rotate_angle * pi / 180.0;
		
		for (int i = 0; i<rot.num_rows(); ++i)
		{
			for (int j = 0; j<rot.num_cols(); ++j)
			{
				if (rotate_dir == 2)
				{
					// Z-axis rotation
					rot(0,0) =   cos(theta);
					rot(0,1) = - sin(theta);
					rot(0,2) =   0.0;
					rot(1,0) =   sin(theta);
					rot(1,1) =   cos(theta);
					rot(1,2) =   0.0;
					rot(2,0) =   0.0;
					rot(2,1) =   0.0;
					rot(2,2) =   1.0;
				}
				else if (rotate_dir == 1)
				{
					// Y-axis rotation
					rot(0,0) =   cos(theta);
					rot(0,1) =   0.0;
					rot(0,2) =   sin(theta);
					rot(1,0) =   0.0;
					rot(1,1) =   1.0;
					rot(1,2) =   0.0;
					rot(2,0) = - sin(theta);
					rot(2,1) =   0.0;
					rot(2,2) =   cos(theta);
				}
				else
				{
					// X-axis rotation
					rot(0,0) =   1.0;
					rot(0,1) =   0.0;
					rot(0,2) =   0.0;
					rot(1,0) =   0.0;
					rot(1,1) =   cos(theta);
					rot(1,2) = - sin(theta);
					rot(2,0) =   0.0;
					rot(2,1) =   sin(theta);
					rot(2,2) =   cos(theta);
				}
			}
		}

		print("Rotating geometry...");
		for (int n = 0; n< geom.nodes.size(); ++n)
		{
			// Get current node coordinates
			auto  xyz_old = geom.nodes[n];
			auto& xyz_new = geom.nodes[n];
			xyz_new[0] = 0.0;
			xyz_new[1] = 0.0;
			xyz_new[2] = 0.0;

			// Rotate coordinate
			for (int i = 0; i<rot.num_rows(); ++i)
			{
				for (int j = 0; j<rot.num_cols(); ++j)
				{
					xyz_new[i] += rot(i,j) * xyz_old[j];
				}
			}
		}
	}

	template <const transform_type ttype, typename real_t>
	requires(ttype == merge)
	static geom::surf_geom_t<real_t> transform_vtk(const utils::file_directory_t& fileStruct)
	{
		print("Merging ",fileStruct.nFiles," vtk files...");

		// Final geometry
		geom::surf_geom_t<real_t> geom;

		int nNodes = 0;
		int nElem  = 0;
		int nVars  = 0;
		
		// Loop files
		for (int ifile = 0; ifile<fileStruct.nFiles; ++ifile)
		{
			// Import geometry
			using geom_t = pptb::geom::surf_geom_t<real_t>;
			geom_t geom_in;
			pptb::io::import_vtk(fileStruct.filenames[ifile], geom_in);

			nNodes += geom_in.nodes.size();
			nElem  += geom_in.connect.size();
			nVars   = 1;
		}

		geom.nodes.resize(nNodes);
		geom.connect.resize(nElem);
		geom.varnames.resize(nVars);
		geom.data.resize(nVars);
		geom.data[0].resize(nElem);

		print("nNodes, nElem, nVars = ",nNodes,nElem,nVars);
		
		nNodes = 0;
		nElem  = 0;
		int nElem2 = 0;
		
		// Loop files
		int offset = 0;
		for (int ifile = 0; ifile<fileStruct.nFiles; ++ifile)
		{
			// Import geometry
			using geom_t = pptb::geom::surf_geom_t<real_t>;
			geom_t geom_in;
			pptb::io::import_vtk(fileStruct.filenames[ifile], geom_in);
			
			// Add nodes
			print("   Adding nodes...");
			for (int n = 0; n<geom_in.nodes.size(); ++n)
			{
				geom.nodes[nNodes] = geom_in.nodes[n];
				++nNodes;
			}

			if (ifile == 1)
			{
				geom.varnames[0] = geom_in.varnames[0];
			}
			
			// Add connectivity
			print("   Adding connectivity...");
			for (int n = 0; n<geom_in.connect.size(); ++n)
			{
				geom.connect[nElem] = geom_in.connect[n];
				// Correct connectivity with offset
				geom.connect[nElem][0] += offset;
				geom.connect[nElem][1] += offset;
				geom.connect[nElem][2] += offset;
				++nElem;
			}

			// Connectivity offset
			offset += geom_in.nodes.size();

			// Add variable data (assumed 1 variable for components)
			print("   Adding variable data...");
			for (int n = 0; n<geom_in.connect.size(); ++n)
			{
				geom.data[0][nElem2] = geom_in.data[0][n];
				++nElem2;
			}
			print("   done...");
		}
			
		return geom;
	}
}
