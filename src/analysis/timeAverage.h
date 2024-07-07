namespace pptb::analysis
{
	template <typename rtype>
	struct timeAverage_t
	{
		using value_type = rtype;
		
		// Member variables
		std::size_t min_nt;
		std::size_t max_nt;
		std::size_t max_vars;

		// Geometry structure to store data
		geom::surf_geom_t<value_type> geom;

		// Constructor
    	timeAverage_t(const auto& minNt_in, const auto& maxNt_in, const auto& maxVars_in, const auto& fileStruct) : min_nt{minNt_in}, max_nt{maxNt_in}, max_vars{maxVars_in}
		{
			print("Beginning time average over range min, max = ",min_nt, max_nt);

			// Import first file to set baseline
			print("Importing first file to initialize...");
			io::import_vtk(fileStruct.filenames[0], geom);

			// Run loop over remaining files
			value_type running_time = 1.0;
			for (int n = 1; n<fileStruct.nFiles; n++)
			{
				print("   Importing file = ",n+1);
				
				// Temporary structure
				geom::surf_geom_t<value_type> tmp_geom;

				// Import this file
				io::import_vtk(fileStruct.filenames[n], tmp_geom, false);

				// Averaging coefficients
				const value_type alpha = value_type(1.0) / (running_time + value_type(1.0));
				const value_type beta  = value_type(1.0) - alpha;
				running_time += value_type(1.0);
				
				// Run time-average into main structure
				print("   Averaging data...");
				for (int ivar = 0; ivar<std::min(geom.varnames.size(), max_vars); ivar++)
				{
					for (int f = 0; f<geom.connect.size(); f++)
					{
						geom.data[ivar][f] *= beta;
						geom.data[ivar][f] += alpha*tmp_geom.data[ivar][f];
					}
				}
			}	
		}

		// Write output VTK
		void export_data()
		{
			io::write_vtk_data<value_type>("timeAverage.vtk", geom, max_vars);
		}
		
	};
}