#pragma once

namespace pptb::analysis
{
	template <typename real_t, typename kernel_t>
	void modify_scalars(geom::surf_geom_t<real_t>& surface, const std::string& varname, const kernel_t& kernel)
	{
		// Loop stored scalars
		for (int v = 0; v<surface.varnames.size(); ++v)
		{
			// Varname
			const std::string var_tmp = surface.varnames[v];
			
			// Is this the one we want?
			if (var_tmp == varname)
			{
				// Get data for this variable
				auto& data_loc = surface.data[v];
				
				// Sweep surface grid
				for (int n = 0; n<surface.connect.size(); n++)
				{
					// Get cell centroid
					const auto centroid = surface.centroid(n);

					// Modify data
					data_loc[n] = kernel(centroid, data_loc[n]);
				}
			}
		}
	}
}
