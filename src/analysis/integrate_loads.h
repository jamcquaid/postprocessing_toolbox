#pragma once

namespace pptb::analysis
{
	template <typename real_t>
	static void integrate_loads(const geom::surf_geom_t<real_t>& geom)
	{
		real_t fz = 0.0;
		for (int i = 0; i<geom.connect.size(); ++i)
		{
			fz -= geom.data[1][i] * geom.area(i) * geom.compute_surface_normal(i)[2];
		}
		print(fz);
	}
}
