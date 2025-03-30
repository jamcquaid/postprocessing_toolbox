#pragma once

#include "spade.h"

namespace pptb::geom
{
	// Generalized geometry structure
	template <typename rtype>
	struct surf_geom_t
	{
		using value_type = rtype;
		using vec_t = spade::ctrs::array<value_type, 3>;

		int maxNpe = 3;
		std::vector<vec_t> nodes;
		std::vector<std::array<std::size_t, 3>> connect;
		std::vector<int> components;
		std::vector<std::string> varnames;
		std::vector<std::vector<value_type>> data;

		// Compute centroid
		inline vec_t centroid(const auto& n) const
		{
			const auto& list = connect[n];
			vec_t xyz{0.0, 0.0, 0.0};
			for (int i = 0; i<3; i++)
			{
				const auto node = nodes[list[i]];
				xyz[0] += node[0] / value_type(3.0);
				xyz[1] += node[1] / value_type(3.0);
				xyz[2] += node[2] / value_type(3.0);
			}
			
			return xyz;
		}
		
		vec_t compute_surface_normal(const auto& idx) const
		{
			const auto& face_list = connect[idx];
			vec_t pp0 = nodes[face_list[0]];
			vec_t pp1 = nodes[face_list[1]];
			vec_t pp2 = nodes[face_list[2]];
			vec_t n = spade::ctrs::cross_prod(pp1-pp0, pp2-pp1);
			n /= spade::ctrs::array_norm(n);
			return n;
		}

		value_type area(const auto& idx) const
		{
			value_type area = 0.0;
			const auto& face = connect[idx];
            const auto& p0   = nodes[face[0]];
            const auto& p1   = nodes[face[1]];
            const auto& p2   = nodes[face[2]];

			// Get cross-products
			const auto tvec1 = p1 - p0;
			const auto tvec2 = p2 - p0;

			// Compute face area
			const auto Avec = spade::ctrs::cross_prod(tvec1, tvec2);
			
            #pragma unroll
            for (int d = 0; d < 3; ++d)
            {
                area += Avec[d] * Avec[d];
            }
			return value_type(0.5) * sqrt(area);
		}
	};
}
