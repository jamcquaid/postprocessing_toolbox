namespace pptb::geom
{
	// Generalized geometry structure
	template <typename rtype>
	struct surf_geom_t
	{
		using value_type = rtype;
		using vec_t = std::array<value_type, 3>;

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
	};
}
