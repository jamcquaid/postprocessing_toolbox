namespace pptb::io
{
	// Generalized geometry structure
	template <typename rtype>
	struct surf_geom_t
	{
		using value_type = rtype;

		std::vector<std::array<rtype, 3>> nodes;
		std::vector<std::array<std::size_t, 3>> connect;
		std::vector<int> components;
	};
}
