#pragma once

namespace pptb::analysis
{
	template <std::floating_point rtype>
	struct ideal_gas_t
	{
		using real_t = rtype;

		// Gas parameters
		const real_t Rgas = 287.15;
		const real_t gamma= 1.4;
		const real_t Pr   = 0.72;
		const real_t Cp   = 1005.0;
		const real_t Cv   = 717.85;
		const real_t c2   = 110.4;

		real_t get_mu(const real_t& T) const
		{
			return real_t(1.45151376745308e-06) * pow(T, real_t(1.5)) / (c2 + T);
		}
	};
}
