#pragma once

#include "spade.h"
#include "gas.h"

namespace pptb::analysis
{
	template <std::floating_point rtype>
	struct solver_inputs_t
	{
		using real_t = rtype;
		real_t Mae;
		real_t Te;
		real_t rhoe;
		real_t Tw = -10;
		real_t xi = 0.25;
	};

	template <typename ctype>
	static std::ostream& operator << (std::ostream& os, const solver_inputs_t<ctype>& data)
	{
		os << " Mae: " << data.Mae << " Te: " << data.Te << " rhoe: " << data.rhoe;
		if (data.Tw<0) os << " Tw: adiabatic";
		else os << " Tw: " << data.Tw;
		os << " xi: " << data.xi;
		
		return os;
	}

	template <std::floating_point rtype>
	struct solver_outputs_t
	{
		using real_t = rtype;
		using vec_t  = spade::ctrs::array<real_t, 5>;
		int npts;
		real_t deta;
		std::vector<vec_t> data;
		std::vector<vec_t> data_dim;
		std::vector<real_t> eta;
		
    	solver_outputs_t(const int& npts_in, const real_t& deta_in) : npts{npts_in}, deta{deta_in}
		{
			data.resize(npts, 0.0);
			data_dim.resize(npts, 0.0);
			eta.resize(npts, 0.0);
			for (int i = 0; i<npts; ++i) eta[i] = i*deta;
		}
	};

	template <std::floating_point rtype>
	struct similarity_solver_t
	{
		using real_t      = rtype;
		using input_type  = solver_inputs_t<real_t>;
		using output_type = solver_outputs_t<real_t>;
		using vec_t       = spade::ctrs::array<real_t, 5>;

		// Inputs
		input_type inputs;

		// Solver settings
		real_t deta      = 0.01;
		int npts         = 3000;
		real_t delta     = 1E-11;
		real_t alpha_ini = 0.1;
		real_t beta_ini  = 3.0;

		// Gas model
		ideal_gas_t<real_t> gas;

		// Convergence settings
		real_t TOL       = 1E-9;
		int maxIter      = 1000;
		
	    similarity_solver_t(const input_type& inputs_in) : inputs{inputs_in} {}

		void set_boundary_condition(output_type& output, const real_t& alpha, const real_t& beta) const
		{
			// Velocity BC
			output.data[0][0] = 0.0;
			output.data[0][1] = 0.0;

			// Temperature BC
			if (inputs.Tw < 0.0)
			{ // Adiabatic
				// Zero gradient at wall
				output.data[0][4] = 0.0;

				// Guessed BC's
				output.data[0][2] = alpha;
				output.data[0][3] = beta;
			}
			else
			{
				// Set wall temperature
				output.data[0][3] = inputs.Tw / inputs.Te;

				// Guessed BC's
				output.data[0][2] = alpha;
				output.data[0][4] = beta;
			}
		}

		real_t comp_eqn1(const real_t& valIn) const
		{
			return valIn;
		}

		real_t comp_eqn2(const real_t& valIn) const
		{
			return valIn;
		}

		real_t comp_eqn3(const real_t& y0, const real_t& y2, const real_t& y3, const real_t& y4) const
		{
			return - y2 * (y4 / (real_t(2.0) * y3) - y4 / (y3 + gas.c2 / inputs.Te)) - y0 * y2 * ((y3 + gas.c2/inputs.Te) / (sqrt(y3) * (real_t(1.0) + gas.c2/inputs.Te)));
		}

		real_t comp_eqn4(const real_t& valIn) const
		{
			return valIn;
		}

		real_t comp_eqn5(const real_t& y0, const real_t& y2, const real_t& y3, const real_t& y4) const
		{
			return - y4*y4 * (real_t(1.0) / (real_t(2.0) * y3) - real_t(1.0) / (y3 + gas.c2/inputs.Te)) - gas.Pr * (y0 * y4 / sqrt(y3)) * (y3 + gas.c2/inputs.Te) /
				(real_t(1.0) + gas.c2/inputs.Te) - (gas.gamma - real_t(1.0)) * gas.Pr * inputs.Mae * inputs.Mae * y2 * y2;
		}
		
		void RK4_solver(output_type& output) const
		{
			vec_t k1, k2, k3, k4;
			for (int i = 0; i<(npts-1); ++i)
			{
				// Get vector for this point
				const auto& y = output.data[i];
				auto& yp1     = output.data[i+1];

				// First substep
				k1[0] = comp_eqn1(y[1]);
				k1[1] = comp_eqn2(y[2]);
				k1[2] = comp_eqn3(y[0], y[2], y[3], y[4]);
				k1[3] = comp_eqn4(y[4]);
				k1[4] = comp_eqn5(y[0], y[2], y[3], y[4]);

				// Second substep
				k2[0] = comp_eqn1(y[1] + real_t(0.5) * deta * k1[1]);
				k2[1] = comp_eqn2(y[2] + real_t(0.5) * deta * k1[2]);
				k2[2] = comp_eqn3(y[0] + real_t(0.5) * deta * k1[0],
								  y[2] + real_t(0.5) * deta * k1[2],
								  y[3] + real_t(0.5) * deta * k1[3],
								  y[4] + real_t(0.5) * deta * k1[4]);
				k2[3] = comp_eqn4(y[4] + real_t(0.5) * deta * k1[4]);
				k2[4] = comp_eqn5(y[0] + real_t(0.5) * deta * k1[0],
								  y[2] + real_t(0.5) * deta * k1[2],
								  y[3] + real_t(0.5) * deta * k1[3],
								  y[4] + real_t(0.5) * deta * k1[4]);

				// Third substep
				k3[0] = comp_eqn1(y[1] + real_t(0.5) * deta * k2[1]);
				k3[1] = comp_eqn2(y[2] + real_t(0.5) * deta * k2[2]);
				k3[2] = comp_eqn3(y[0] + real_t(0.5) * deta * k2[0],
								  y[2] + real_t(0.5) * deta * k2[2],
								  y[3] + real_t(0.5) * deta * k2[3],
								  y[4] + real_t(0.5) * deta * k2[4]);
				k3[3] = comp_eqn4(y[4] + real_t(0.5) * deta * k2[4]);
				k3[4] = comp_eqn5(y[0] + real_t(0.5) * deta * k2[0],
								  y[2] + real_t(0.5) * deta * k2[2],
								  y[3] + real_t(0.5) * deta * k2[3],
								  y[4] + real_t(0.5) * deta * k2[4]);

				// Fourth substep
				k4[0] = comp_eqn1(y[1] + deta * k3[1]);
				k4[1] = comp_eqn2(y[2] + deta * k3[2]);
				k4[2] = comp_eqn3(y[0] + deta * k3[0],
								  y[2] + deta * k3[2],
								  y[3] + deta * k3[3],
								  y[4] + deta * k3[4]);
				k4[3] = comp_eqn4(y[4] + deta * k3[4]);
				k4[4] = comp_eqn5(y[0] + deta * k3[0],
								  y[2] + deta * k3[2],
								  y[3] + deta * k3[3],
								  y[4] + deta * k3[4]);
				
				// Update next grid point
				yp1 = y + real_t(1./6.) * deta * (k1 + real_t(2.0) * k2 + real_t(2.0) * k3 + k4);
			}
		}

		void dimensionalize_profile(output_type& output) const
		{
			print("Dimensionalizing profile...");

			// Edge velocity
			const real_t ae = sqrt(gas.gamma * gas.Rgas * inputs.Te);
			const real_t ue = ae * inputs.Mae;
			const real_t mue= gas.get_mu(inputs.Te);
			const real_t se = inputs.rhoe * ue * mue * inputs.xi;
			
			for (int i = 0; i<npts; ++i)
			{
				const auto& sol = output.data[i];
				auto& data      = output.data_dim[i];

				if (i==0) data[0] = 0.0; // Y-coordinate
				else
				{
					const real_t alpha = sol[3] * inputs.rhoe * ue / sqrt(real_t(2.0) * se);
					data[0]            = output.data_dim[i-1][0] + deta / alpha; // Y-coordinate
				}
				data[1] = inputs.rhoe * (real_t(1.0) / sol[3]); // density profile
				data[2] = ue * sol[1]; // U-velocity
				data[4] = inputs.Te * sol[3]; // Temperature profile

				//const real_t rho     = data[1]; // local density
				//const real_t deta_dx = output.eta[i] / dx;
				//data[3] = - (real_t(1.0) / rho) * ((real_t(1.0) / sqrt(real_t(2.0) * se)) * sol[0] * mue * inputs.rhoe * ue + sqrt(real_t(2.0) * se) * sol[1] * deta_dx); // V-velocity
				data[3] = real_t(0.5) * sqrt(mue * ue / (inputs.rhoe * inputs.xi)) * (output.eta[i] * sol[1] - sol[0]);
			}
		}

		void write_output(const output_type& output) const
		{
			print("Writing output...");
			//const std::string fname = "similarityProfile_Ma"+std::to_string(inputs.Mae)+"_Te"+std::to_string(inputs.Te)+"_rhoe"+std::to_string(inputs.rhoe)+"_Tw"+std::to_string(inputs.Tw)+"_xi"+std::to_string(inputs.xi)+".dat";
			const std::string fname = "tmp.dat";

			// Write file
			std::ofstream mf(fname);
			mf << "eta f f' f'' g g' rho u v w T\n";
			for (int i = 0; i<npts; ++i)
			{
				mf << output.eta[i] << " ";
				for (int n = 0; n<vec_t::size(); ++n)
				{
					mf << output.data[i][n] << " ";
				}
				for (int n = 0; n<vec_t::size(); ++n)
				{
					mf << output.data_dim[i][n] << " ";
				}
				mf << std::endl;
			}
		}
		
		void solve() const
		{
			output_type output(npts, deta);

			print("Solving system for ", inputs);

			// Initialize error
			real_t error1 = 100.0;
			real_t error2 = 100.0;
			int iter      = 0;
			real_t alpha  = alpha_ini;
			real_t beta   = beta_ini;

			// Solver loop
			while (((spade::utils::abs(error1) > TOL) || (spade::utils::abs(error2) > TOL)) && iter++<maxIter)
			{
				print("Iter, error1, error2 = ", iter, error1, error2);

				// BC
				set_boundary_condition(output, alpha, beta);

				// RK4 solver
				RK4_solver(output);

				// Error
				error1 = (output.data[npts-1][1] - 1.0);
				error2 = (output.data[npts-1][3] - 1.0);

				// Store baseline
				const real_t y2_old = output.data[npts-1][1];
				const real_t y4_old = output.data[npts-1][3];

				// Perturb BC
				set_boundary_condition(output, alpha+delta, beta);

                // RK4 solver
				RK4_solver(output);

				// Store first perturbation
				const real_t y2_n1 = output.data[npts-1][1];
				const real_t y4_n1 = output.data[npts-1][3];

				// Perturb BC
				set_boundary_condition(output, alpha, beta+delta);

                // RK4 solver
				RK4_solver(output);

				// Store first perturbation
				const real_t y2_n2 = output.data[npts-1][1];
				const real_t y4_n2 = output.data[npts-1][3];

				// Update solution
				const real_t p11  = (y2_n1 - y2_old) / delta;
				const real_t p21  = (y4_n1 - y4_old) / delta;
				const real_t p12  = (y2_n2 - y2_old) / delta;
				const real_t p22  = (y4_n2 - y4_old) / delta;
				const real_t r1   = 1.0 - y2_old;
				const real_t r2   = 1.0 - y4_old;
				const real_t idet = real_t(1.0) / (p11 * p22 - p12 * p21);
				alpha += (p22 * r1 - p12 * r2) * idet;
				beta  += (p11 * r2 - p21 * r1) * idet;
			}

			// Dimensionalize profile
			dimensionalize_profile(output);

			// Write output
			write_output(output);
		}
	};
}
