#include "pptb.h"
#include "spade.h"
#include "scidf.h"

using real_t = double;

int main(int argc, char** argv)
{
	// Set solver parameters
	pptb::analysis::solver_inputs_t<real_t> inputs;
	inputs.Mae = 6.0;
	inputs.Te  = 500.0;
	inputs.rhoe= 0.01;
	inputs.Tw  = 300.0;

	// Initialize solver
	pptb::analysis::similarity_solver_t<real_t> similarity(inputs);

	// Solve system
	similarity.solve();
	
	return 0;
}
