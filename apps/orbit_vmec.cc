
#include <cmath>
#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <numbers>

#include <gyronimo/version.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>

#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <gyronimo/dynamics/classical_boris.hh>
#include <gyronimo/dynamics/curvilinear_boris.hh>
#include <gyrotests/odeint_lorentz.hh>

#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <gsl/gsl_errno.h>

#include <argh.h>
#include <gyrotests/guiding_centre_conversion.hh>

void print_help() {
	std::cout << "orbit_vmec, powered by ::gyronimo::v"
		<< gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
	std::cout << "usage: orbit_vmec [options] vmec_netcdf_file\n";
	std::cout <<
		"reads a vmec output file, prints the required orbit to stdout.\n";
	std::cout << "options:\n";
	std::cout << "  -stepper=sss   Name of the stepper algorithm:\n";
	std::cout << "      + guiding_centre (default)\n";
	std::cout << "      + classical_boris\n";
	std::cout << "      + curvilinear_boris\n";
	std::cout << "      + lorentz\n";
	std::cout << "      + boris_advanced\n";
	std::cout << "  -lref=val      Reference length   (in m,   default 1).\n";
	std::cout << "  -vref=val      Reference velocity (in m/s, default 1).\n";
	std::cout << "  -mass=val      Particle mass (in m_proton, default 1).\n";
	std::cout << "  -charge=val    Particle charge (in q_proton, default 1).\n";
	std::cout << "  -seed=val      Seed for the RNG (default CPU time)\n";
	std::cout << "      + randomize 'flag' using: -flag=random \n";
	std::cout << "      + randomizable: flux, zeta, theta \n";
	std::cout << "  -flux=val      Initial toroidal flux (vmec, default random).\n";
	std::cout << "  -zeta=val      Initial zeta  (vmec in rad, default random).\n";
	std::cout << "  -theta=val     Initial theta (vmec in rad, default random).\n";
	std::cout << "  -energy=val    Energy value (in eV, default 10000).\n";
	std::cout << "  -lambda=val    Lambda value, signed as v_par (default random).\n";
	std::cout << "  -pitch=val     Pitch angle, between 1 and -1 (default random).\n";
	std::cout << "  -gyrophase=val Gyrophase angle (default random).\n";
	std::cout << "  -Xgc           Initialize from Guiding-Centre position.\n";
	std::cout << "  -Xfo           Initialize from Full-Orbit position.\n";
	std::cout << "  -tfinal=val    Time limit (in Tref, default 500).\n";
	std::cout << "  -nsamples=val  Number of orbit samples (default 512).\n";
	std::cout << "Note: lambda=magnetic_moment_si*B_axis_si/energy_si.\n";
	std::exit(0);
}

// ODEInt observer object to print diagnostics at each time step.
class orbit_observer {
public:
	orbit_observer(const gyronimo::morphism_vmec* morph, 
				const gyronimo::guiding_centre* gc)
			: morph_(morph), gc_pointer_(gc) {
		std::cout.precision(16);
		std::cout.setf(std::ios::scientific);
	};
	void operator()(const gyronimo::guiding_centre::state& s, double t) {
		gyronimo::IR3 q = gc_pointer_->get_position(s);
		gyronimo::IR3 x = (*morph_)(q);
        double Epar = gc_pointer_->energy_parallel(s);
        double Eper = gc_pointer_->energy_perpendicular(s, t);
		std::cout << t << ' '
			<< x[gyronimo::IR3::u] << ' '
			<< x[gyronimo::IR3::v] << ' '
			<< x[gyronimo::IR3::w] << ' '
			<< q[gyronimo::IR3::u] << ' '
			<< q[gyronimo::IR3::v] << ' '
			<< q[gyronimo::IR3::w] << ' '
			<< Epar+Eper << ' '
			<< Eper << ' '
			<< Epar << ' '
			<< gc_pointer_->mu_tilde() << '\n';
	};
private:
	const gyronimo::morphism_vmec* morph_;
	const gyronimo::guiding_centre* gc_pointer_;
};

// ODEInt observer object to print diagnostics at each time step.
template<typename lorentz_handler>
class full_orbit_observer {
public:
	full_orbit_observer(const gyronimo::morphism_vmec* morph, const lorentz_handler* lo)
			: morph_(morph), lo_(lo) {
		std::cout.precision(16);
		std::cout.setf(std::ios::scientific);
	};
	~full_orbit_observer() {};
	void operator()(const lorentz_handler::state& s, double t) {
		gyronimo::IR3 q = lo_->get_position(s);
		gyronimo::IR3 x = (*morph_)(q);
		double Bmag = lo_->magnetic_field()->magnitude(q, t);
		double Eperp = lo_->energy_perpendicular(s, t);
		std::cout << t << ' '
			<< x[gyronimo::IR3::u] << ' '
			<< x[gyronimo::IR3::v] << ' '
			<< x[gyronimo::IR3::w] << ' '
			<< q[gyronimo::IR3::u] << ' '
			<< q[gyronimo::IR3::v] << ' '
			<< q[gyronimo::IR3::w] << ' '
            << lo_->energy_kinetic(s) << ' '
			<< Eperp << ' '
			<< lo_->energy_parallel(s, t) << ' '
			<< Eperp/Bmag << '\n';
	};
private:
	const gyronimo::morphism_vmec* morph_;
	const lorentz_handler* lo_;
};

template<typename stepperclass, typename observer_t>
void integrate_trajectory(const stepperclass *stepper, 
		typename stepperclass::state init, double ti,
		size_t Nsteps, double timestep, observer_t &observer) {
			
	auto begin = std::chrono::high_resolution_clock::now();  // starts ticking...

	// run simulation
	double time = ti;
	typename stepperclass::state s = init;
	observer(s, time);
    bool escaped = false;
	for(size_t i = 1; i <= Nsteps; ++i) {

        try {
		    s = stepper->do_step(s, time, timestep);
        } catch(std::domain_error &e) {
            escaped = true;
            break;
        }
		time = ti + i * timestep;
		observer(s, time);
	}
    if(escaped) std::cerr << "Particle escaped confinement" << std::endl;

	auto end = std::chrono::high_resolution_clock::now();  // stops ticking...
	auto elapsed_mseconds =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cerr << elapsed_mseconds.count()/1000.0 << " sec.\n";


	return;
}

void new_handler(const char *reason, const char *file, int line, int gsl_errno) {
	if(gsl_errno == GSL_EDOM) {
		throw std::domain_error("Outside domain of interpolation.");
	}
	std::cout << "gsl: " << file << ":" << line << ": ERROR: " << reason << "\n";
	std::cout << "Custom GSL error handler invoked.\nAborted\n";
	std::exit(1);
	return;
}



int main(int argc, char* argv[]) {
	auto command_line = argh::parser(argv);
	if (command_line[{"h", "help"}]) print_help();
	if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
		std::cout << "orbit_vmec: no vmec mapping file provided; -h for help.\n";
		std::exit(1);
	}
    
	gsl_set_error_handler_off();
	gsl_set_error_handler(&new_handler);
	
	/* Initialize helena equilibrium */
	gyronimo::parser_vmec vmap(command_line[1]);
	gyronimo::cubic_gsl_factory ifactory;
	gyronimo::morphism_vmec m(&vmap, &ifactory);
	gyronimo::metric_vmec g(&m);
	gyronimo::equilibrium_vmec veq(&g, &ifactory);

	size_t seed; command_line("seed", 0) >> seed;
	if(seed == 0) {
		auto time = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::hours>(time.time_since_epoch());
		seed = duration.count();
	}
	std::mt19937 rng(seed);
	std::uniform_real_distribution uniform_01(0.0, 1.0);
	std::uniform_real_distribution uniform_angle(0.0, 2*std::numbers::pi);

	/* Reads parameters from the command line: */
	std::string steppername; command_line("stepper", "guiding_centre") >> steppername;
	double Lref; command_line("lref", 1.0) >> Lref;  // Lref in m.
	double Vref; command_line("vref", 1.0) >> Vref;  // Vref in m/s.
	double mass;   command_line("mass",   1.0) >> mass;   // m_proton units.
	double charge; command_line("charge", 1.0) >> charge; // q_electron units.

	double flux; std::string flux_str; command_line("flux", "random") >> flux_str;
	if(flux_str != "random") command_line("flux", 0.2) >> flux; // normalized units.
	else flux = uniform_01(rng);

	double zeta; std::string zeta_str; command_line("zeta", "random") >> zeta_str;
	if(zeta_str != "random") command_line("zeta", 0.0) >> zeta; // angle units.
	else zeta = uniform_angle(rng);

	double theta; std::string theta_str; command_line("theta", "random") >> theta_str;
	if(theta_str != "random") command_line("theta", 0.0) >> theta; // angle units.
	else theta = uniform_angle(rng);

	gyronimo::IR3 pos = {flux, zeta, theta};
	double energy; command_line("energy", 10000.0) >> energy; // energy in eV.

	double pitch; std::string pitch_str; command_line("pitch", "random") >> pitch_str;
	double lambda; std::string lambda_str; command_line("lambda", "random") >> lambda_str;
	if(pitch_str != "random") {
		command_line("pitch", 0.5) >> pitch;
		lambda = (1 - pitch*pitch) / veq.magnitude(pos, 0.0);
	} else if (lambda_str != "random") {
		command_line("lambda", 0.5) >> lambda;
		double vpp_sign = std::copysign(1.0, lambda);
		lambda = std::abs(lambda);
		pitch = vpp_sign*std::sqrt(1 - lambda*veq.magnitude(pos, 0.0));
	} else {
		pitch = 2.0 * uniform_01(rng) - 1.0;
		lambda = (1 - pitch*pitch) / veq.magnitude(pos, 0.0);
	}

	double gyrophase; std::string phase_str; command_line("gyrophase", "random") >> phase_str;
	if(phase_str != "random") command_line("gyrophase", 0.0) >> gyrophase; // angle units.

	double Tfinal; command_line("tfinal", 500.0) >> Tfinal;
	size_t nsamples; command_line("nsamples", 512) >> nsamples;
	double dt = Tfinal/nsamples;

	/* Computes normalisation constants: */
	double Uref = 0.5*gyronimo::codata::m_proton*mass*Vref*Vref/gyronimo::codata::e; // Uref in eV
	double energy_tilde = energy / Uref;
	double qom = charge / mass;

	/* Print output header */
	std::cout << "# orbit_vmec, powered by ::gyronimo:: v"
		<< gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
	std::cout << "# args: ";
	for(int i = 1; i < argc; i++) std::cout << argv[i] << " ";
	std::cout << std::endl;
	std::cout << "# l_ref = " << Lref << " [m];";
	std::cout << " v_ref = " << Vref << " [m/s];";
	std::cout << " u_ref = " << Uref << " [eV];";
	std::cout << " energy = " << energy << " [eV]." << "\n";
	std::cout << "# vars: t x y z s phi theta Ekin/Uref Eperp/Uref Epar/Uref mu*Bref/Uref\n";

	/* Integrate trajectory */
	if(steppername == "guiding_centre") {
 
		/* Guiding centre initial position */
		gyronimo::IR3 Xgc = convert_to_GC(Lref, Vref, qom, pos, energy_tilde, 
			pitch, gyrophase, &veq, &m, command_line[{"Xfo"}]);
		double mu_tilde = lambda * energy_tilde;

		gyronimo::guiding_centre gc(Lref, Vref, qom, mu_tilde, &veq, nullptr);
		gyronimo::guiding_centre::state initial_state = gc.generate_state(
			Xgc, energy_tilde, (pitch > 0 ? 
			gyronimo::guiding_centre::plus : gyronimo::guiding_centre::minus), 0);

		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples, using RK4.
		orbit_observer observer(&m, &gc);
		boost::numeric::odeint::runge_kutta4<gyronimo::guiding_centre::state>
			integration_algorithm;
		
		auto begin = std::chrono::high_resolution_clock::now();  // starts ticking...
		boost::numeric::odeint::integrate_const(
			integration_algorithm, gyronimo::odeint_adapter(&gc),
			initial_state, 0.0, Tfinal, dt, observer);
		auto end = std::chrono::high_resolution_clock::now();  // stops ticking...
		auto elapsed_mseconds =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cerr << elapsed_mseconds.count()/1000.0 << " sec.\n";

		return 0;
	}

	auto [q0, v0] = convert_to_FO(Lref, Vref, qom, pos, energy_tilde, 
		pitch, gyrophase, &veq, &m, command_line[{"Xgc"}]);

	if(steppername == "classical_boris") {

		gyronimo::classical_boris boris(Lref, Vref, qom, &veq, nullptr);
		gyronimo::classical_boris::state initial_state = 
			boris.half_back_step(q0, v0, 0.0, dt);

		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples, using RK4.
		full_orbit_observer<gyronimo::classical_boris> observer(&m, &boris);
		integrate_trajectory(&boris, initial_state, 0.0, nsamples, dt, observer);

	} else if(steppername == "curvilinear_boris") {

		gyronimo::curvilinear_boris boris(Lref, Vref, qom, &veq, nullptr);
		gyronimo::curvilinear_boris::state initial_state = 
			boris.half_back_step(q0, v0, 0.0, dt);

		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples, using RK4.
		full_orbit_observer<gyronimo::curvilinear_boris> observer(&m, &boris);
		integrate_trajectory(&boris, initial_state, 0.0, nsamples, dt, observer);

	} else if(steppername == "lorentz") {

		toolbox::rk4_lorentz rk4(Lref, Vref, qom, &veq, nullptr);
		toolbox::rk4_lorentz::state initial_state = 
			rk4.generate_initial_state(q0, v0, 0.0, dt);
		
		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples, using RK4.
		full_orbit_observer<toolbox::rk4_lorentz> observer(&m, &rk4);
		integrate_trajectory(&rk4, initial_state, 0.0, nsamples, dt, observer);

	} else {
		std::cout << "orbit_vmec: stepper name not recognized; -h for help.\n";
		std::exit(1);
	}

	return 0;
}