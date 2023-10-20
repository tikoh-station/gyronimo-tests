
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
#include <gyronimo/fields/equilibrium_helena.hh>
#include <gyronimo/fields/equilibrium_circular.hh>
#include <gyronimo/interpolators/bicubic_gsl.hh>

#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <gyronimo/dynamics/classical_boris.hh>
#include <gyronimo/dynamics/curvilinear_boris.hh>
#include <gyrotests/odeint_lorentz.hh>

#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <argh.h>
#include <gyrotests/guiding_centre_conversion.hh>


void print_help() {
	std::cout << "orbit_helena, powered by gyronimo-v"
		<< gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
	std::cout << "usage: orbit_helena [options] hmap\n";
	std::cout << "reads an HELENA hmap and prints the required orbit to stdout.\n";
	std::cout << "options:\n";
	std::cout << "  -stepper=val   Name of the stepper algorithm:\n";
	std::cout << "      + guiding_centre (default)\n";
	std::cout << "      + classical_boris\n";
	std::cout << "      + curvilinear_boris\n";
	std::cout << "      + lorentz\n";
	std::cout << "  -lref=val      Reference length   (in SI, default 1).\n";
	std::cout << "  -vref=val      Reference velocity (in SI, default 1).\n";
	std::cout << "  -mass=val      Particle mass   (in m_proton, default 1).\n";
	std::cout << "  -charge=val    Particle charge (in q_proton, default 1).\n";
	std::cout << "  -s=val         Flux surface label (default random).\n";
	std::cout << "  -pphi=val      Pphi value (in eV.s, default random).\n";
	std::cout << "  -chi=val       Poloidal angle (in rad, default random).\n";
	std::cout << "  -energy=val    Energy (in eV, default 2e4 eV).\n";
	std::cout << "  -pitch=val     Pitch angle (between [-1,1], default random).\n";
	std::cout << "  -lambda=val    Lambda value, signed as Vpp (default random).\n";
	std::cout << "  -gyrophase=val Gyrophase angle (default random).\n";
	std::cout << "  -Xgc           Initialize from Guiding-Centre position.\n";
	std::cout << "  -Xfo           Initialize from Full-Orbit position.\n";
	std::cout << "  -tfinal=val    Time limit (in Tref, default 1).\n";
	std::cout << "  -nsamples=val  Number of time samples (default 512).\n";
	std::cout << "  -seed=val      Random generator seed (default 0 sets CPU time).\n";
	std::exit(0);
}

// ODEInt observer object to print diagnostics at each time step.
class orbit_observer {
public:
	orbit_observer(double pstar, double vstar,
				const gyronimo::morphism* morph, 
				const gyronimo::guiding_centre* gc)
			: pstar_(pstar), vstar_(vstar), morph_(morph), 
			magnetic_field_(gc->magnetic_field()), gc_pointer_(gc) {
		std::cout.precision(16);
		std::cout.setf(std::ios::scientific);
	};
	void operator()(const gyronimo::guiding_centre::state& s, double t) {
		gyronimo::IR3 q = gc_pointer_->get_position(s);
		gyronimo::IR3 x = (*morph_)(q);
		double v_parallel = gc_pointer_->get_vpp(s);
		double bphi = magnetic_field_->covariant_versor(q, t)[gyronimo::IR3::w];
		double flux = q[gyronimo::IR3::u]*q[gyronimo::IR3::u];
        double eper = gc_pointer_->energy_perpendicular(s, t);
        double epar = gc_pointer_->energy_parallel(s);
		std::cout << t << ' '
			<< x[gyronimo::IR3::u] << ' '
			<< x[gyronimo::IR3::v] << ' '
			<< x[gyronimo::IR3::w] << ' '
			<< q[gyronimo::IR3::u] << ' '
			<< q[gyronimo::IR3::v] << ' '
			<< q[gyronimo::IR3::w] << ' '
			<< -pstar_*flux + vstar_*v_parallel*bphi << ' '
			<< eper+epar << ' '
			<< eper << ' '
			<< epar << ' '
			<< gc_pointer_->mu_tilde() << '\n';
	};
private:
	double pstar_, vstar_;
	const gyronimo::morphism* morph_;
	const gyronimo::IR3field* magnetic_field_;
	const gyronimo::guiding_centre* gc_pointer_;
};

// ODEInt observer object to print diagnostics at each time step.
template<typename lorentz_handler>
class full_orbit_observer {
public:
	full_orbit_observer(double pstar, double vstar,
				const gyronimo::morphism* morph, 
				const lorentz_handler* lo)
			: pstar_(pstar), vstar_(vstar), morph_(morph), 
			magnetic_field_(lo->magnetic_field()), lo_(lo) {
		std::cout.precision(16);
		std::cout.setf(std::ios::scientific);
	};
	~full_orbit_observer() {};
	void operator()(const lorentz_handler::state& s, double t) {
		gyronimo::IR3 q = lo_->get_position(s);
		gyronimo::IR3 v_con = lo_->get_velocity(s);
		gyronimo::IR3 v = morph_->from_contravariant(v_con, q);
		double v_cov_phi = morph_->to_covariant(v, q)[gyronimo::IR3::w];
		gyronimo::IR3 x = (*morph_)(q);
		double Bmag = magnetic_field_->magnitude(q, t);
		double flux = q[gyronimo::IR3::u]*q[gyronimo::IR3::u];
		double Eperp = lo_->energy_perpendicular(s, t);
		std::cout << t << ' '
			<< x[gyronimo::IR3::u] << ' '
			<< x[gyronimo::IR3::v] << ' '
			<< x[gyronimo::IR3::w] << ' '
			<< q[gyronimo::IR3::u] << ' '
			<< q[gyronimo::IR3::v] << ' '
			<< q[gyronimo::IR3::w] << ' '
			<< -pstar_*flux + vstar_*v_cov_phi << ' '
			<< lo_->energy_kinetic(s) << ' '
			<< Eperp << ' '
			<< lo_->energy_parallel(s, t) << ' '
			<< Eperp/Bmag << '\n';
	};
private:
	double pstar_, vstar_;
	const gyronimo::morphism* morph_;
	const gyronimo::IR3field* magnetic_field_;
	const lorentz_handler* lo_;
};

double match_target(double target, std::function<double(double)> &functional) {
	auto orbit = [&functional, target](double s) {
		return functional(s) - target;};
	auto s_grid = gyronimo::linspace<std::valarray<double>>(0.0, 1.0, 1024);
	auto bracketing_iterator = std::adjacent_find(
		std::begin(s_grid), std::end(s_grid),
		[&orbit](double x, double y){return orbit(x)*orbit(y) < 0.0;});
	if(bracketing_iterator == std::end(s_grid)) {
		std::cerr << "# solution not found.\n";
		std::exit(1);
	}
	auto root_interval = boost::math::tools::bisect(
		[&orbit](double s){return orbit(s);},
		*bracketing_iterator, *(bracketing_iterator + 1),
		[](double a, double b){return std::abs(b - a) < 1.0e-09;});
	return 0.5*(root_interval.first + root_interval.second);
}

template<typename stepperclass, typename observer_t>
void integrate_trajectory(const stepperclass *stepper, 
		typename stepperclass::state init, double ti,
		size_t Nsteps, double timestep, observer_t &observer) {

	auto begin = std::chrono::high_resolution_clock::now();  // starts ticking...

	// run simulation
	double time = ti;
	typename stepperclass::state s = init;
	observer(s, time);
	for(size_t i = 1; i <= Nsteps; ++i) {

		s = stepper->do_step(s, time, timestep);
		time = ti + i * timestep;
		observer(s, time);
	}

	auto end = std::chrono::high_resolution_clock::now();  // stops ticking...
	auto elapsed_mseconds =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cerr << elapsed_mseconds.count()/1000.0 << " sec.\n";
	return;
}

int main(int argc, char* argv[]) {
	auto command_line = argh::parser(argv);
	if (command_line[{"h", "help"}]) print_help();
	if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
		std::cout << "orbit_helena: no helena mapping file provided; -h for help.\n";
		std::exit(1);
	}

	/* Initialize helena equilibrium */
	gyronimo::parser_helena hmap(command_line[1]);
	gyronimo::bicubic_gsl_factory ifactory(false,
		(hmap.is_symmetric() ? 0 : 9), (hmap.is_symmetric() ? 9 : 0));
	gyronimo::morphism_cache<gyronimo::morphism_helena> m(&hmap, &ifactory);
	gyronimo::metric_cache<gyronimo::metric_helena> g(&m, &ifactory);
	gyronimo::IR3_field_c1_cache<gyronimo::equilibrium_helena> heq(&g, &ifactory);

	/* Reads parameters from the command line: */
    size_t seed; command_line("seed", 0) >> seed;
	if(seed == 0) {
		auto time = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::hours>(time.time_since_epoch());
		seed = duration.count();
	}
	std::mt19937 rng(seed);
	std::uniform_real_distribution uniform_01(0.0, 1.0);
	std::uniform_real_distribution uniform_angle(0.0, 2*std::numbers::pi);

	std::string steppername; command_line("stepper", "guiding_centre") >> steppername;
	double Lref;   command_line("lref",   1.0) >> Lref;   // SI units.
	double Vref;   command_line("vref",   1.0) >> Vref;   // SI units.
	double mass;   command_line("mass",   1.0) >> mass;   // m_proton units.
	double charge; command_line("charge", 1.0) >> charge; // q_electron units.
	double energy; command_line("energy", 3.5e6) >> energy; // energy in eV.

    // Uref in eV
    double moe = gyronimo::codata::m_proton / gyronimo::codata::e;
	double Uref = 0.5*moe*mass*Vref*Vref;
	double energy_tilde = energy / Uref;
	double vel_tilde = std::sqrt(energy_tilde); // particle velocity (in Vref)
	double qom = charge / mass;
    double pstar = charge * hmap.cpsurf()*heq.B0()*heq.R0()*heq.R0();
	double vstar = Vref*mass*moe;
	double vdagger = Vref*mass*moe*vel_tilde;

    double chi; std::string chi_str; command_line("chi", "random") >> chi_str;
	if(chi_str != "random") command_line("chi", 0.0) >> chi; // normalized units.
	else chi = uniform_angle(rng);
    
    double s; std::string s_str; command_line("s", "random") >> s_str;
    double pphi; std::string pphi_str; command_line("pphi", "random") >> pphi_str; // pphi in eV.s.
	double pitch; std::string pitch_str; command_line("pitch", "random") >> pitch_str;
	double lambda; std::string lambda_str; command_line("lambda", "random") >> lambda_str;
	if(s_str != "random" || pphi_str == "random") {
        if(s_str != "random") command_line("s", -1.0) >> s;
        else s = uniform_01(rng);
        if(pitch_str != "random" || lambda_str == "random") {
            if(pitch_str != "random") command_line("pitch", 0.5) >> pitch;
            else pitch = 2.0 * uniform_01(rng)-1.0;
            lambda = (1-pitch*pitch)/heq.magnitude({s, chi, 0}, 0.0);
        } else if (lambda_str != "random") {
            command_line("lambda", 0.5) >> lambda;
            double vpp_sign = std::copysign(1.0, lambda);
            lambda = std::abs(lambda);
            pitch = vpp_sign*std::sqrt(1 - lambda*heq.magnitude({s, chi, 0}, 0.0));
        }
        double bphi = heq.covariant_versor({s, chi, 0.0}, 0.0)[gyronimo::IR3::w];
        pphi = -pstar*s*s + vdagger*pitch*bphi;
    } else {
        command_line("pphi", 0.0) >> pphi;
        if(pitch_str != "random" || lambda_str == "random") {
            if(pitch_str != "random") command_line("pitch", 0.5) >> pitch;
            else pitch = 2.0 * uniform_01(rng)-1.0;

            std::function<double(double)> pphi_functional = [&](double s_) {
                double bphi = heq.covariant_versor({s_, chi, 0.0}, 0.0)[gyronimo::IR3::w];
                double pphi_ = -pstar*s_*s_ + vdagger*pitch*bphi;
                return pphi_;
            };

            s = match_target(pphi, pphi_functional);
            lambda = (1-pitch*pitch)/heq.magnitude({s, chi, 0}, 0.0);

        } else {
            command_line("lambda", 0.5) >> lambda;
            double vpp_sign_ = std::copysign(1.0, lambda);
            lambda = std::abs(lambda);

            std::function<double(double)> pphi_functional = [&](double s_) {
                double Btilde = heq.magnitude({s_, chi, 0.0}, 0.0);
                double bphi = heq.covariant_versor({s_, chi, 0.0}, 0.0)[gyronimo::IR3::w];
                double pphi_ = -pstar*s_*s_ + vpp_sign_*vdagger*bphi*std::sqrt(1-lambda*Btilde);
                return pphi_;
            };

            s = match_target(pphi, pphi_functional);
            pitch = vpp_sign_*std::sqrt(1 - lambda*heq.magnitude({s, chi, 0}, 0.0));
        }
    }
	double vpp_sign = std::copysign(1.0, pitch);
    gyronimo::IR3 pos = {s, chi, 0.0};

    double gyrophase; std::string gyrophase_str; command_line("gyrophase", "random") >> gyrophase_str;
	if(gyrophase_str != "random") command_line("gyrophase", 0.0) >> gyrophase;
	else gyrophase = uniform_angle(rng);

	double Tfinal; command_line("tfinal", 1.0) >> Tfinal; // in Tref units
	size_t nsamples; command_line("nsamples", 512) >> nsamples;
	double dt = Tfinal / nsamples;

	/* Print output header */
	std::cout << "# orbit_helena, powered by ::gyronimo:: v"
		<< gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
	std::cout << "# args: ";
	for(int i = 1; i < argc; i++) std::cout << argv[i] << " ";
	std::cout << std::endl;
	std::cout << "# l_ref = " << Lref << " [m];";
	std::cout << " v_ref = " << Vref << " [m/s];";
	std::cout << " u_ref = " << Uref << " [eV];";
	std::cout << " energy = " << energy << " [eV];";
	std::cout << " pphi = " << pphi << " [eV.s];";
	std::cout << " vpp_sign = " << vpp_sign << " ;";
	std::cout << " lambda = " << lambda << " \n";
	std::cout << "# vars: t x y z s chi phi Pphi/e Ekin/Uref Eperp/Uref Epar/Uref mu*Bref/Uref\n";


	/* Integrate trajectory */
	if(steppername == "guiding_centre") {

	    gyronimo::IR3 Xgc = convert_to_GC(Lref, Vref, qom, pos, energy_tilde, 
			pitch, gyrophase, &heq, &m, command_line[{"Xfo"}]);
	    double mu_tilde = lambda * energy_tilde;
    
		gyronimo::guiding_centre gc(Lref, Vref, qom, mu_tilde, &heq, nullptr);
		gyronimo::guiding_centre::state initial_state = gc.generate_state(
			Xgc, energy_tilde, (vpp_sign > 0 ? 
			gyronimo::guiding_centre::plus : gyronimo::guiding_centre::minus), 0);

		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples, using RK4.
		orbit_observer observer(pstar, vstar, &m, &gc);
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
		pitch, gyrophase, &heq, &m, command_line[{"Xgc"}]);

	if(steppername == "classical_boris") {

		gyronimo::classical_boris boris(Lref, Vref, qom, &heq, nullptr);
		gyronimo::classical_boris::state initial_state = 
			boris.half_back_step(q0, v0, 0.0, dt);
		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples.
		full_orbit_observer<gyronimo::classical_boris> observer(pstar, vstar, &m, &boris);
		integrate_trajectory(&boris, initial_state, 0.0, nsamples, dt, observer);

	} else if(steppername == "curvilinear_boris") {

		gyronimo::curvilinear_boris boris(Lref, Vref, qom, &heq, nullptr);
		gyronimo::curvilinear_boris::state initial_state = 
			boris.half_back_step(q0, v0, 0.0, dt);

		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples.
		full_orbit_observer<gyronimo::curvilinear_boris> observer(pstar, vstar, &m, &boris);
		integrate_trajectory(&boris, initial_state, 0.0, nsamples, dt, observer);

	} else if(steppername == "lorentz") {

		toolbox::rk4_lorentz rk4(Lref, Vref, qom, &heq, nullptr);
		toolbox::rk4_lorentz::state initial_state = 
			rk4.generate_initial_state(q0, v0, 0.0, dt);
		
		// integrates for t in [0,Tfinal], with dt=Tfinal/nsamples.
		full_orbit_observer<toolbox::rk4_lorentz> observer(pstar, vstar, &m, &rk4);
		integrate_trajectory(&rk4, initial_state, 0.0, nsamples, dt, observer);

	} else {
		std::cout << "orbit_helena: stepper name not recognized; -h for help.\n";
		std::exit(1);
	}

	return 0;
}