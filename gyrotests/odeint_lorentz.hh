#ifndef TOOLBOX_ODEINT_LORENTZ
#define TOOLBOX_ODEINT_LORENTZ

#include <gyronimo/core/error.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <boost/numeric/odeint.hpp>

namespace toolbox {

template<typename odeint_stepper>
class odeint_lorentz {
public:
	using state = gyronimo::lorentz::state;

	odeint_lorentz(const double Lref, const double Vref, const double qom, 
		const gyronimo::IR3field *electric, const gyronimo::IR3field *magnetic);
	~odeint_lorentz() {};

	state do_step(const state &in, const double t, const double dt) const;

	state generate_state(const gyronimo::IR3 &position, const gyronimo::IR3 &velocity) const;
	gyronimo::IR3 get_position(const state& s) const;
	gyronimo::IR3 get_velocity(const state& s) const;

	double energy_kinetic(const state& s) const;
	double energy_parallel(const state& s, const double &time) const;
	double energy_perpendicular(const state& s, const double &time) const;

	state generate_initial_state(const gyronimo::IR3 &q, 
		const gyronimo::IR3 &cartesian_velocity, 
		const double &tinit, const double &dt) const;

	//! Returns the pointer to the electric field.
	const gyronimo::IR3field* electric_field() const {return lo_.electric_field();};
	//! Returns the pointer to the magnetic field.
	const gyronimo::IR3field* magnetic_field() const {return lo_.magnetic_field();};
	//! Returns the pointer to the morphism.
	const gyronimo::morphism* my_morphism() const {return my_morphism_;};

	//! Returns the charge over mass `qom`.
	double qom() const {return qom_;};
	//! Returns the reference gyration frequency scale `Oref`.
	double Oref() const {return lo_.Oref_tilde();};
	//! Returns the reference length scale `Lref`.
	double Lref() const {return lo_.Lref();};
	//! Returns the reference time scale `Tref`.
	double Tref() const {return lo_.Tref();};
	//! Returns the reference velocity scale `Vref`.
	double Vref() const {return lo_.Vref();};

private:
	const double qom_;
	const gyronimo::metric_connected *metric_;
	const gyronimo::morphism *my_morphism_;
	const gyronimo::lorentz lo_;
	const gyronimo::odeint_adapter<gyronimo::lorentz> sys_;

};	// end class odeint_lorentz

using rk4_lorentz = odeint_lorentz<boost::numeric::odeint::runge_kutta4<gyronimo::lorentz::state>>;
using rkck54_lorentz = odeint_lorentz<boost::numeric::odeint::runge_kutta_cash_karp54<gyronimo::lorentz::state>>;



template<typename odeint_stepper>
odeint_lorentz<odeint_stepper>::odeint_lorentz(const double Lref, const double Vref, const double qom, 
		const gyronimo::IR3field *electric, const gyronimo::IR3field *magnetic) 
			: qom_(qom), metric_(magnetic ? dynamic_cast<const gyronimo::metric_connected*>(magnetic->metric()) : nullptr),
			my_morphism_(metric_ ? metric_->my_morphism() : nullptr),
			lo_(Lref, Vref, qom, electric, magnetic), sys_(&lo_) {

	if(!metric_) gyronimo::error(__func__, __FILE__, __LINE__, 
		" field metric must inherit from 'metric_connected'.", 1);
}

//! Performs a time step `dt` update to the state `in` and returns the result.
template<typename odeint_stepper>
odeint_lorentz<odeint_stepper>::state odeint_lorentz<odeint_stepper>::do_step(
		const state &in, const double t, const double dt) const {
	odeint_stepper stepper;
	odeint_lorentz<odeint_stepper>::state out(in);
	stepper.do_step(sys_, out, t, dt);
	return out;
}

//! Returns the vector position of the state normalized to `Lref`.
template<typename odeint_stepper>
gyronimo::IR3 odeint_lorentz<odeint_stepper>::get_position(
		const state& s) const {
	return lo_.get_position(s);
}
//! Returns the vector velocity of the state normalized to `Vref`.
template<typename odeint_stepper>
gyronimo::IR3 odeint_lorentz<odeint_stepper>::get_velocity(
		const state& s) const {
	return lo_.get_velocity(s);
}
//! Returns the kinetic energy of the state, normalized to `Uref`.
template<typename odeint_stepper>
double odeint_lorentz<odeint_stepper>::energy_kinetic(
		const state& s) const {
	return lo_.energy_kinetic(s);
}
template<typename odeint_stepper>
double odeint_lorentz<odeint_stepper>::energy_parallel(
		const state& s, const double &time) const {
	return lo_.energy_parallel(s, time);
}
template<typename odeint_stepper>
double odeint_lorentz<odeint_stepper>::energy_perpendicular(
		const state& s, const double &time) const {
	return lo_.energy_perpendicular(s, time);
}

//! Returns the `exact_stepper::state` from a normalized point in phase-space.
template<typename odeint_stepper>
odeint_lorentz<odeint_stepper>::state odeint_lorentz<odeint_stepper>::generate_state(
		const gyronimo::IR3 &position, const gyronimo::IR3 &velocity) const {
	return lo_.generate_state(position, velocity);
}

//! Creates the first `exact_stepper::state` from a normalized point in cartesian phase-space.
template<typename odeint_stepper>
odeint_lorentz<odeint_stepper>::state odeint_lorentz<odeint_stepper>::generate_initial_state(
		const gyronimo::IR3 &q, const gyronimo::IR3 &cartesian_velocity, 
		const double &tinit, const double &dt) const {
	gyronimo::IR3 vel = my_morphism_->to_contravariant(cartesian_velocity, q);
	return lo_.generate_state(q, vel);
}

} // end namespace toolbox

#endif // TOOLBOX_ODEINT_LORENTZ