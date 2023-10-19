#include "guiding_centre_conversion.hh"

#include <cmath>
#include <gyronimo/core/codata.hh>

gyronimo::IR3 convert_to_GC(double Lref, double Vref, double qom, 
	gyronimo::IR3 pos, double energy_tilde, double pitch, 
	double gyrophase, const gyronimo::IR3field *magnetic_field, 
	const gyronimo::morphism *morph, bool from_FO) {

	if(!from_FO) return pos;

	double sin_phi = std::sin(gyrophase);
	double cos_phi = std::cos(gyrophase);

	gyronimo::IR3 Bversor_con = magnetic_field->contravariant_versor(pos, 0.0);
	double Bmag = magnetic_field->magnitude(pos, 0.0);
	double Bref = magnetic_field->m_factor();
	gyronimo::IR3 Bversor = morph->from_contravariant(Bversor_con, pos);
	double bx = Bversor[gyronimo::IR3::u];
	double by = Bversor[gyronimo::IR3::v];
	double bz = Bversor[gyronimo::IR3::w];
	double r = std::sqrt(bx*bx+by*by);
	double ca = bz, sa = r;
	double cb = by/r, sb = bx/r;
	gyronimo::IR3 rho_versor = {
		 cos_phi*cb + sin_phi*ca*sb,
		-cos_phi*sb + sin_phi*ca*cb,
		-sin_phi*sa
	};
	double rho_ref = (gyronimo::codata::m_proton/gyronimo::codata::e)*Vref/Bref;
	double rho_mag = rho_ref * std::sqrt(energy_tilde) / (std::abs(qom) * Bmag);
	gyronimo::IR3 rho = rho_mag * rho_versor;

    gyronimo::IR3 posGC = pos;
	if(from_FO) posGC = morph->translation(pos, (-1)*rho);

	return posGC;
}

std::tuple<gyronimo::IR3,gyronimo::IR3> convert_to_FO(double Lref, double Vref, double qom, 
	gyronimo::IR3 pos, double energy_tilde, double pitch, 
	double gyrophase, const gyronimo::IR3field *magnetic_field, 
	const gyronimo::morphism *morph, bool from_GC) {

	double sin_phi = std::sin(gyrophase);
	double cos_phi = std::cos(gyrophase);

	gyronimo::IR3 Bversor_con = magnetic_field->contravariant_versor(pos, 0.0);
	double Bmag = magnetic_field->magnitude(pos, 0.0);
	double Bref = magnetic_field->m_factor();
	gyronimo::IR3 Bversor = morph->from_contravariant(Bversor_con, pos);
	double bx = Bversor[gyronimo::IR3::u];
	double by = Bversor[gyronimo::IR3::v];
	double bz = Bversor[gyronimo::IR3::w];
	double r = std::sqrt(bx*bx+by*by);
	double ca = bz, sa = r;
	double cb = by/r, sb = bx/r;
	gyronimo::IR3 rho_versor = {
		 cos_phi*cb + sin_phi*ca*sb,
		-cos_phi*sb + sin_phi*ca*cb,
		-sin_phi*sa
	};
	gyronimo::IR3 vperp_versor = {
		 sin_phi*cb - cos_phi*ca*sb,
		-sin_phi*sb - cos_phi*ca*cb,
		cos_phi*sa
	};
	double rho_ref = (gyronimo::codata::m_proton/gyronimo::codata::e)*Vref/Bref;
	double rho_mag = rho_ref * std::sqrt(energy_tilde) / (std::abs(qom) * Bmag);
	gyronimo::IR3 rho = rho_mag * rho_versor;

	double vel_tilde = std::sqrt(energy_tilde);
	double vpar_tilde = vel_tilde * pitch;
	double vperp_tilde = vel_tilde * std::sqrt(1-pitch*pitch);
	double sign = std::copysign(1.0, qom);

	gyronimo::IR3 posFO = pos;
	if(from_GC) posFO = morph->translation(pos, rho);
	gyronimo::IR3 velFO = vpar_tilde*Bversor + (sign*vperp_tilde)*vperp_versor;

	return std::tuple<gyronimo::IR3, gyronimo::IR3>(posFO, velFO);
}