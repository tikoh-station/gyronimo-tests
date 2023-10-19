#ifndef GUIDING_CENTRE_CONVERSION
#define GUIDING_CENTRE_CONVERSION

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/morphism.hh>

gyronimo::IR3 convert_to_GC(double Lref, double Vref, double qom, 
	gyronimo::IR3 pos, double energy, double pitch, 
	double gyrophase, const gyronimo::IR3field *magnetic_field, 
	const gyronimo::morphism *morph, bool from_FO);

std::tuple<gyronimo::IR3,gyronimo::IR3> convert_to_FO(double Lref, double Vref, double qom, 
	gyronimo::IR3 pos, double energy_tilde, double pitch, 
	double gyrophase, const gyronimo::IR3field *magnetic_field, 
	const gyronimo::morphism *morph, bool from_GC);

#endif // GUIDING_CENTRE_CONVERSION