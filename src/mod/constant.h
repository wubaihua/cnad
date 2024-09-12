#ifndef CONSTANT_H
#define CONSTANT_H

#include <complex.h>

// Constants
extern const double pi               ;
extern const double complex im       ;
extern const double hbar             ;
extern const double kb               ;
extern const double au2ev            ;
extern const double lightspeed_c     ;
extern const double vac_permittivity ;
// Unit conversion constants
// Energy converters
extern const double au_2_wn          ;  // au to wavenumber
extern const double au_2_eV          ;  // au to electron Volt
extern const double au_2_kcalpmol    ;  // au to kcal per mole
extern const double au_2_kJpmol      ;  // au to kJ per mole
// Time converters
extern const double au_2_fs          ;  // au to femtosecond
extern const double au_2_ps          ;  // au to picosecond
// Mass converter
extern const double amu_2_au         ;  // amu to au
// Length converter
extern const double au_2_angstrom    ; // au to Angstrom
#endif // CONSTANT_H
