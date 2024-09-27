#ifndef CONSTANT_H
#define CONSTANT_H

#include <complex.h>

// Constants
const double pi = 3.1415926535897932;

const double complex im = 1.0 * I;

const double hbar = 1.0;

const double kb = 1.3806E-5 / 4.35974;

const double au2ev = 27.211386245988;

const double lightspeed_c = 299792458.0;

const double vac_permittivity = 8.854187817E-12;

// Unit conversion constants
// Energy converters
const double au_2_wn = 219474.6313702;      // au to wavenumber
const double au_2_eV = 27.21138602;         // au to electron Volt
const double au_2_kcalpmol = 627.509474;    // au to kcal per mole
const double au_2_kJpmol = 2625.499638;     // au to kJ per mole

// Time converters
const double au_2_fs = 2.418884326505E-02;  // au to femtosecond
const double au_2_ps = 2.418884326505E-05;  // au to picosecond

// Mass converter
const double amu_2_au = 1.82289E3;          // amu to au

// Length converter
const double au_2_angstrom = 0.52917721067; // au to Angstrom

#endif // CONSTANT_H
