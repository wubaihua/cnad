# cnad

## OC equilibrium density

For `msmodelname: "oc"`, `if_mqcequden: 2` calculates the unnormalized
electronic equilibrium matrix `exp(-beta H_e(R))`. With
`if_mqcequden: 1`, the matrix also includes the scalar nuclear kinetic and
potential energy: `exp(-beta (H_e(R) + H_n(R,P) I))`. The resulting
trajectory estimator is written to the `.mqcequden` output file. As in the
SBM implementation, this is not normalized by the matrix trace.

Set the temperature with the JSON key `temperature` (not `tempurature`).
The thermodynamic `beta` is distinct from the OC hopping parameter `beta_oc`.

The x86 CMake configuration can enable the OC equilibrium-density test with
`CNAD_BUILD_TESTS=ON`; run the `oc_mqcequden` test through CTest.

The OC model uses a sparse derivative transformation for adiabatic forces and
nonadiabatic couplings, while retaining the shared eigenstate matching logic.
The `oc_nacv` CTest compares it against the dense transformation for both
force types, including a near-degenerate pair of states.