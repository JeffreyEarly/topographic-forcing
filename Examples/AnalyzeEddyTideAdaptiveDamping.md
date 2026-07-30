# Adaptive damping is not local viscous dissipation

`AnalyzeEddyTideAdaptiveDamping` plots the exact local work of WaveVortexModel's adaptive spectral damping:

\[
P_d = uF_u + vF_v + wF_w + N^2\eta_{\mathrm{true}}F_\eta,
\qquad
\epsilon_{\mathrm{eff}}=-\rho_0P_d.
\]

The damping operator is negative definite in spectral space, so the spatial integral of \(P_d\) is non-positive. Its inverse transform is nonlocal, however, and \(P_d\) can have either sign at an individual point. Therefore \(\epsilon_{\mathrm{eff}}\) is an exact signed energy-removal field, not the non-negative viscous dissipation density used by [Shakespeare (2023)](https://doi.org/10.1175/JPO-D-23-0127.1).

A non-negative localization could be formed by applying the square root of the negative damping operator, \(\sqrt{-D}\), to each energy variable and then constructing its quadratic energy density. That field could preserve the integrated loss, but it would remain a derived localization of a nonlocal closure rather than physical viscosity. It is not included in this example.
