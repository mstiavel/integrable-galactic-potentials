Intermediate files for "Integrable galactic potentials beyond Stäckel: an analytic search"
==========================================================================================

These are the pickles (.pkl) and arrays (.npy) that the scripts in 02_axis_direct_method,
03_scalefree_direct_method and 04_transition_2d_direct_method write and later scripts read.
They are provided so the slow steps need not be rerun. They are NOT needed to use the
repository: tools/generate_intermediates.py regenerates all of them.

Layout: each file sits under the directory of the script that wrote it. Unpack at the
repository root and the files land where the scripts expect them:

    tar xzf intermediates-<version>.tar.gz

MANIFEST.txt lists, for every generating command, the exit status, wall time, and the SHA-256
of each output, plus the Python / sympy / numpy / scipy versions used. Check integrity with

    shasum -a 256 -c SHA256SUMS

Caveats
* Security: unpickling can execute arbitrary code. Load these files only if you trust their
  source (this release of this repository) and the checksums match; otherwise regenerate.
* Compatibility: pickles of sympy objects are tied to the sympy version that wrote them
  (see MANIFEST.txt). With a different sympy version, regenerate rather than load.
* Not included: 03_scalefree_direct_method/fast_deg6_d0_mu1.pkl (sf_fast.py 6 0 1 exceeded the
  3 h limit). first_deg6_d0.pkl from sf_first.py 6 0 1, which is included, gives the same
  decisive condition.
* 02_axis_direct_method/so_num.py (numerical solves, 10-30 min, parameters chosen
  interactively) and the 05 Birkhoff logs are not included.
