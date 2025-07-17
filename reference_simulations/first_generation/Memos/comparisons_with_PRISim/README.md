# Comparison with PRISim

The PRISim simulator was designed with the intent of simulating wide field,
high bandwidth interferometers specifically targeting 21 cm instruments. It has mainly been used
in the context of a delay-spectrum style analysis, though images have been made and are known to
be roughly correct. The data for this comparison can be found on the NRAO computing cluster at
`/lustre/aoc/projects/hera/djacobs/prisim_ref/`.

- In the frequency domain, PRISim simulations agree at 10^-5 with pyuvsim, with quantization type
errors at 10^-6.
- The difference shows a floor in the delay spectrum at 1e-6 which is above the expected Blackman
Harris floor at 1e-10
- In the image domain a phase offset at roughly the psf size and can be traced to differences in
uvw calculation.  Antenna positions agree to tolerance which suggests uvw difference is related
to phasing.  PRISim and pyvusim use different codes for phasing and time calculation. The PRISim
phasing code is not covered with analytic tests.
