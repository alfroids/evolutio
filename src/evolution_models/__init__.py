"""
Models of allele-frequency change under selection, genetic drift,
migration, and population-size bottlenecks - split into a one-locus
module (drift, selection, migration, bottleneck, freely combinable) and a
two-loci module (deterministic phase portrait plus a stochastic drift
trajectory).

This package has no sympy or scipy dependency: the two-locus model's
symbolic derivation was done once by tools/derive_symbolic_model.py and
its result is hardcoded in _generated_symbolic.py, and equilibria are
found with a small Newton's-method solver instead of scipy.optimize.
"""

from .one_locus import (
	plot_heterozygosity,
	plot_one_locus_evolution,
	plot_one_locus_selection,
	plot_replicate_variance,
	simulate_evolution,
)
from .two_loci import (
	classify_equilibrium,
	delta_field,
	find_equilibria,
	plot_two_loci_selection,
	plot_two_loci_selection_with_drift,
	simulate_two_loci_drift_with_selection,
)

__all__ = [
	"classify_equilibrium",
	"delta_field",
	"find_equilibria",
	"plot_heterozygosity",
	"plot_one_locus_evolution",
	"plot_one_locus_selection",
	"plot_replicate_variance",
	"plot_two_loci_selection",
	"plot_two_loci_selection_with_drift",
	"simulate_evolution",
	"simulate_two_loci_drift_with_selection",
]
