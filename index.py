import sys

sys.path.insert(0, "src")
sys.path.insert(0, ".")

import panel as pn

from tabs import (
	hardy_weinberg_equilibrium,
	one_locus_drift,
	one_locus_drift_with_bottleneck,
	one_locus_drift_with_selection,
	one_locus_selection,
	two_loci_selection,
)

pn.extension("plotly")

app = pn.Tabs(
	hardy_weinberg_equilibrium.build_tab(),
	one_locus_selection.build_tab(),
	two_loci_selection.build_tab(),
	one_locus_drift.build_tab(),
	one_locus_drift_with_bottleneck.build_tab(),
	one_locus_drift_with_selection.build_tab(),
	margin=(15, 30),
)
app.servable(title="Evolutio")
