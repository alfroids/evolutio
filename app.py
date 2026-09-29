"""
Interactive web app for the evolutionary models in evolution_models.

Each tab is built by its own build_tab() function in tabs/, which creates
all of its widgets, panes, and callbacks as local variables and returns
(title, layout). This file only wires those tabs together - it holds no
widget or pane instances of its own, so there is no way for two tabs to
accidentally end up sharing the same underlying Bokeh model (the bug
behind a tab appearing during the loading skeleton and then vanishing
once rendering completes).
"""

import sys

# Plain relative paths, not something built from __file__: under panel
# convert's pyodide-worker target, this script is exec'd from a string
# rather than run as a real file, so __file__ is never defined there. The
# --resources bundle unpacks preserving this project's layout relative to
# the working directory, so relative paths work in both that environment
# and a normal `panel serve` run from the project root.
sys.path.insert(0, "src")
sys.path.insert(0, ".")

import panel as pn

from tabs import (
	one_locus_drift,
	one_locus_drift_with_bottleneck,
	one_locus_drift_with_selection,
	one_locus_selection,
	two_loci_selection,
)

pn.extension("plotly")

app = pn.Tabs(
	one_locus_selection.build_tab(),
	two_loci_selection.build_tab(),
	one_locus_drift.build_tab(),
	one_locus_drift_with_bottleneck.build_tab(),
	one_locus_drift_with_selection.build_tab(),
)
app.servable(title="Evolutio")
