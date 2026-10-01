import numpy as np
import panel as pn

from evolution_models import plot_two_loci_selection


def build_tab():
	# A_genos = ("aa", "Aa", "AA")
	# B_genos = ("BB", "Bb", "bb")
	w_inputs = [
		[
			pn.widgets.FloatInput(
				name="w_aaBB", value=1.0, start=0.0, step=0.05, width=95
			),
			pn.widgets.FloatInput(
				name="w_AaBB", value=1.1, start=0.0, step=0.05, width=95
			),
			pn.widgets.FloatInput(
				name="w_AABB", value=1.0, start=0.0, step=0.05, width=95
			),
		],
		[
			pn.widgets.FloatInput(
				name="w_aaBb", value=0.9, start=0.0, step=0.05, width=95
			),
			pn.widgets.FloatInput(
				name="w_AaBb", value=0.8, start=0.0, step=0.05, width=95
			),
			pn.widgets.FloatInput(
				name="w_AABb", value=1.0, start=0.0, step=0.05, width=95
			),
		],
		[
			pn.widgets.FloatInput(
				name="w_aabb", value=0.9, start=0.0, step=0.05, width=95
			),
			pn.widgets.FloatInput(
				name="w_Aabb", value=1.0, start=0.0, step=0.05, width=95
			),
			pn.widgets.FloatInput(
				name="w_AAbb", value=1.1, start=0.0, step=0.05, width=95
			),
		],
	]
	w_grid = pn.GridBox(*[w for row in w_inputs for w in row], ncols=3)
	update_button = pn.widgets.Button(name="Update plot", button_type="primary")

	def fix_float_display(event):
		for row in w_inputs:
			for w in row:
				w.value = round(w.value, 6)

	for row in w_inputs:
		for w in row:
			w.param.watch(fix_float_display, "value")

	def make_plot(w):
		try:
			fig, _ = plot_two_loci_selection(w=w)
		except ValueError as error:
			return pn.pane.Markdown(f"**Invalid parameters:** {error}")

		return fig

	plot_output = pn.Column(
		make_plot(w=np.array([[cell.value for cell in row] for row in w_inputs]))
	)

	def on_update(event):
		plot_output.loading = True
		try:
			plot_output[:] = [
				make_plot(
					w=np.array([[cell.value for cell in row] for row in w_inputs])
				)
			]
		finally:
			plot_output.loading = False

	update_button.on_click(on_update)

	layout = pn.Row(
		pn.Column(
			"### Parameters",
			"**Fitness matrix (AxB)**",
			w_grid,
			update_button,
			width=450,
		),
		pn.Column(
			"### Two-loci selection",
			plot_output,
		),
	)
	return "Two-loci selection", layout
