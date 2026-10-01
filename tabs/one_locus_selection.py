import panel as pn

from evolution_models import plot_one_locus_selection


def build_tab():
	T_slider = pn.widgets.EditableIntSlider(
		name="Generations (T)", fixed_start=1, end=500, step=1, value=100
	)
	f_A_slider = pn.widgets.EditableFloatSlider(
		name="Starting frequency of A (f_A)",
		fixed_start=0.0,
		fixed_end=1.0,
		step=0.01,
		value=0.25,
	)
	w_AA_input = pn.widgets.FloatInput(name="w_AA", value=1.0, start=0.0, width=95)
	w_Aa_input = pn.widgets.FloatInput(name="w_Aa", value=0.9, start=0.0, width=95)
	w_aa_input = pn.widgets.FloatInput(name="w_aa", value=0.8, start=0.0, width=95)
	w_grid = pn.GridBox(w_AA_input, w_Aa_input, w_aa_input, ncols=1)
	update_button = pn.widgets.Button(name="Update plot", button_type="primary")

	def fix_float_display(event):
		w_AA_input.value = round(w_AA_input.value, 6)
		w_Aa_input.value = round(w_Aa_input.value, 6)
		w_aa_input.value = round(w_aa_input.value, 6)

	w_AA_input.param.watch(fix_float_display, "value")
	w_Aa_input.param.watch(fix_float_display, "value")
	w_aa_input.param.watch(fix_float_display, "value")

	def make_plot(T, f_A, w_AA, w_Aa, w_aa):
		try:
			fig = plot_one_locus_selection(T=T, f_A=f_A, w=(w_AA, w_Aa, w_aa))
		except ValueError as error:
			return pn.pane.Markdown(f"**Invalid parameters:** {error}")

		return fig

	plot_output = pn.Column(
		make_plot(
			T=T_slider.value,
			f_A=f_A_slider.value,
			w_AA=w_AA_input.value,
			w_Aa=w_Aa_input.value,
			w_aa=w_aa_input.value,
		)
	)

	def on_update(event):
		plot_output.loading = True
		try:
			plot_output[:] = [
				make_plot(
					T=T_slider.value,
					f_A=f_A_slider.value,
					w_AA=w_AA_input.value,
					w_Aa=w_Aa_input.value,
					w_aa=w_aa_input.value,
				)
			]
		finally:
			plot_output.loading = False

	update_button.on_click(on_update)

	layout = pn.Row(
		pn.Column(
			"### Parameters",
			T_slider,
			f_A_slider,
			"**Fitness (w_AA, w_Aa, w_aa)**",
			w_grid,
			update_button,
			width=450,
		),
		pn.Column(
			"### One-locus selection",
			plot_output,
		),
	)
	return "One-locus selection", layout
