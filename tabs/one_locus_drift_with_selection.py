import panel as pn

from evolution_models import (
	plot_heterozygosity,
	plot_one_locus_evolution,
	plot_replicate_variance,
	simulate_evolution,
)


def build_tab():
	T_slider = pn.widgets.EditableIntSlider(
		name="Generations (T)", start=1, end=500, step=1, value=100
	)
	N_slider = pn.widgets.EditableIntSlider(
		name="Population size (N)", start=1, end=500, step=1, value=50
	)
	f_A_slider = pn.widgets.EditableFloatSlider(
		name="Starting frequency of A (f_A)", start=0.0, end=1.0, step=0.01, value=0.25
	)
	R_slider = pn.widgets.EditableIntSlider(
		name="Number of replicates", start=1, end=50, step=1, value=20
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

	def make_plots(T, N, f_A, R, w_AA, w_Aa, w_aa):
		try:
			sims = simulate_evolution(T=T, N=N, f_A=f_A, R=R, w=(w_AA, w_Aa, w_aa))
			fig1 = plot_one_locus_evolution(sims)
			fig2 = plot_heterozygosity(sims)
			fig3 = plot_replicate_variance(sims)
		except ValueError as error:
			return pn.pane.Markdown(f"**Invalid parameters:** {error}")

		return fig1, fig2, fig3

	plot_output = pn.Column(
		*make_plots(
			T=T_slider.value,
			N=N_slider.value,
			f_A=f_A_slider.value,
			R=R_slider.value,
			w_AA=w_AA_input.value,
			w_Aa=w_Aa_input.value,
			w_aa=w_aa_input.value,
		)
	)

	def on_update(event):
		plot_output.loading = True
		try:
			plot_output[:] = [
				*make_plots(
					T=T_slider.value,
					N=N_slider.value,
					f_A=f_A_slider.value,
					R=R_slider.value,
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
			N_slider,
			f_A_slider,
			R_slider,
			w_grid,
			update_button,
			width=450,
		),
		pn.Column(
			"### One-locus drift with selection",
			plot_output,
		),
	)
	return "One-locus drift with selection", layout
