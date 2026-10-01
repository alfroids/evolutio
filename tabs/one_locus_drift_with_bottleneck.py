import panel as pn

from evolution_models import (
	plot_heterozygosity,
	plot_one_locus_evolution,
	plot_replicate_variance,
	simulate_evolution,
)


def build_tab():
	T_slider = pn.widgets.EditableIntSlider(
		name="Generations (T)", fixed_start=1, end=500, step=1, value=100
	)
	T_bottleneck_slider = pn.widgets.IntRangeSlider(
		name="Bottleneck interval (start .. end)",
		start=1,
		end=500,
		step=1,
		value=(30, 60),
	)
	N_slider = pn.widgets.EditableIntSlider(
		name="Population size (N)", fixed_start=1, end=500, step=1, value=50
	)
	N_bottleneck_slider = pn.widgets.EditableIntSlider(
		name="Population size during bottleneck",
		fixed_start=1,
		end=100,
		step=1,
		value=5,
	)
	f_A_slider = pn.widgets.EditableFloatSlider(
		name="Starting frequency of A (f_A)", start=0.0, end=1.0, step=0.01, value=0.25
	)
	R_slider = pn.widgets.EditableIntSlider(
		name="Number of replicates", fixed_start=1, end=50, step=1, value=20
	)
	update_button = pn.widgets.Button(name="Update plot", button_type="primary")

	def make_plots(T, N, f_A, R, bottleneck):
		try:
			sims = simulate_evolution(T=T, N=N, f_A=f_A, R=R, bottleneck=bottleneck)
			fig1 = plot_one_locus_evolution(sims, bottleneck=bottleneck)
			fig2 = plot_heterozygosity(sims, bottleneck=bottleneck)
			fig3 = plot_replicate_variance(sims, bottleneck=bottleneck)
		except ValueError as error:
			return pn.pane.Markdown(f"**Invalid parameters:** {error}")

		return fig1, fig2, fig3

	plot_output = pn.Column(
		*make_plots(
			T=T_slider.value,
			N=N_slider.value,
			f_A=f_A_slider.value,
			R=R_slider.value,
			bottleneck=(*T_bottleneck_slider.value, N_bottleneck_slider.value),
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
					bottleneck=(*T_bottleneck_slider.value, N_bottleneck_slider.value),
				)
			]
		finally:
			plot_output.loading = False

	update_button.on_click(on_update)

	layout = pn.Row(
		pn.Column(
			"### Parameters",
			T_slider,
			T_bottleneck_slider,
			N_slider,
			N_bottleneck_slider,
			f_A_slider,
			R_slider,
			update_button,
			width=450,
		),
		pn.Column(
			"### One-locus drift with bottleneck",
			plot_output,
		),
	)
	return "One-locus drift with bottleneck", layout
