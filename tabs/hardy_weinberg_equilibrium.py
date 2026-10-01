import panel as pn

from evolution_models import plot_hardy_weinberg_equilibrium


def build_tab():
	f_A_slider = pn.widgets.EditableFloatSlider(
		name="Frequency of A (f_A)",
		fixed_start=0.0,
		fixed_end=1.0,
		step=0.01,
		value=0.5,
	)
	AA_indicator = pn.indicators.Number(
		label="AA",
		value=0.25,
		format="{value:.4f}",
		default_color="red",
		font_size="24pt",
		title_size="16pt",
	)
	Aa_indicator = pn.indicators.Number(
		label="Aa",
		value=0.5,
		format="{value:.4f}",
		default_color="green",
		font_size="24pt",
		title_size="16pt",
	)
	aa_indicator = pn.indicators.Number(
		label="aa",
		value=0.25,
		format="{value:.4f}",
		default_color="blue",
		font_size="24pt",
		title_size="16pt",
	)
	geno_indicators = pn.Row(AA_indicator, Aa_indicator, aa_indicator)

	plot_output = pn.Column(plot_hardy_weinberg_equilibrium(ind_fA=0.5))

	def update_geno_indicators(event):
		plot_output.loading = True

		f = event.new

		AA_indicator.value = f * f
		Aa_indicator.value = 2 * f * (1 - f)
		aa_indicator.value = (1 - f) * (1 - f)

		try:
			plot_output[:] = [plot_hardy_weinberg_equilibrium(ind_fA=f)]
		finally:
			plot_output.loading = False

	f_A_slider.param.watch(update_geno_indicators, "value")

	layout = pn.Row(
		pn.Column(
			"### Parameters",
			f_A_slider,
			"#### Expected frequencies",
			geno_indicators,
			"### Qui-squared test",
			"coming soon...",
			width=450,
		),
		pn.Column(
			"### Hardy-Weinberg equilibrium",
			plot_output,
		),
	)
	return "Hardy-Weinberg equilibrium", layout
