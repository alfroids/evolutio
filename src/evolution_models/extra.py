import numpy as np
import plotly.graph_objects as go


def plot_hardy_weinberg_equilibrium(
	ind_fA: float | None = None, resolution: int = 201
) -> go.Figure:
	fA = np.linspace(0.0, 1.0, resolution)

	fig = go.Figure()

	fig.add_trace(
		go.Scatter(
			x=fA,
			y=fA * fA,
			mode="lines",
			line={"color": "red", "width": 2},
			name="AA genotype",
		)
	)

	fig.add_trace(
		go.Scatter(
			x=fA,
			y=2 * fA * (1 - fA),
			mode="lines",
			line={"color": "green", "width": 2},
			name="Aa genotype",
		)
	)

	fig.add_trace(
		go.Scatter(
			x=fA,
			y=(1 - fA) * (1 - fA),
			mode="lines",
			line={"color": "blue", "width": 2},
			name="aa genotype",
		)
	)

	if ind_fA is not None:
		fig.add_vline(
			x=ind_fA,
			line_width=2,
			line_color="black",
			line_dash="dash",
		)

	fig.update_layout(
		xaxis_title="A allele frequency", yaxis_title="genotype frequency"
	)

	return fig
