from __future__ import annotations

import numpy as np
import numpy.typing as npt
import plotly.express as px
import plotly.graph_objects as go

RNG = np.random.default_rng()

FloatArray = npt.NDArray[np.float64]


def _one_locus_selection_step(
	p: float | FloatArray, w: tuple[float, float, float]
) -> float | FloatArray:
	w_AA, w_Aa, w_aa = w
	q = 1.0 - p
	numerator = p * p * w_AA + p * q * w_Aa
	denominator = numerator + p * q * w_Aa + q * q * w_aa
	return numerator / denominator


def plot_one_locus_selection(
	T: int, f_A: float, w: tuple[float, float, float]
) -> go.Figure:
	w_AA, w_Aa, w_aa = w
	if T < 1:
		raise ValueError("T must be at least 1.")
	if not 0.0 <= f_A <= 1.0:
		raise ValueError("f_A must be between 0 and 1.")
	if min(w_AA, w_Aa, w_aa) < 0.0:
		raise ValueError("Fitness values must be non-negative.")
	if w_AA + w_Aa + w_aa <= 0.0:
		raise ValueError("At least one fitness value must be positive.")

	freq_A = np.zeros(T + 1)
	freq_A[0] = f_A
	for t in range(T):
		freq_A[t + 1] = _one_locus_selection_step(freq_A[t], (w_AA, w_Aa, w_aa))

	data = {"generation": np.arange(T + 1), "A": freq_A, "a": 1.0 - freq_A}
	fig = px.line(
		data,
		x="generation",
		y=["A", "a"],
		labels={"variable": "allele", "value": "frequency"},
	)

	fig.update_yaxes(range=[0, 1])
	fig.update_layout(
		xaxis_title="generation", yaxis_title="frequency", legend_title="allele"
	)

	return fig


def simulate_evolution(
	T: int,
	N: int,
	f_A: float,
	R: int = 1,
	P: int = 1,
	w: tuple[float, float, float] | None = None,
	m: float = 0.0,
	bottleneck: tuple[int, int, int] | None = None,
	deterministic: bool = False,
) -> FloatArray:
	if T < 1:
		raise ValueError("T must be at least 1.")
	if R < 1 or P < 1:
		raise ValueError("R and P must be at least 1.")
	if not 0.0 <= f_A <= 1.0:
		raise ValueError("f_A must be between 0 and 1.")
	if not 0.0 <= m <= 1.0:
		raise ValueError("m must be between 0 and 1.")

	if bottleneck is not None:
		generation_start, generation_end, N_bottleneck = bottleneck
		if not 1 <= generation_start <= generation_end <= T:
			raise ValueError("bottleneck window must satisfy 1 <= start <= end <= T.")
		N_schedule: int | FloatArray = np.full(T + 1, N, dtype=int)
		N_schedule[generation_start : generation_end + 1] = N_bottleneck
	else:
		N_schedule = np.full(T + 1, N, dtype=int)

	if N_schedule.shape != (T + 1,):
		raise ValueError("N must be a scalar or an array of length T + 1.")

	counts = np.zeros((T + 1, P, R), dtype=int)
	counts[0] = np.round(f_A * N_schedule[0])

	for t in range(T):
		freq = counts[t] / N_schedule[t]

		if w is not None:
			freq = _one_locus_selection_step(freq, w)

		if deterministic:
			counts[t + 1] = np.round(freq * N_schedule[t + 1]).astype(int)
			continue

		if P > 1 and m > 0.0:
			n_migrants = round(m * N_schedule[t + 1])
			other_deme_mean = (freq.sum(axis=0, keepdims=True) - freq) / (P - 1)
			other_deme_mean = np.clip(other_deme_mean, 0.0, 1.0)
			counts[t + 1] = RNG.binomial(
				N_schedule[t + 1] - n_migrants, np.clip(freq, 0.0, 1.0)
			) + RNG.binomial(n_migrants, other_deme_mean)
		else:
			counts[t + 1] = RNG.binomial(N_schedule[t + 1], np.clip(freq, 0.0, 1.0))

	return counts / N_schedule[:, None, None]


def plot_one_locus_evolution(
	sims: FloatArray,
	bottleneck: tuple[int, int, int] | None = None,
) -> go.Figure:
	if bottleneck is not None:
		generation_start, generation_end, _ = bottleneck
		if not 1 <= generation_start <= generation_end <= sims.shape[0] - 1:
			raise ValueError("bottleneck window must satisfy 1 <= start <= end <= T.")

	_, P, R = sims.shape

	fig = go.Figure()

	if bottleneck is not None:
		fig.add_vrect(
			x0=generation_start,
			x1=generation_end,
			fillcolor="yellow",
			opacity=0.2,
			line_width=0,
		)

	for deme in range(P):
		for replicate in range(R):
			trace_name = (
				f"deme {deme}, replicate {replicate}"
				if P > 1
				else f"replicate {replicate}"
			)
			fig.add_trace(
				go.Scatter(
					y=sims[:, deme, replicate],
					mode="lines",
					line={"color": "grey", "width": 1},
					opacity=0.5,
					name=trace_name,
					showlegend=False,
				)
			)

	fig.add_trace(
		go.Scatter(
			y=sims.mean(axis=(1, 2)),
			mode="lines",
			line={"color": "red", "width": 2},
			name="mean",
		)
	)

	fig.add_hline(y=sims[0, 0, 0], line={"color": "blue", "width": 2, "dash": "dash"})

	fig.update_yaxes(range=[0, 1])
	fig.update_layout(xaxis_title="generation", yaxis_title="frequency of A")

	return fig


def plot_heterozygosity(
	sims: FloatArray,
	bottleneck: tuple[int, int, int] | None = None,
) -> go.Figure:
	if bottleneck is not None:
		generation_start, generation_end, _ = bottleneck
		if not 1 <= generation_start <= generation_end <= sims.shape[0] - 1:
			raise ValueError("bottleneck window must satisfy 1 <= start <= end <= T.")

	_, P, R = sims.shape

	fig = go.Figure()

	if bottleneck is not None:
		fig.add_vrect(
			x0=generation_start,
			x1=generation_end,
			fillcolor="yellow",
			opacity=0.2,
			line_width=0,
		)

	heterozygosity = 2 * sims * (1 - sims)

	for deme in range(P):
		for replicate in range(R):
			trace_name = (
				f"deme {deme}, replicate {replicate}"
				if P > 1
				else f"replicate {replicate}"
			)
			fig.add_trace(
				go.Scatter(
					y=heterozygosity[:, deme, replicate],
					mode="lines",
					line={"color": "grey", "width": 1},
					opacity=0.5,
					name=trace_name,
					showlegend=False,
				)
			)

	fig.add_trace(
		go.Scatter(
			y=heterozygosity.mean(axis=(1, 2)),
			mode="lines",
			line={"color": "red", "width": 2},
			name="mean",
		)
	)

	fig.update_yaxes(range=[0, 0.5])
	fig.update_layout(xaxis_title="generation", yaxis_title="heterozygosity")

	return fig


def plot_replicate_variance(
	sims: FloatArray,
	bottleneck: tuple[int, int, int] | None = None,
) -> go.Figure:
	if bottleneck is not None:
		generation_start, generation_end, _ = bottleneck
		if not 1 <= generation_start <= generation_end <= sims.shape[0] - 1:
			raise ValueError("bottleneck window must satisfy 1 <= start <= end <= T.")

	_, P, R = sims.shape

	fig = go.Figure()

	if bottleneck is not None:
		fig.add_vrect(
			x0=generation_start,
			x1=generation_end,
			fillcolor="yellow",
			opacity=0.2,
			line_width=0,
		)

	if P == 1:
		variance = np.var(sims, axis=(1, 2))

		fig.add_trace(
			go.Scatter(
				y=variance,
				mode="lines",
				line={"color": "red", "width": 2},
				name="variance across replicates",
			)
		)

	else:
		variance = np.var(sims, axis=1)

		for replicate in range(R):
			fig.add_trace(
				go.Scatter(
					y=variance[:, replicate],
					mode="lines",
					line={"color": "grey", "width": 1},
					opacity=0.5,
					name=f"replicate {replicate}",
					showlegend=False,
				)
			)

		fig.add_trace(
			go.Scatter(
				y=variance.mean(axis=1),
				mode="lines",
				line={"color": "red", "width": 2},
				name="mean variance across demes",
			)
		)

	fig.update_layout(
		xaxis_title="generation", yaxis_title="variance in the frequency of A"
	)

	return fig
