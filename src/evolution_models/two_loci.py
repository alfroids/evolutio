from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np
import numpy.typing as npt
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go

from ._generated_symbolic import delta_fA_fn, delta_fB_fn, jacobian_fn
from .one_locus import RNG, _one_locus_selection_step

FloatArray = npt.NDArray[np.float64]
StabilityLabel = Literal["stable", "unstable", "saddle"]


def delta_field(
	f_A: float | FloatArray, f_B: float | FloatArray, w: FloatArray
) -> tuple[float | FloatArray, float | FloatArray]:
	w_flat = w.flatten()
	return delta_fA_fn(f_A, f_B, *w_flat), delta_fB_fn(f_A, f_B, *w_flat)


def _newton_solve(
	x0: FloatArray,
	w_flat: FloatArray,
	max_iterations: int = 50,
	tol: float = 1e-10,
) -> tuple[FloatArray, bool]:
	x = np.asarray(x0, dtype=float)
	for _ in range(max_iterations):
		residual = np.array(
			[delta_fA_fn(x[0], x[1], *w_flat), delta_fB_fn(x[0], x[1], *w_flat)]
		)
		jacobian = jacobian_fn(x[0], x[1], *w_flat)
		try:
			step = np.linalg.solve(jacobian, residual)
		except np.linalg.LinAlgError:
			return x, False
		x = x - step
		if np.linalg.norm(step) < tol:
			return x, True
	return x, False


def find_equilibria(
	w: FloatArray,
	seed_points: Sequence[tuple[float, float]] = (),
	n_backup_seeds: int = 5,
	dedupe_tol: float = 1e-4,
) -> FloatArray:
	w_flat = w.flatten()

	corners = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]
	backup_grid = np.linspace(0.05, 0.95, n_backup_seeds)
	candidates = (
		list(seed_points) + corners + [(x, y) for x in backup_grid for y in backup_grid]
	)

	found_points = []
	for x0, y0 in candidates:
		solution, converged = _newton_solve(np.array([x0, y0]), w_flat)
		if not converged:
			continue
		if not (-1e-6 <= solution[0] <= 1 + 1e-6 and -1e-6 <= solution[1] <= 1 + 1e-6):
			continue
		solution = np.clip(solution, 0.0, 1.0)
		residual = np.array(
			[
				delta_fA_fn(solution[0], solution[1], *w_flat),
				delta_fB_fn(solution[0], solution[1], *w_flat),
			]
		)
		if not np.allclose(residual, 0.0, atol=1e-7):
			continue
		found_points.append(solution)

	unique_points: list = []
	for point in found_points:
		if not any(np.linalg.norm(point - u) < dedupe_tol for u in unique_points):
			unique_points.append(point)
	return np.array(unique_points)


def classify_equilibrium(point: FloatArray, w: FloatArray) -> StabilityLabel:
	jacobian = jacobian_fn(point[0], point[1], *w.flatten())
	eigenvalue_moduli = np.abs(np.linalg.eigvals(np.eye(2) + jacobian))
	if np.all(eigenvalue_moduli < 1 - 1e-6):
		return "stable"
	if np.all(eigenvalue_moduli > 1 + 1e-6):
		return "unstable"
	return "saddle"


_EQUILIBRIUM_MARKER_STYLE = {
	"stable": {"symbol": "circle", "color": "#2ecc71"},
	"unstable": {"symbol": "x", "color": "#e74c3c"},
	"saddle": {"symbol": "diamond", "color": "#f1c40f"},
}


def _build_two_loci_phase_portrait(
	w: FloatArray,
	heatmap_resolution: int,
	streamline_resolution: int,
	previous_equilibria: FloatArray | None,
) -> tuple[go.Figure, FloatArray]:
	heatmap_grid = np.linspace(0.0, 1.0, heatmap_resolution)
	f_A_grid, f_B_grid = np.meshgrid(heatmap_grid, heatmap_grid)
	delta_A_grid, delta_B_grid = delta_field(f_A_grid, f_B_grid, w)
	magnitude = np.sqrt(delta_A_grid**2 + delta_B_grid**2)

	fig = px.imshow(
		magnitude,
		x=heatmap_grid,
		y=heatmap_grid,
		origin="lower",
		color_continuous_scale="Viridis",
	)

	stream_grid = np.linspace(0.01, 0.99, streamline_resolution)
	f_A_stream, f_B_stream = np.meshgrid(stream_grid, stream_grid)
	delta_A_stream, delta_B_stream = delta_field(f_A_stream, f_B_stream, w)

	stream_fig = ff.create_streamline(
		stream_grid,
		stream_grid,
		delta_A_stream,
		delta_B_stream,
		arrow_scale=0.03,
		line={"color": "white", "width": 1},
	)
	for trace in stream_fig.data:
		trace.showlegend = False
		fig.add_trace(trace)

	seed_points = (
		[tuple(point) for point in previous_equilibria]
		if previous_equilibria is not None
		else ()
	)
	equilibria = find_equilibria(w, seed_points=seed_points)
	stability_labels = np.array(
		[classify_equilibrium(point, w) for point in equilibria]
	)

	for label, marker_spec in _EQUILIBRIUM_MARKER_STYLE.items():
		mask = stability_labels == label
		if not np.any(mask):
			continue
		fig.add_trace(
			go.Scatter(
				x=equilibria[mask, 0],
				y=equilibria[mask, 1],
				mode="markers",
				marker={"size": 12, **marker_spec},
				name=label,
			)
		)

	fig.update_layout(
		xaxis_title="frequency of A",
		yaxis_title="frequency of B",
		legend={
			"title": "equilibrium",
			"orientation": "h",
			"yanchor": "bottom",
			"y": -0.25,
			"xanchor": "center",
			"x": 0.5,
		},
		coloraxis_colorbar={"title": "magnitude of change"},
	)
	return fig, equilibria


def plot_two_loci_selection(
	w: FloatArray,
	heatmap_resolution: int = 201,
	streamline_resolution: int = 40,
	previous_equilibria: FloatArray | None = None,
) -> tuple[go.Figure, FloatArray]:
	w_internal = np.flip(w, axis=1).T
	return _build_two_loci_phase_portrait(
		w_internal, heatmap_resolution, streamline_resolution, previous_equilibria
	)


def simulate_two_loci_drift_with_selection(
	T: int, N: int, f_A: float, f_B: float, w: FloatArray
) -> FloatArray:
	if T < 1:
		raise ValueError("T must be at least 1.")
	if N < 1:
		raise ValueError("N must be at least 1.")
	if not (0.0 <= f_A <= 1.0 and 0.0 <= f_B <= 1.0):
		raise ValueError("f_A and f_B must be between 0 and 1.")

	counts = np.zeros((T + 1, 2))
	counts[0] = np.round([f_A * N, f_B * N])

	for t in range(T):
		p_A = counts[t, 0] / N
		q_A = 1.0 - p_A
		p_B = counts[t, 1] / N
		q_B = 1.0 - p_B

		w_AA = w[0, 2] * p_B * p_B + w[1, 2] * 2 * p_B * q_B + w[2, 2] * q_B * q_B
		w_Aa = w[0, 1] * p_B * p_B + w[1, 1] * 2 * p_B * q_B + w[2, 1] * q_B * q_B
		w_aa = w[0, 0] * p_B * p_B + w[1, 0] * 2 * p_B * q_B + w[2, 0] * q_B * q_B

		w_BB = w[0, 2] * p_A * p_A + w[0, 1] * 2 * p_A * q_A + w[0, 0] * q_A * q_A
		w_Bb = w[1, 2] * p_A * p_A + w[1, 1] * 2 * p_A * q_A + w[1, 0] * q_A * q_A
		w_bb = w[2, 2] * p_A * p_A + w[2, 1] * 2 * p_A * q_A + w[2, 0] * q_A * q_A

		next_p_A = _one_locus_selection_step(p_A, (w_AA, w_Aa, w_aa))
		next_p_B = _one_locus_selection_step(p_B, (w_BB, w_Bb, w_bb))

		counts[t + 1] = RNG.binomial(N, [next_p_A, next_p_B])

	return counts / N


def plot_two_loci_selection_with_drift(
	T: int,
	N: int,
	f_A: float,
	f_B: float,
	w: FloatArray,
	heatmap_resolution: int = 201,
	streamline_resolution: int = 40,
	previous_equilibria: FloatArray | None = None,
) -> tuple[go.Figure, FloatArray]:
	trajectory = simulate_two_loci_drift_with_selection(T, N, f_A, f_B, w)
	w_internal = np.flip(w, axis=1).T

	fig, equilibria = _build_two_loci_phase_portrait(
		w_internal, heatmap_resolution, streamline_resolution, previous_equilibria
	)

	fig.add_trace(
		go.Scatter(
			x=trajectory[:, 0],
			y=trajectory[:, 1],
			mode="lines",
			line={"color": "red", "width": 4},
			opacity=0.5,
			name="trajectory",
			legend="legend2",
		)
	)
	fig.add_trace(
		go.Scatter(
			x=[trajectory[0, 0]],
			y=[trajectory[0, 1]],
			mode="markers",
			marker={"symbol": "circle", "size": 10, "color": "red"},
			name="start",
			legend="legend2",
		)
	)
	fig.add_trace(
		go.Scatter(
			x=[trajectory[-1, 0]],
			y=[trajectory[-1, 1]],
			mode="markers",
			marker={"symbol": "square", "size": 10, "color": "red"},
			name="end",
			legend="legend2",
		)
	)
	fig.update_layout(
		legend2={
			"orientation": "h",
			"yanchor": "bottom",
			"y": -0.4,
			"xanchor": "center",
			"x": 0.5,
		}
	)
	return fig, equilibria
