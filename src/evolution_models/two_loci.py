"""
Two-locus selection model: the deterministic phase portrait (change field
and equilibria), and a stochastic drift trajectory overlaid on it.

The one-generation frequency-change functions and their Jacobian are
imported from _generated_symbolic, which was produced once by
tools/derive_symbolic_model.py from a symbolic (sympy) derivation. Neither
sympy nor scipy is a dependency of this module - equilibria are found with
a small Newton's-method solver using the same analytic Jacobian that
sympy derived, rather than scipy.optimize.fsolve.
"""

from __future__ import annotations

from typing import Literal, Optional, Sequence, Tuple

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
) -> Tuple[float | FloatArray, float | FloatArray]:
    """One-generation change in two-locus allele frequencies.

    f_A, f_B: current frequency of allele A and allele B (scalars or arrays
        of matching shape).
    w: 3x3 fitness matrix in the model's internal convention - rows indexed
        by genotype at locus A (AA, Aa, aa), columns by genotype at locus B
        (BB, Bb, bb).
    Returns (delta_f_A, delta_f_B), the change in each frequency.
    """
    w_flat = w.flatten()
    return delta_fA_fn(f_A, f_B, *w_flat), delta_fB_fn(f_A, f_B, *w_flat)


def _newton_solve(
    x0: FloatArray,
    w_flat: FloatArray,
    max_iterations: int = 50,
    tol: float = 1e-10,
) -> Tuple[FloatArray, bool]:
    """Newton's method root finder for (delta_fA, delta_fB) = 0, using the
    analytic Jacobian.

    x0: starting point (f_A, f_B).
    w_flat: flattened 3x3 fitness matrix (9 values), model convention.
    max_iterations: iteration cap; the search gives up (reports failure)
        rather than looping indefinitely on a starting point that does not
        converge.
    tol: the step is considered converged once its L2 norm drops below
        this value.
    Returns (solution, converged).
    """
    x = np.asarray(x0, dtype=float)
    for _ in range(max_iterations):
        residual = np.array(
            [delta_fA_fn(x[0], x[1], *w_flat), delta_fB_fn(x[0], x[1], *w_flat)]
        )
        jacobian = jacobian_fn(x[0], x[1], *w_flat)
        try:
            step = np.linalg.solve(jacobian, residual)
        except np.linalg.LinAlgError:
            return x, False  # singular Jacobian at this point - give up
        x = x - step
        if np.linalg.norm(step) < tol:
            return x, True
    return x, False


def find_equilibria(
    w: FloatArray,
    seed_points: Sequence[Tuple[float, float]] = (),
    n_backup_seeds: int = 5,
    dedupe_tol: float = 1e-4,
) -> FloatArray:
    """Root-find equilibria of the two-locus selection model for matrix w.

    w: 3x3 fitness matrix in the model's internal convention (see
        delta_field).
    seed_points: previously found equilibria, used to warm-start the solver
        when w has only changed slightly (e.g. after a slider drag), so
        convergence is fast instead of a cold search from scratch.
    n_backup_seeds: a backup grid of n_backup_seeds^2 starting points is
        always searched too, so equilibria that newly appear as w crosses a
        bifurcation are still found even without a useful warm start.
    dedupe_tol: minimum distance between two roots for them to be treated as
        distinct equilibria.
    Returns an array of shape (n_equilibria, 2) with unique equilibria in
    [0, 1] x [0, 1].
    """
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
    """Stability of an equilibrium of the two-locus selection model.

    point: (f_A, f_B) coordinates of the equilibrium.
    w: 3x3 fitness matrix in the model's internal convention, used to
        evaluate the local Jacobian.
    Classifies stability from the linearized discrete map f' = f + delta_f(f):
    "stable" if all eigenvalue moduli of (I + J) are below 1, "unstable" if
    all are above 1, otherwise "saddle".
    """
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
    previous_equilibria: Optional[FloatArray],
) -> Tuple[go.Figure, FloatArray]:
    """Heatmap + streamlines + equilibria for the two-locus selection field.

    w: 3x3 fitness matrix in the model's internal convention.
    heatmap_resolution: number of grid points per axis for the background
        heatmap of the magnitude of change.
    streamline_resolution: number of grid points per axis for the vector
        field used to draw direction-of-change streamlines.
    previous_equilibria: equilibria from a previous call, passed through to
        find_equilibria as warm-start seeds.
    Returns (fig, equilibria), with frequency of A on the x axis and
    frequency of B on the y axis.
    """
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

    # Streamlines are drawn on a coarser grid and avoid the exact
    # boundaries, where the vector field is degenerate.
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
    previous_equilibria: Optional[FloatArray] = None,
) -> Tuple[go.Figure, FloatArray]:
    """Deterministic two-locus phase portrait: change field and equilibria.

    w: 3x3 fitness matrix in table layout - rows indexed by genotype at
        locus B (BB, Bb, bb) and columns by genotype at locus A in reverse
        order (aa, Aa, AA), i.e. the layout produced by writing out a table
        with A genotypes decreasing left to right and B genotypes
        decreasing top to bottom. Converted internally to the model's
        convention before use.
    heatmap_resolution: grid density for the background heatmap.
    streamline_resolution: grid density for the streamlines.
    previous_equilibria: equilibria from a previous call, used to warm-start
        the equilibrium search.
    Returns (fig, equilibria), with frequency of A on the x axis and
    frequency of B on the y axis.
    """
    w_internal = np.flip(w, axis=1).T
    return _build_two_loci_phase_portrait(
        w_internal, heatmap_resolution, streamline_resolution, previous_equilibria
    )


def simulate_two_loci_drift_with_selection(
    T: int, N: int, f_A: float, f_B: float, w: FloatArray
) -> FloatArray:
    """Stochastic simulation of two-locus selection and drift.

    T: number of generations to simulate.
    N: population size (a single effective size applies to both loci, since
        they are assumed to be in linkage equilibrium).
    f_A, f_B: starting frequency of allele A and allele B.
    w: 3x3 fitness matrix in the same table layout as plot_two_loci_selection
        - rows indexed by genotype at locus B (BB, Bb, bb), columns by
        genotype at locus A in reverse order (aa, Aa, AA).
    Returns an array of shape (T + 1, 2) with columns (f_A, f_B) over time.

    At each generation, genotype fitness at one locus is obtained by
    averaging over the current genotype frequencies at the other locus -
    this is what linkage equilibrium buys us, since haplotype frequencies
    never need to be tracked directly - selection is then applied at each
    locus independently, and the next generation is sampled by binomial
    draws.
    """
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

        # Genotype fitness at locus A, averaged over locus B's genotype
        # frequencies (and symmetrically for locus B).
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
    previous_equilibria: Optional[FloatArray] = None,
) -> Tuple[go.Figure, FloatArray]:
    """Two-locus phase portrait with one stochastic trajectory overlaid.

    T, N, f_A, f_B: passed to simulate_two_loci_drift_with_selection.
    w: 3x3 fitness matrix in table layout, as described in
        simulate_two_loci_drift_with_selection and plot_two_loci_selection.
    heatmap_resolution, streamline_resolution, previous_equilibria: passed
        to the phase-portrait construction (see plot_two_loci_selection).
    Returns (fig, equilibria), with the simulated (f_A, f_B) trajectory
    drawn as a thick, semi-transparent red line on top of the deterministic
    phase portrait, its starting point marked with a circle and its ending
    point with a square. The trajectory's legend entries sit on their own
    row below the equilibrium legend.
    """
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
