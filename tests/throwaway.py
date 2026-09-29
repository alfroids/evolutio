import numpy as np
import panel as pn
import plotly.graph_objects as go

pn.extension("plotly")

slider = pn.widgets.FloatSlider(name="Frequency", start=0.5, end=5, step=0.5, value=1)


def make_plot(frequency):
    x = np.linspace(0, 10, 200)
    y = np.sin(frequency * x)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode="lines"))
    fig.update_layout(xaxis_title="x", yaxis_title="sin(frequency * x)")
    return fig


plot_pane = pn.pane.Plotly(pn.bind(make_plot, slider))

app = pn.Column(slider, plot_pane)
app.servable()
