# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import plotly.plotly as py
from plotly.graph_objs import *

DoFs = [8958, 35460, 56286, 120264]

MPRGP = Scatter(
    x=DoFs,
    y=[0.9009429812, 6.372233868, 12.80960178, 24.31423187]
)
ICA = Scatter(
    x=DoFs,
    y=[3.460378885, 33.74056625, 79.34896088, 233.6007996]
)
data = Data([MPRGP, ICA])
plot_url = py.plot(data, filename='mprgp_ica_comp')
