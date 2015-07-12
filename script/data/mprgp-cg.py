# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import plotly.plotly as py
from plotly.graph_objs import *

def loadData(filename):
    f = open(filename).read().split(",")
    data = []
    for d in f:
        data.append(float(d))
    return data

# main
mprgp = loadData("mprgp.txt")[0:100]
cg = loadData("cg.txt")[0:100]

trace1 = Scatter(
    x=range(0,len(mprgp)),
    y=mprgp
)
trace2 = Scatter(
    x=range(0,len(mprgp)),
    y=cg
)
data = Data([trace1, trace2])
plot_url = py.plot(data, filename='mprgp-cg')
