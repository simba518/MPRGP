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
general_mprgp = loadData("dec_solving_time.txt")
decoupled_mprgp = loadData("gen_solving_time.txt")
print general_mprgp[0:10]
print decoupled_mprgp[0:10]

trace1 = Scatter(
    x=range(0,len(general_mprgp)),
    y=general_mprgp
)
trace2 = Scatter(
    x=range(0,len(decoupled_mprgp)),
    y=decoupled_mprgp
)
data = Data([trace1, trace2])
plot_url = py.plot(data, filename='decoupled_general_comp')
