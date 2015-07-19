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

mprgp_0 = mprgp[0]
mprgp_opt = mprgp[-1]
cg_0 = cg[0]
cg_opt = cg[-1]
for i in range(0,len(mprgp)):
    mprgp[i] = abs(mprgp[i]-mprgp_opt)/abs(mprgp_0-mprgp_opt)
    cg[i] = abs(cg[i]-cg_opt)/abs(cg_0-cg_opt)

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
