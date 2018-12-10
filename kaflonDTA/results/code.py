# Python code for dominique

# Extraction
species = []
amino = []
da = []
ami = {' H ': 0,
       ' F ': 1,
       ' Y ': 2,
       ' C ': 3,
       ' N ': 4,
       ' D ': 5,
       ' E ': 6,
       ' Q ': 7,
       ' K ': 8}
X = 440
minv = 7
# replace the folder by other for Fungi etc
datatot = pd.concat([pd.read_csv('parameterFitsBacteria/nonlinfit' + str(i) + '.txt', header=None) for i in range(minv, 16)])
datatot = datatot[datatot[4] < 0.00999]
datatot = datatot[datatot[5] < 0.00999]
data = datatot.drop([4, 5], axis=1)
data = np.array(data)
arr = np.zeros((X, len(ami), 16 - minv))
previ = -1
j = 0
for line in data:
  for i in range(X):
    if line[0] == i + 1:
      for k, v in ami.iteritems():
        if line[1] == k:
          arr[i, v, j] = man(line[2], line[3])  # eucli(line[2],line[3])
          break
      break
  if i - previ < 0:
    j += 1
  previ = i
arr = np.ma.masked_equal(arr, 0)
arr = arr.mean(2)

# replace the number by the one corresponding to the required AA.
amiC = arr[:, 3]
regC = np.array(amiC[np.logical_not(np.ma.getmask(amiC))])
indC = np.argsort(regC)
regC = regC[indC]
variationC = (np.array(arr.var(1))[np.logical_not(np.ma.getmask(amiC))])[indC]


# bokeh plots
from math import pi
from bokeh.plotting import figure, show, output_file, output_notebook
output_notebook()
np.savetxt("CBacteriaMean.csv", np.vstack((indC, regC, variationC)).T,
           delimiter=",")

TOOLS = "pan,wheel_zoom,box_zoom,reset,save"

p = figure(tools=TOOLS, plot_width=1000, title="C aa and variation across aa for Bacterias")
p.xaxis.major_label_orientation = pi / 4
p.grid.grid_line_alpha = 0.3

p.vbar(x=range(len(indC)), top=regC + (variationC / 2), bottom=regC - (variationC / 2), width=0.5, fill_color="#D5E1DD")

output_file("C-Bacterias-manhattan-full-ordered.html", title="C aa and variation across aa for Bacterias")
show(p)  # open a browser
