#!/usr/bin/env python3

# arguments:: project vtu
# extracts flow parameters for a number of points
# from a vtu file

import vtk
import sys
#import math as *
import matplotlib.pyplot as plt
import numpy as np
#import math
#from scipy.special import erfc
#from scipy import interpolate
from scipy.interpolate import interp1d
import os

vtk_interval = 0 # set this number to read from whatever file you like
#SETTINGS ANALYTICAL SOLUTION
nodes = 5000
Tolerance_L1_NORM = 0.01
filename1 = 'analyt' # it is for the analytical solution
vtu_number1 = '0' # it is for analytical number

# serial
reader1 = vtk.vtkXMLUnstructuredGridReader()
reader1.SetFileName(filename1+'_'+str(vtu_number1)+'.vtu')

ugrid1 = reader1.GetOutputPort()


################# Initial and last coordinate of the probe ###################
# this is a line passes through the domain and extracts the data for comarision
x0 = 0.0
x1 = 1.0

y0 = 0.0333#0.0001 # 1.0/float(NUMBER)
y1 = y0 #<==Temporary, it can handle different values

z0 = 0.0
z1 = 0.0
#Resolution of the probe
resolution = nodes

hx = (x1 - x0) / resolution
hy = (y1 - y0) / resolution
hz = (z1 - z0) / resolution

detector = []
Experimental_X = []
for i in range(resolution+1):
    detector.append([hx * i + x0, hy * i + y0, hz * i + z0])
    Experimental_X.append(hx * i + x0)


################### Create the probe line on Analytical curve ################

points1 = vtk.vtkPoints()
points1.SetDataTypeToDouble()

for i in range(len(detector)-1):
    points1.InsertNextPoint(detector[i][0], detector[i][1], detector[i][2])

detectors1 = vtk.vtkPolyData()
detectors1.SetPoints(points1)    


probe1 = vtk.vtkProbeFilter()
probe1.SetInputConnection(ugrid1)


probe1.SetSourceConnection(ugrid1)
probe1.SetInputData(detectors1)
probe1.Update()

data1 = probe1.GetOutput()

#NAME OF THE VARIABLE YOU WANT TO EXTRACT DATA FROM
data_name1 = 'Tracer'  # dataname in your vtk files
FS1=[]
for j in range(points1.GetNumberOfPoints()):
    #print data.GetPointData().GetScalars(data_name)
    FS1.append(  data1.GetPointData().GetScalars(data_name1).GetTuple(j)) # values of your unknown
    
#Convert tuple to array
Analytical_Y = []
for item in FS1:
    Analytical_Y.extend(item)
Analytical_Y.extend(item)

#Create spline curve
f = interp1d(Experimental_X, Analytical_Y,kind ='cubic') # can be linear

########################### Experimental data #################################
path = os.getcwd()
#binpath = path + '/runmycase'
#os.system('rm -f ' + path+ '/*.vtu')
#os.system(binpath)

#RETRIEVE AUTOMATICALLY THE LAST VTU FILE
AutoNumber = 0
for files in os.listdir(path):
    if files.endswith(".vtu"):
        pos = files.rfind('_')
        pos2 = files.rfind('.')
        AutoFile = files[:pos]
        AutoNumber = max(AutoNumber, int(files[pos+1:pos2])-vtk_interval) #here 5 equals to vtk_interval in the main code


AutomaticFile = AutoFile # it gives Tracer
AutomaticVTU_Number = AutoNumber # it gives last vtu file number


################################AUTOMATIC STUFF###############################

try:
    filename   = sys.argv[1]
    vtu_number = int(sys.argv[2])
except:
    filename = AutomaticFile
    vtu_number = int(AutomaticVTU_Number)


# parallel
#reader = vtk.vtkXMLPUnstructuredGridReader()
#reader.SetFileName(filename+'_'+str(vtu_number)+'.pvtu')

# serial
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(filename+'_'+str(vtu_number)+'.vtu')

ugrid = reader.GetOutputPort()
#ugrid.Update()


################### Create the probe line on experimental curve################

points = vtk.vtkPoints()
points.SetDataTypeToDouble()

for i in range(len(detector)-1):
    points.InsertNextPoint(detector[i][0], detector[i][1], detector[i][2])

detectors = vtk.vtkPolyData()
detectors.SetPoints(points)    


probe = vtk.vtkProbeFilter()
probe.SetInputConnection(ugrid)


probe.SetSourceConnection(ugrid)
probe.SetInputData(detectors)
probe.Update()

data = probe.GetOutput()
###############################################################################

#NAME OF THE VARIABLE YOU WANT TO EXTRACT DATA FROM
data_name = 'Tracer'  # dataname in your vtk files
FS=[]
for j in range(points.GetNumberOfPoints()):
    #print data.GetPointData().GetScalars(data_name)
    FS.append(  data.GetPointData().GetScalars(data_name).GetTuple(j)) # values of your unknown

#Calculate analytical solution
Experimental_X = np.linspace(0, 1.0, num=nodes, endpoint=True, retstep=False, dtype=None)
Analytical_X = Experimental_X

#Analytical_Y = analytical_solution(Experimental_X)

#Convert tuple to array
Experimental_Y = []
for item in FS:
    Experimental_Y.extend(item)
    


############################# L-norm ##########################################
L1_sum = 0.0
L2_sum = 0.0
N_shock = 0
Infinite_Norm = 0.0
infy_norm = []

for i in range(len(Experimental_X)):
    if (i==0):#The first position is exact, so no need to interpolate
        offset = Analytical_Y[i] - Experimental_Y[i]
        L1_sum = L1_sum + abs(offset)
        continue
    
    position = Experimental_X[i]
#    x = getAnalytical_interpolated( Analytical_X, Analytical_Y, position)
    x = f(position)
    if (x==-1):
        print('The size of the Experimental and Analytical experiments are different')
        quit

    if (abs(x - Experimental_Y[i])> Infinite_Norm):
        Infinite_Norm = abs(x - Experimental_Y[i])

    offset = x - Experimental_Y[i]
    L2_sum = L2_sum + offset**2
    L1_sum = L1_sum + abs(offset)
    infy_norm.append(abs(offset))
     
       
L1_norm= L1_sum / len(Experimental_X) 
L2_norm = L2_sum**0.5  / len(Experimental_X)
Infinite_Norm = max(infy_norm)

Passed = True
if (L1_norm > Tolerance_L1_NORM): Passed = False


####################### Plot the results in 2d ##############################
showPlot = True
if (Passed): 
    print('3D BL works OK')
else:
    print('3D BL does NOT work')

print( "L1_norm:", L1_norm)
print( "L2_norm:", L2_norm)
print( "Infinite_Norm:", Infinite_Norm)


Analytical_Y = Analytical_Y[:-1]
if (showPlot):
    fig, ax = plt.subplots()
#    x = []
#    Analytical = []
#    for i in range(len(detector)):
#        x.append(float(detector[i][0]))
#        Analytical.append(float(FS[i][0]))
    line = plt.Line2D(Experimental_X, Analytical_Y, color='black', linewidth=2, label='Analytical')
    line2 = plt.Line2D(Experimental_X, FS, color='blue', linewidth=2, label='Experimental')
    #line.text.set_color('red')
    #line.text.set_fontsize(16)
    ax.add_line(line)
    ax.add_line(line2)
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.legend()
    plt.show()
