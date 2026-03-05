"""

Created by Robert Wolle, 9/19
Program that reads in data and fits it to Lorentzian peaks

"""


import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.



def make_DataArray(delimiter, data_columns, data_start, file_path):
    with open(file_path, "r") as file:
        data = file.read()

    data_values = []
    last_break = data_start
    for char in range(data_start,len(data)):
        if data[char] == delimiter or data[char] == '\n':
            if data[last_break:char] == '':     # Empty data blocks are saved as 0
                data_values.append(0)
            else:
                data_values.append(float(data[last_break:char]))
            last_break = char+1

    # We can use the number of data points divided by number of columns to get the number of rows of data
    data_rows = int(len(data_values)/data_columns)

    data_organized = []
    for row in range(data_rows):
        data_organized.append(data_values[row*data_columns:(row+1)*data_columns])
    return data_organized

def nm_to_eV(wavelength):
    return 1240/wavelength


# User inputs:
delimiter = '\t'
data_columns = 3
data_start = 521     #character position at which the data begins
file_path = "/workspaces/Projects/Physics/PeakFitting/WSe2_MoSe2_NbSe2_100uW_PL_Flake3Test8.vdat"


DataArray = np.array(make_DataArray(delimiter, data_columns, data_start, file_path))
x_axis = nm_to_eV(DataArray[:,0])
y_axis = DataArray[:,2]


x_region = [1.2,1.4]

plt.plot(x_axis,y_axis,zorder=2)
#plt.gca().invert_xaxis()
#plt.xlim(x_region)
#plt.ylim(500,1000)
plt.title('PL Intensity, on NbSe2 10s exposure')
plt.xlabel('Energy (eV)')
plt.ylabel('Counts')

plt.savefig('/workspaces/Projects/Physics/PeakFitting/output.png')





