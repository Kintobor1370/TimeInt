#**********************************************************************************************
#*                                    PLOT VISUALISER FILE                                   **
#*              The code below uses data from 'VisData.csv' file to make plots.              **
#**********************************************************************************************

import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import numpy as np
from array import array


colors = [
    'blue',
    'orange',
    'limegreen',
    'gold',
    'purple'
]


#_______________________________EXTRACTING DATA FROM VISDATA.CSV_______________________________
with open('Data\VisData.csv', 'r') as VisDataFile:
    FileReader = csv.reader(VisDataFile)

    ExperimentNumber = float(next(FileReader)[1])
    Title = next(FileReader)[1]
    next(FileReader)
    NumberOfPlots = int(next(FileReader)[1])
    next(FileReader)
    next(FileReader)
    Size = int(next(FileReader)[1])
    next(FileReader)
    next(FileReader)
    Labels = next(FileReader)
    
    x = array('d')
    y = [0] * NumberOfPlots
    for i in range(NumberOfPlots):
        y[i] = array('d')
    
    for new_line in FileReader:
        x.append(float(new_line[0]))
        for i in range(NumberOfPlots):
            y[i].append(float(new_line[i+1]))


#_____________________________________MAKING A SINGLE PLOT_____________________________________
plt.plot(x, y[0], color=colors[0])
plt.title(Title)
plt.xlabel(Labels[0])
plt.ylabel(Labels[1])

plt.xlim(x[0], x[Size-1])
if ExperimentNumber == 14.4:
    plt.ylim(0, 350)
elif ExperimentNumber == 14.6:
    eq = np.full(Size, 200)
    plt.ylim(0, 400)
    plt.plot(x, eq, linestyle = '--', color = 'grey')

plt.show()


#_____________________________________MAKING SEVERAL PLOTS_____________________________________
if NumberOfPlots > 1:
    for i in range(NumberOfPlots):
        plt.plot(x, y[i], color=colors[i])
    plt.title(Title)
    plt.xlabel(Labels[0])

    patches = [0] * NumberOfPlots
    for i in range(NumberOfPlots):
        new_patch = mp.Patch(color=colors[i], label=Labels[i+1])
        patches[i] = new_patch
    plt.legend(handles=patches)

    plt.xlim(x[0], x[Size-1])
    if ExperimentNumber == 14.4:
        plt.ylim(0, 350)
    elif ExperimentNumber == 14.6:
        plt.ylim(0, 400)
        plt.plot(x, eq, linestyle = '--', color = 'grey')

    plt.show()


#_______________________________MAKING AN ERROR PERCENTAGE PLOT________________________________
if NumberOfPlots == 2:
    er = [0] * Size
    for i in range(Size):
        er[i] = abs(y[0][i] - y[1][i])
        er[i] = er[i] / y[1][i]
        er[i] = er[i] * 100


    plt.plot(x, er, color='blue')
    plt.title(Title)
    plt.xlabel(Labels[0])
    plt.ylabel('Error (in %)')
    plt.xlim(x[0], x[Size-1])
    plt.ylim(0, 10)
    plt.show()


#________________________________MAKING PAIR CORRELATION PLOT__________________________________
with open('Data\CVisData.csv', 'r') as CVisDataFile:
    CReader = csv.reader(CVisDataFile)
    next(CReader)
    next(CReader)
    next(CReader)
    NumberOfPlotsC = int(next(CReader)[1])
    next(CReader)
    next(CReader)
    sizeC = int(next(CReader)[1])
    next(CReader)
    next(CReader)
    next(CReader)
    
    Cx = array('d')
    Cy = [0] * NumberOfPlotsC
    for i in range(NumberOfPlotsC):
        Cy[i] = array('d')
    
    for new_line in CReader:
        Cx.append(float(new_line[0]))
        for i in range(NumberOfPlotsC):
            Cy[i].append(float(new_line[i+1]))

length = NumberOfPlotsC
N_eq = y[0][Size-1]
dist = array('d')
pair_cor = array('d')
print(N_eq)
for i in range(length):
    dist.append(i)
    C_eq = Cy[i][sizeC-1]
    print(C_eq)
    pair_cor.append(C_eq / (N_eq * N_eq))

plt.plot(dist, pair_cor, color='blue')
plt.title("Pair correlation")
plt.ylabel("C / N^2")
plt.xlim(dist[0], dist[length-1])

print('\n============================PLOTS CREATED!=============================\n')
plt.show()
