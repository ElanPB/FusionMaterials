import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mendeleev
import sys
import os

# Function to parse data from txt files and return as titled data frames
def parse(isotope):
    # Read in file
    file = open('../RawData/' + isotope + '_ENDF.txt', 'r')
    lines = file.read().split('\n')

    # Define variables used in loop
    evalSet = ['ENDF', 'JEFF']
    endfNames = ('X(MeV)', 'Y(barns)')
    exforNames = ('X(MeV)', '+-dX(MeV)', 'Y(barns)', '+-dY(barns)', 'Year', 'Author(s)', 'EXFOR-ID')
    delim = ((0, 14), (14, 27), (27, 40), (40, 53), (56, 61), (62, 82), (85, 93))
    reactions = {}
    name = ''
    evaluated = None
    data = []
    col = []

    # Loop through all lines in the txt file
    for i, e in enumerate(lines):
        # Look for // to signify the end of each data set
        if e.startswith('//'):
            if len(data[0]) > 2:
                col = exforNames
            else:
                col = endfNames
            # When end is found, save data from previous set as pandas Data Frame and reset variables
            if name in reactions.keys():
                reactions[name].append(pd.DataFrame(data, columns=col))
            else:
                reactions[name] = pd.DataFrame(data, columns=col)
            name = ''
            evaluated = None
            data.clear()
        
        # Filter out commented lines
        elif e.startswith('#'):
            # Search for name of data set
            if e.startswith('#name:'):
                line = strToArray(e)
                for x in line:
                    for y in evalSet:
                        if x.startswith(y):
                            evaluated = y
                            break
                    if evaluated != None:
                        break
                name = isotope + '(' + e.split('(')[1].split(')')[0] + ')'
                if evaluated != None:
                    name += '_' + evaluated
        
        # If not commented and not the end of a data set, it is data
        else:
            # Convert data from sting to array
            info = []
            for x in delim:
                sec = e[x[0]:x[1]]
                try:
                    sec = float(sec)
                except:
                    sec = sec.strip()
                if sec != '':
                    info.append(sec)

            # Save data row of data to data set
            data.append(info)
    
    return reactions    
# Function to convert a string/line of data to an array
def strToArray(text):
    # Split up line by commas and spaces
    text = text.replace(',', ' ')
    out = text.split(' ')

    # Filter out blank and useless elements
    i = 0
    while i < len(out):
        if out[i] in ['', '#', '##']:
            del out[i]
        else:
            i += 1
    return out

def plot(name, data, chi_squared):
    # Declaring variables for plot
    fig = plt.figure()
    fig.suptitle(name)
    fig.set_size_inches(8,5)
    ax = plt.subplot(111)

    # Variables used for deciding log vs linear
    cutoffRatio = 1000
    cutoffMax = 10
    x_min = np.inf
    y_min = np.inf
    x_max = 0
    y_max = 0
    num = 0

    # Looping through all data sets (ENDF, EXFOR, ...)
    for n, df in data.items():
        # Plotting data
        # Plotting evaluated data
        if len(n.split('_')) == 2:
            num += 1
            ax.plot(df['X(MeV)'], df['Y(barns)'], label = n.split('_')[1])
        
        # Plotting experimental data
        else:
            # Plotting by EXFOR ID
            for key, grp in df.groupby(['EXFOR-ID']):
                num += 1
                k = key
                if type(key) == float:
                    k = str(int(key))
                ax.errorbar(grp['X(MeV)'], grp['Y(barns)'], xerr = grp['+-dX(MeV)'],\
                             yerr = grp['+-dY(barns)'], ls = 'none', marker = 'o',\
                             mfc = 'white', markersize = 4 , label=k)
        
        #Finding extremes for purpose of determining log vs linear plotting
        x_min = min(x_min, df['X(MeV)'].min())
        y_min = min(y_min, df['Y(barns)'].min())
        x_max = max(x_max, df['X(MeV)'].max())
        y_max = max(y_max, df['Y(barns)'].max())
    
    # Plot formatting
    box = ax.get_position()
    numCol = num // 18 + 1
    ax.set_position([box.x0, box.y0, box.width * (5-numCol) / 5, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = numCol)
    for i in chi_squared.values():
        if i != -1:
            fig.text(.5, .01, 'Chi-Squared: ' + str(i), ha='center')

    # Scaling log vs linear
    try:
        if x_max/x_min > cutoffRatio:
            plt.xscale('log')
    except ZeroDivisionError:
        plt.xscale('log')

    if y_max > cutoffMax:
        try:
            if y_max/y_min > cutoffRatio:
                plt.yscale('log')
        except ZeroDivisionError:
            plt.yscale('log')

    return fig

# Determine Chi-Squared value for evaluation compared to experimental
def chiSquared(data):
    chi_squared = {}

    for n, ev in data.items(): # Loop through data sets for a given reaction type
        if len(n.split('_')) == 2: # If data set is an evaluated data set (not experimental)
            exfor_chi = {} # Create directory for Chi-Squared values of experimental data sets
            for t, ex in data.items():
                if len(t.split('_')) == 1: # Find experimental data sets
                    for key, grp in ex.groupby(['EXFOR-ID']): # Calculate Chi-Squared for EXFOR data by ID
                        chi = 0
                        for i in grp.index: # Calculate value inside sum for each point
                            y = lerp(ev, grp['X(MeV)'][i])
                            chi += ((grp['Y(barns)'][i] - y)**2)/(grp['+-dY(barns)'][i]**2)
                        exfor_chi[str(key)] = [chi, grp.shape[0]] # Save Chi-Squared value and number of points by ID
            chi_squared[n] = exfor_chi # Save Chi-Squared values
    
    return chi_squared

def filterChi(vals):
    

# Take in Chi-Squared values for multiple data sets and combine into 
def combineChi(vals):
    if not vals:
        return -1
    
    combined = {}
    for k, i in vals.items():
        tot = 0
        points = 0
        for j in i.values():
            tot += j[0]
            points += j[1]
        chi = tot/(points-1)
        combined[k] = chi
    return combined

def lerp(points, x):
    closest = points.iloc[(points['X(MeV)']-x).abs().argsort()[:1]].index.tolist()[0]
    if points['X(MeV)'][closest] < x:
        low = closest
        high = closest + 1
    else:
        low = closest - 1
        high = closest
    slope = (points['Y(barns)'][high] - points['Y(barns)'][low])/(points['X(MeV)'][high] - points['X(MeV)'][low])
    lerp = slope * (x - points['X(MeV)'][low]) + points['Y(barns)'][low]
    return lerp

# Sort reactions into map of types of reactions (i.e. n,p) and matching data sets by name
def sortReactions(reactions):
    types = {} # map of types of reactions

    # Sort through reaction data sets and categorize by name
    for i in reactions.keys():
        type = i.split('(')[1].split(')')[0] # determine type of reaction
        # If the reactions has been recorded with another dataset, append to current list
        if type in types.keys():
            types[type].append(i)
        # If first reaction of its kind make new entry in map for given type or reaction
        else:
            types[type] = [i]

    return types

if __name__ == '__main__':
    isotope = sys.argv[1]
    path = '../ProcessedData/'+isotope
    path = (os.path.abspath(path))
    os.makedirs(path, exist_ok = True)
    reacts = parse(isotope)
    types = sortReactions(reacts)
    for k, v in types.items():
        chi = combineChi(chiSquared({key: reacts[key] for key in v}))
        fig = plot(k, {key: reacts[key] for key in v}, chi)
        fig.savefig('../ProcessedData/'+isotope+'/'+k+'.png')
