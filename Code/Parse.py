import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Function to parse data from txt files and return as titled data frames
def parse(isotope):
    # Read in file
    file = open('../RawData/' + isotope + '.txt', 'r')
    lines = file.read().split('\n')

    # Define variables used in loop
    evalSet = ['ENDF', 'JEFF']
    reactions = {}
    name = ''
    evaluated = None
    data = []
    col = []

    # Loop through all lines in the txt file
    for i, e in enumerate(lines):
        # Look for // to signify the end of each data set
        if e.startswith('//'):
            # When end is found, save data from previous set as pandas Data Frame and reset variables
            reactions[name] = pd.DataFrame(data, columns=col)
            name = ''
            evaluated = None
            data.clear()
            col.clear()
        
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

            # Search for data set coloumn names
            elif (i < len(lines) - 1 and len(lines[i + 1]) != 0 and lines[i + 1][0] != '#'):
                col = strToArray(e)
                for j, x in enumerate(strToArray(lines[i - 1])):
                    col[j] = x + '(' + col[j] + ')'
        
        # If not commented and not the end of a data set, it is data
        else:
            # Convert data from sting to array
            line = strToArray(e)
            info = []
            for x in range(len(line)):
                # Convert elements of array to numbers when possible
                try:
                    info.append(float(line[x]))
                except ValueError:
                    info.append(line[x])

            # Search for remaining string values in array and combine adjacent ones
            x = 0
            while x < len(info):
                if type(info[x]) == str:
                    y = x + 1
                    while (y < len(info) and type(info[y]) == str):
                        info[x] = info[x] + ' ' + info[y]
                        del info[y]
                x += 1
            
            # Save data row of data to data set
            data.append(info[:len(col)])
    
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

def plot(name, data):
    fig = plt.figure()
    fig.suptitle(name)
    for n, df in data.items():
        if len(n.split('_')) == 2:
            plt.plot(df['X(MeV)'], df['Y(barns)'], label = n.split('_')[1])
        else:
            plt.errorbar(df['X(MeV)'], df['Y(barns)'], xerr = df['+-dX(MeV)'], yerr = df['+-dY(barns)'], ls='none', label = df['EXFOR-ID'][0], marker = 'o', mfc='white')
    plt.legend(loc = 'best')
    return fig

def sortReactions(reactions):
    types = {}
    for i in reactions.keys():
        type = i.split('(')[1].split(')')[0]
        if type in types.keys():
            types[type].append(i)
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
        fig = plot(k, {key: reacts[key] for key in v})
        fig.savefig('../ProcessedData/'+isotope+'/'+k+'.png')
