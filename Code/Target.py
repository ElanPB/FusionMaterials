import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import radioactivedecay as rd
import sys
import os
import Channel

class Target:
    evalSet = ['ENDF', 'JEFF'] # List of evaluated data sets to read in
    particles = {'g' : [0, 0], # [A, Z] for all particle that we will see
                 'n' : [1, 0],
                 'p' : [1, 1],
                 'd' : [2, 1],
                 't' : [3, 1],
                 'a' : [4, 2]}
    evalNames = ('X(MeV)', 'Y(barns)') # Name for columns for evaluated data
    experNames = ('X(MeV)', '+-dX(MeV)', 'Y(barns)', '+-dY(barns)', 'Year', 'Author(s)', 'EXFOR-ID') # Name for columns for experimental data
    delim = ((0, 14), (14, 27), (27, 40), (40, 53), (56, 61), (62, 82), (85, 93)) # Tuple of ranges for each column of data

    def __init__(self, isotope):
        self.isotope = isotope
        self.ele = rd.Nuclide(isotope)
        self.A = self.ele.A
        self.Z = self.ele.Z
        self.files = self.__getFiles()
        self.channels = {}
        self.__readFiles()
        self.chiSquared = {}
        self.outliers = {}
        self.plots = []
        self.__analyzeData()

    def mkdir(self):
        path = '../ProcessedData/'+self.isotope
        os.makedirs(path, exist_ok = True)

    def __getFiles(self):
        files = []
        for i in os.listdir('../RawData/'):
            if i.startswith(self.isotope):
                files.append(i)
        return files
    
    def parse(self, text):
        dataSets = {}

        name = ''
        evaluated = None
        read = False
        data = []
        col = []

        if not read:
            for i in text:
                # Look for // to signify the end of each data set
                if i.startswith('//'):
                    if len(data[0]) > 2:
                        col = self.experNames
                    else:
                        col = self.evalNames
                    # When end is found, save data from previous set as pandas Data Frame and reset variables
                    if name in dataSets.keys():
                        dataSets[name].append(pd.DataFrame(data, columns=col))
                    else:
                        dataSets[name] = pd.DataFrame(data, columns=col)
                    name = ''
                    evaluated = None
                    read = False
                    data.clear()
                
                # Filter out commented lines
                elif i.startswith('#'):
                    # Search for name of data set
                    if i.startswith('#name:'):
                        line = self.__strToArray(i)
                        for x in line:
                            for y in self.evalSet:
                                if x.startswith(y):
                                    evaluated = y
                                    break
                            if evaluated != None:
                                break
                        name = self.isotope + '(' + i.split('(')[1].split(')')[0] + ')'
                        if evaluated != None:
                            name += '_' + evaluated
                        if name in dataSets.keys():
                            read = True
                
                # If not commented and not the end of a data set, it is data
                else:
                    # Convert data from sting to array
                    info = []
                    for x in self.delim:
                        sec = i[x[0]:x[1]]
                        try:
                            sec = float(sec)
                        except:
                            sec = sec.strip()
                        if sec != '':
                            info.append(sec)

                    # Save data row of data to data set
                    data.append(info)
        
        return dataSets

    def __readFiles(self):
        prefix = '../RawData/'
        for f in self.files:
            file = open(prefix + f, 'r')
            lines = file.read().split('\n')
            self.sortReactions(self.parse(lines))

    def __analyzeData(self):
        for k in self.channels.keys():
            self.chi(k)
            self.filterChi(k)
            self.plot(k)
        self.saveFigs()

    def __strToArray(self, text):
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
    
    def chi(self, channel):
        data = self.channels[channel]
        chiSquared = {}
        for n, ev in data.items(): # Loop through data sets for a given reaction type
            if len(n.split('_')) == 2: # If data set is an evaluated data set (not experimental)
                exfor_chi = {} # Create directory for Chi-Squared values of experimental data sets
                for t, ex in data.items():
                    if len(t.split('_')) == 1: # Find experimental data sets
                        for key, grp in ex.groupby(['EXFOR-ID']): # Calculate Chi-Squared for EXFOR data by ID
                            chi = 0
                            for i in grp.index: # Calculate value inside sum for each point
                                y = self.__lerp(ev, grp['X(MeV)'][i])
                                chi += ((grp['Y(barns)'][i] - y)**2)/(grp['+-dY(barns)'][i]**2)
                            exfor_chi[str(key)] = [chi, grp.shape[0]] # Save Chi-Squared value and number of points by ID
                chiSquared[n] = exfor_chi # Save Chi-Squared values
        
        self.chiSquared[channel] = chiSquared
    
    def filterChi(self, channel):
        vals = self.chiSquared[channel]
        dataSets = {}
        if not vals:
            return {}
        for key, values in vals.items():
            outliers = self.findOutliers(values)
            dataSets[key] = outliers
        self.outliers[channel] = dataSets
    
    def findOutliers(self, data):
        outliers = []
        median = self.weightedMedian(data.values())
        offMedian = []
        for i in data.values():
            offMedian.append([abs(i[0] - median), i[1]])
        MAD = self.weightedMedian(offMedian)
        for k, i in data.items():
            val = 0.6745 * (i[0] - median) / MAD
            if val > 3.5:
                outliers.append(k)
        return outliers

    def weightedMedian(self, pairs):
        sortedNums = sorted(pairs, key=lambda x:x[0])
        if len(sortedNums) < 1:
            return -1
        front = sortedNums[0][1]
        end = sortedNums[-1][1]
        while len(sortedNums) > 1:
            if front > end:
                front -= end
                sortedNums.pop()
                end = sortedNums[-1][1]
            elif end > front:
                end -= front
                sortedNums.pop(0)
                front = sortedNums[0][1]
            else:
                if len(sortedNums) == 2:
                    return ((sortedNums[0][0] + sortedNums[1][0]) / 2)
                else:
                    sortedNums.pop(0)
                    sortedNums.pop()
                    front = sortedNums[0][1]
                    end = sortedNums[-1][1]
        return sortedNums[0][0]

    def combineChi(self, channel):
        vals = self.chiSquared[channel]
        if not vals:
            return {}
        
        combined = {}
        for k, i in vals.items():
            tot = 0
            points = 0
            for key, j in i.items():
                if not (key in self.outliers[channel][k]):
                    tot += j[0]
                    points += j[1]
            chi = tot/(points-1)
            combined[k] = chi
        return combined
    
    # Linear interpolation for a given x value
    def __lerp(self, points, x):
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
    
    def sortReactions(self, reactions):
        # Sort through reaction data sets and categorize by name
        for n, df in reactions.items():
            type = n.split('(')[1].split(')')[0] # determine type of reaction
            # If the reactions has been recorded with another dataset, append to current list
            if type in self.channels.keys():
                self.channels[type][n] = [df]
            # If first reaction of its kind make new entry in map for given type or reaction
            else:
                self.channels[type] = {n: df}
    
    def plot(self, channel):
        # Declaring variables for plot
        fig = plt.figure()
        fig.suptitle(channel)
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
        for n, df in self.channels[channel].items():
            # Plotting data
            # Plotting evaluated data
            if len(n.split('_')) == 2:
                num += 1
                plt.plot(df['X(MeV)'], df['Y(barns)'], label = n.split('_')[1])
            
            # Plotting experimental data
            else:
                # Plotting by EXFOR ID
                for key, grp in df.groupby(['EXFOR-ID']):
                    num += 1
                    k = key
                    if type(key) == float:
                        k = str(int(key))
                    plt.errorbar(grp['X(MeV)'], grp['Y(barns)'], xerr = grp['+-dX(MeV)'],\
                                yerr = grp['+-dY(barns)'], ls = 'none', marker = 'o',\
                                mfc = 'white', markersize = 4 , label=k)
            
            #Finding extremes for purpose of determining log vs linear plotting
            x_min = min(x_min, df['X(MeV)'].min())
            y_min = min(y_min, df['Y(barns)'].min())
            x_max = max(x_max, df['X(MeV)'].max())
            y_max = max(y_max, df['Y(barns)'].max())
        
        # Plot formatting
        box = ax.get_position()
        numCol = num//18 + 1
        ax.set_position([box.x0, box.y0, box.width * (5 - numCol) / 5, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = numCol)
        chi_squared = self.combineChi(channel)
        for n, i in chi_squared.items():
            print(i)
            if i != -1:
                fig.text(.5, .01, 'Chi-Squared (' + n + '): ' + str(i), ha='center')

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
        
        self.plots.append(fig)

    def saveFig(self):
        for i in self.plots:
            i.savefig('../ProcessedData/'+isotope+'/'+k+'.png')

    def getDaughter(self, channel):
        proj, eject = channel.split(',')
        dA = self.A + self.particles[proj][0]
        dZ = self.Z + self.particles[proj][1]

        if not(eject.lower() in ['el', 'inl']):
            out = eject.split('+')
            for i in out:
                print(i)
                num = 0
                while i[0].isnumeric():
                    num = num * 10 + int(i[0])
                    i = i[1:]
                if num == 0:
                    num = 1
                dA -= num * self.particles[i.lower()][0]
                dZ -= num * self.particles[i.lower()][1]

        ID = dZ * 1e7 + dA * 1e4
        name = self.translate(ID)
        return(name)
    
    def translate(self, ID):
        name = rd.Nuclide(ID).nuclide.split('-')
        name = name[1] + name[0]
        return name
    
    def getIsomers(self, iso):
        isomers = set([])
        while True:
            try:
                rd.Nuclide(iso)
            except:
                break
            isomers.add(iso)
            iso += 1
        return isomers

    def radProds(self, name):
        iso = rd.Nuclide(name)
        rad = self.getIsomers(iso.id)
        for i in iso.progeny():
            if i != 'SF':
                rad = rad.union(self.radProds(i))
        return rad
    
if __name__ == '__main__':
    isotope = sys.argv[1]
    