import pandas as pd
import radioactivedecay as rd
import matplotlib.pyplot as plt
import numpy as np

class Channel:
    
    particles = {'g' : [0, 0],
                 'n' : [1, 0],
                 'p' : [1, 1],
                 'd' : [2, 1],
                 't' : [3, 1],
                 'he3' : [3, 2],
                 'a' : [4, 2]}

    def __init__(self, title, data):
        self._data = {title : data}
        self._tar, self._proj, self._eject = self.components(title)
        self._name = self._tar + '(' + self._proj + ',' + self._eject + ')'
        self._ele = self.getEle()
        self._A = self._ele.A
        self._Z = self._ele.Z
        self._daughter = self.findDaughter()

    def exists(self, title):
        return([self._proj, self._eject] == self.components(title))
    
    def getEle(self):
        tar = self._tar
        print(tar)
        while tar.startswith('0'):
            tar = tar[1:]
        return rd.Nuclide(tar)
            
    # Append new data to data dictionary
    def addData(self, name, data):
        print(self._data.keys())
        if name in self._data.keys():
            self._data[name].append(data)
        else:
            self._data[name] = data
        print(self._data.keys())
    
    # Call Chi-Squared related function and plot data
    def analyze(self):
        self._chiSquared = self.chi()
        self._outliers = self.filterChi()
        self._chiVals = self.combineChi()
        self._plot = self.plot()
    
    # Find the target, projectiles, and ejectiles for this reaction channel
    def components(self, title):
        for i in ['(', ',', ')']:
            title = title.replace(i, '_')
        return title.split('_')[:3]

    # Use components to find the daughter of given reaction channel
    def findDaughter(self):
        dA = self._A + self.particles[self._proj.lower()][0]
        dZ = self._Z + self.particles[self._proj.lower()][1]

        if not(self._eject.lower() in ['el', 'inl']):
            out = self._eject.split('+')
            for i in out:
                num = 0
                while i[0].isnumeric():
                    num = num * 10 + int(i[0])
                    i = i[1:]
                if num == 0:
                    num = 1
                dA -= num * self.particles[i.lower()][0]
                dZ -= num * self.particles[i.lower()][1]

        ID = int(dZ * 1e7 + dA * 1e4)
        name = self.translate(ID)
        return(name)
    
    # Change formate of isotope name from ID from radioactivedecay to ZZZEEN
    def translate(self, ID):
        name = rd.Nuclide(ID).nuclide.split('-')
        name = name[1] + name[0]
        return name
    
    # Find all isomers of a given isotope
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

    # Find all isotopes below the given one in a decay chain along with potential isomers
    def radProds(self, name):
        iso = rd.Nuclide(name)
        rad = self.getIsomers(iso.id)
        for i in iso.progeny():
            if i != 'SF':
                rad = rad.union(self.radProds(i))
        return rad
    
    # Find the Chi-Squared value for all EXFOR-ID with all evaluated data sets
    def chi(self):
        chiSquared = {}
        for n, ev in self._data.items(): # Loop through data sets for a given reaction type
            if len(n.split('_')) == 2: # If data set is an evaluated data set (not experimental)
                exfor_chi = {} # Create directory for Chi-Squared values of experimental data sets
                for t, ex in self._data.items():
                    if len(t.split('_')) == 1: # Find experimental data sets
                        for key, grp in ex.groupby(['EXFOR-ID']): # Calculate Chi-Squared for EXFOR data by ID
                            chi = 0
                            for i in grp.index: # Calculate value inside sum for each point
                                y = self.lerp(ev, grp['X(MeV)'][i])
                                chi += ((grp['Y(barns)'][i] - y)**2)/(grp['+-dY(barns)'][i]**2)
                            exfor_chi[str(key)] = [chi, grp.shape[0]] # Save Chi-Squared value and number of points by ID
                chiSquared[n] = exfor_chi # Save Chi-Squared values
        return chiSquared
    
    # Helper function for chi, linear interpolates for a given x value
    def lerp(self, points, x):
        closest = points.iloc[(points['X(MeV)']-x).abs().argsort()[:1]].index.tolist()[0]
        if points['X(MeV)'][closest] < x:
            low = closest
            high = closest + 1
        else:
            low = closest - 1
            high = closest
        print(len(points['Y(barns)']))
        print(high)
        print(self._name)
        print('************************************************')
        slope = (points['Y(barns)'][high] - points['Y(barns)'][low])/(points['X(MeV)'][high] - points['X(MeV)'][low])
        lerp = slope * (x - points['X(MeV)'][low]) + points['Y(barns)'][low]
        return lerp
  
    # Find all EXFOR-IDs for each evaluated data sets which are not outliers
    def filterChi(self):
        vals = self._chiSquared
        dataSets = {}
        if not vals:
            return {}
        for key, values in vals.items():
            outliers = self.findOutliers(values)
            dataSets[key] = outliers
        return dataSets
    
    # Helper function for filterChi, takes in list of Chi-Squared values and find which are outliers
    def findOutliers(self, values):
        outliers = []
        median = self.weightedMedian(values.values())
        if median == -1:
            return values
        offMedian = []
        for i in values.values():
            print('**************************************')
            print(i)
            print('**************************************')
            offMedian.append([abs(i[0] - median), i[1]])
        MAD = self.weightedMedian(offMedian)
        for k, i in values.items():
            val = 0.6745 * (i[0] - median) / MAD
            if val > 3.5:
                outliers.append(k)
        return outliers

    # Helper function for findOutliers
    def weightedMedian(self, pairs):
        if len(pairs) < 1:
            return -1
        sortedNums = sorted(pairs, key=lambda x:x[0])
        print(pairs)
        print(sortedNums)
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
                sortedNums.pop(0)
                sortedNums.pop()
                front = sortedNums[0][1]
                end = sortedNums[-1][1]
        return sortedNums[0][0]

    # Find combined Chi-Squared value for a all evaluated set without outliers
    def combineChi(self):
        vals = self._chiSquared
        if not vals:
            return {}
        
        combined = {}
        for k, i in vals.items():
            tot = 0
            points = 0
            for key, j in i.items():
                if not (key in self._outliers[k]):
                    tot += j[0]
                    points += j[1]
            chi = tot/(points-1)
            combined[k] = chi
        return combined
    
    # Create the cross section plot for a given channel
    def plot(self):
        # Declaring variables for plot
        fig = plt.figure()
        fig.suptitle(self._name)
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
        for n, df in self._data.items():
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
            x_max = min(x_min, (df['X(MeV)'] > 0).min())
            y_max = min(y_min, (df['Y(barns)'] > 0).min())
            x_max = max(x_max, df['X(MeV)'].max())
            y_max = max(y_max, df['Y(barns)'].max())
        
        # Plot formatting
        box = ax.get_position()
        numCol = num//18 + 1
        ax.set_position([box.x0, box.y0, box.width * (5 - numCol) / 5, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = numCol)
        for n, i in self._chiVals.items():
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
        
        return fig
