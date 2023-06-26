import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import radioactivedecay as rd
import sys
import os
from Channel import Channel

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
        self.mkdir()
        self.ele = self.getEle()
        self.A = self.ele.A
        self.Z = self.ele.Z
        self.files = self.__getFiles()
        self.channels = []
        self.__readFiles()
        self.__analyzeData()

    def getEle(self):
        iso = self.isotope
        print(iso)
        while iso.startswith('0'):
            iso = iso[1:]
        return rd.Nuclide(iso)

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
        for i in self.channels:
            i.analyze()
        self.saveFig()

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
      
    def sortReactions(self, reactions):
        # Sort through reaction data sets and categorize by name
        for n, df in reactions.items():
            kind = n.split('(')[1].split(')')[0] # determine kind of reaction
            # If the reactions has been recorded with another dataset, append to current list
            exists = False
            if kind.endswith('TOT') or kind.endswith('X'):
                continue
            for c in self.channels:
                if c.exists(kind):
                    c.addData(n, df)
                    exists = True
                    break
            if not exists:
                self.channels.append(Channel(n, df))

    def saveFig(self):
        for i in self.channels:
            print(i._name)
            i._plot.savefig('../ProcessedData/'+isotope+'/'+i._name+'.png')

if __name__ == '__main__':
    isotope = sys.argv[1]
    Target(isotope)
    