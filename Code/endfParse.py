import pandas as pd

# Function to parse data from txt files
def parse(isotope):
    # Read in file
    file = open('../RawData/' + isotope + '.txt', 'r')
    lines = file.read().split('\n')

    # Define variables used in loop
    reactions = []
    name = ''
    data = []
    col = []

    # Loop through all lines in the txt file
    for i, e in enumerate(lines):
        # Look for // to signify the end of each data set
        if e.startswith('//'):
            # When end is found, save data from previous set as pandas Data Frame and reset variables
            df = pd.DataFrame(data, columns=col)
            reactions.append([name, df])
            name = ''
            data.clear()
            col.clear()
        
        # Filter out commented lines
        elif e.startswith('#'):
            # Search for name of data set
            if e.startswith('#name:'):
                name = isotope + ' (' + e.split('(')[1].split(')')[0] + ')'
            
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

parse('158Gd')