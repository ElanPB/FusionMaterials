import pandas as pd

def parse(isotope):
    file = open('../RawData/' + isotope + '.txt', 'r')
    lines = file.readlines()
    reactions = []
    name = ''
    data = []
    col = []
    for i, e in enumerate(lines):
        if e.startswith('//'):
            print(data[0])
            df = pd.DataFrame(data, columns=col)
            reactions.append([name, df])
            name = ''
            data.clear()
            col.clear()
        elif e.startswith('#'):
            if e.startswith('#name:'):
                name = isotope + ' (' + e.split('(')[1].split(')')[0] + ')'
            elif (i < len(lines) - 1 and lines[i + 1][0] != '#'):
                col = strToArray(e)
                for j, x in enumerate(strToArray(lines[i - 1])):
                    col[j] = x + '(' + col[j] + ')'
        else:
            line = strToArray(e)
            info = []
            for i in range(len(line)):
                try:
                    info.append(float(line[i]))
                except ValueError:
                    info.append(line[i])
            data.append(info)
    print(len(reactions))
    for i in range(min(5,len(reactions))):
        print(i)
        reactions[i][1].style

def strToArray(text):
    text = text.replace(',', ' ')
    out = text.split(' ')
    i = 0
    while i < len(out):
        if out[i] in ['', '#', '##']:
            del out[i]
        else:
            i += 1
    return out

parse('158Gd')