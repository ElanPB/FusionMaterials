import pandas as pd

def parse(isotope):
    file = open('../RawData/' + isotope + '.txt', 'r')
    lines = file.read().split('\n')
    reactions = []
    name = ''
    data = []
    col = []
    for i, e in enumerate(lines):
        if e.startswith('//'):
            df = pd.DataFrame(data, columns=col)
            reactions.append([name, df])
            name = ''
            data.clear()
            col.clear()
        elif e.startswith('#'):
            if e.startswith('#name:'):
                name = isotope + ' (' + e.split('(')[1].split(')')[0] + ')'
            elif (i < len(lines) - 1 and len(lines[i + 1]) != 0 and lines[i + 1][0] != '#'):
                col = strToArray(e)
                for j, x in enumerate(strToArray(lines[i - 1])):
                    col[j] = x + '(' + col[j] + ')'
        else:
            line = strToArray(e)
            info = []
            for x in range(len(line)):
                try:
                    info.append(float(line[x]))
                except ValueError:
                    info.append(line[x])
            x = 0
            while x < len(info):
                if type(info[x]) == str:
                    y = x + 1
                    while (y < len(info) and type(info[y]) == str):
                        info[x] = info[x] + ' ' + info[y]
                        del info[y]
                x += 1

            data.append(info[:len(col)])
    for i in range(len(reactions)):
        print(reactions[i][1])

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