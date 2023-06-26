import radioactivedecay as rd

def reaction(target, proj, eject):
    ele = rd.Nuclide(target)
    A = ele.A
    Z = ele.Z

    particles = {'g' : [0, 0],
                 'n' : [1, 0],
                 'p' : [1, 1],
                 'd' : [2, 1],
                 't' : [3, 1],
                 'a' : [4, 2]}
    
    dA = A + particles[proj][0]
    dZ = Z + particles[proj][1]

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
            dA -= num * particles[i.lower()][0]
            dZ -= num * particles[i.lower()][1]

    ID = dZ * 1e7 + dA * 1e4
    name = translate(ID)
    return(name)

def getIsomers(iso):
    isomers = set([])
    while True:
        try:
            rd.Nuclide(iso)
        except:
            break
        isomers.add(iso)
        iso += 1
    return isomers


def radProds(name):
    iso = rd.Nuclide(name)
    rad = getIsomers(iso.id)
    for i in iso.progeny():
        if i != 'SF':
            rad = rad.union(radProds(i))
    return rad

def translate(ID):
    name = rd.Nuclide(ID).nuclide.split('-')
    name = name[1] + name[0]
    return name

def order(name):
    return rd.Nuclide(name).id

prods = radProds('226Ra')
names = [translate(i) for i in prods]
names.sort(reverse=True, key=order)
print(names)