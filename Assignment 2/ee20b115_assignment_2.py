import sys                                                                           # importing necessary libraries
import cmath
from numpy import *

VS = "V"
CS = "I"

class resistor:                                                                     # Classes for each circuit component
    def __init__(self, name, n1, n2, val):
        self.name = name
        self.value = float(val)
        self.node1 = n1
        self.node2 = n2

class inductor:
    def __init__(self, name, n1, n2, val):
        self.name = name
        self.value = float(val)
        self.node1 = n1
        self.node2 = n2

class capacitor:
    def __init__(self, name, n1, n2, val):
        self.name = name
        self.value = float(val)
        self.node1 = n1
        self.node2 = n2

class voltageSource:
    def __init__(self, name, n1, n2, val, phase=0):
        self.name = name
        self.value = float(val)
        self.node1 = n1
        self.node2 = n2
        self.phase = float(phase)

class currentSource:
    def __init__(self, name, n1, n2, val, phase=0):
        self.name = name
        self.value = float(val)
        self.node1 = n1
        self.node2 = n2
        self.phase = float(phase)

n = len(sys.argv)                                
if (n == 1):
    print("\nThe name of the netlist file is not provided.\n")                        # if only name of python file is provided, then this message
elif(n == 2):
    print("\nThe name of the netist file provided "+ sys.argv[1]+".\n")               # using 2nd argument as name of netlist file         

    cktfile = sys.argv[1]
    freq = 50                                                                         # initialising ac frequency
    component = { "R":[],"C":[],"L":[],VS:[],CS:[]}
    cktnode = []

    filelines = []
    with open (cktfile,"r") as f:
        for line in f.readlines():
            filelines.append(line.split('#')[0].split('\n')[0])
            if(line[:3] == '.ac'):
                freq = float(line.split()[2])                                         # getting circuit frequency and defining omega
                w = 2*pi*freq

        start = filelines.index(".circuit")
        end = filelines.index(".end")
        cktbody = filelines[start+1:end]
        for line in cktbody:

            token = line.split()                                                      # extracting each part of line as a token
            if token[1] not in cktnode:
                cktnode.append(token[1])
            if token[2] not in cktnode:
                cktnode.append(token[2])

            if token[0][0] == "R":                                                                                       # extracting resistor,capacitor and inductor info
                component["R"].append(resistor(token[0],token[1],token[2],token[3]))
            elif token[0][0] == "C":
                component["C"].append(capacitor(token[0],token[1],token[2],token[3]))
            elif token[0][0] == "L":
                component["L"].append(inductor(token[0],token[1],token[2],token[3]))
 
            elif token[0][0] == VS:                                                                                     # extracting voltage source and current source info
                if len(token) == 5: 
                    component[VS].append(voltageSource(token[0],token[1],token[2],float(token[4])))
                elif len(token) == 6:
                    if freq == 50:
                        sys.exit("frequency not given")
                    component[VS].append(voltageSource(token[0], token[1],token[2],float(token[4])/2,token[5]))
            elif token[0][0] == CS:
                if len(token) == 5: 
                    component[CS].append(currentSource(token[0],token[1],token[2],float(token[4])))
                elif len(token) == 6:
                    if freq == 50:
                        sys.exit("frequency not given")
                    component[CS].append(currentSource(token[0],token[1],token[2],float(token[4])/2,token[5]))

        nodenum = {cktnode[i]:i for i in range(len(cktnode))}                                                           # dictionary with node names and their numbers                                    
        numnode = len(cktnode)
        numVS = len(component[VS])

        M = zeros((numnode+numVS, numnode+numVS),complex)                                                               # creating matrices M and b          
        b = zeros((numnode+numVS,),complex)
                    
        for r in component["R"]:                                                                                        # equations for resistors,capacitors and inductors respectively
            if r.node1 != 'GND':
                M[nodenum[r.node1]][nodenum[r.node1]] += 1/r.value
                M[nodenum[r.node1]][nodenum[r.node2]] -= 1/r.value
            if r.node2 != 'GND':
                M[nodenum[r.node2]][nodenum[r.node1]] -= 1/r.value
                M[nodenum[r.node2]][nodenum[r.node2]] += 1/r.value
        for c in component["C"]:
            if c.node1 != 'GND':
                M[nodenum[c.node1]][nodenum[c.node1]] += complex(0, w*c.value)
                M[nodenum[c.node1]][nodenum[c.node2]] -= complex(0, w*c.value)
            if c.node2 != 'GND':
                M[nodenum[c.node2]][nodenum[c.node1]] -= complex(0, w*c.value)
                M[nodenum[c.node2]][nodenum[c.node2]] += complex(0, w*c.value)
        for l in component["L"]:
            if l.node1 != 'GND':
                M[nodenum[l.node1]][nodenum[l.node1]] += complex(0, -1.0/(w*l.value))
                M[nodenum[l.node1]][nodenum[l.node2]] -= complex(0, -1.0/(w*l.value))
            if l.node2 != 'GND':
                M[nodenum[l.node2]][nodenum[l.node1]] -= complex(0, -1.0/(w*l.value))
                M[nodenum[l.node2]][nodenum[l.node2]] += complex(0, -1.0/(w*l.value))
                           
        for i in range(len(component[VS])):                                                                                  # voltage source equations 

            if component[VS][i].node1 != 'GND':                                                                           
                M[nodenum[component[VS][i].node1]][numnode+i] = 1.0
            if component[VS][i].node2 != 'GND':
                M[nodenum[component[VS][i].node2]][numnode+i] = -1.0
                       
            M[0][0] = 1.0                                                                                                   # ground and other equations              
            M[numnode+i][nodenum[component[VS][i].node1]] = -1.0
            M[numnode+i][nodenum[component[VS][i].node2]] = +1.0
            b[numnode+i] = cmath.rect(component[VS][i].value, component[VS][i].phase*pi/180)

        for i in component[CS]:                                                                                             # current source equations 
            if i.node1 != 'GND':
                b[nodenum[i.node1]] = -1*i.value
            if i.node2 != 'GND':
                b[nodenum[i.node2]] = i.value
                 
x = linalg.solve(M, b)

print("The voltage amplitude at each node and current amplitude throgh current sources are:\n")
print(x)
                
