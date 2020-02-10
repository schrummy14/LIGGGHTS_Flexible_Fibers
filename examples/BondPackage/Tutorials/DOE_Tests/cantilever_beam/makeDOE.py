
from numpy import linspace
from numpy import meshgrid

def make_DOE_txt(names, variables, outFile = "DOE_Params.txt"):
    with open (outFile, "w") as f:
        wrtieVar2File(f, range(1,len(variables[0])+1), names[0], 'int')
        for k in range(1,len(names)):
            wrtieVar2File(f, variables[k-1], names[k], 'double')

def wrtieVar2File(f, values, varName, valType='double'):
    f.write("variable %s universe &\n" % varName)

    if valType=='int':
        for k in range(len(values)-1):
            f.write("%i &\n" % values[k])
        f.write("%i\n\n" % values[-1])
    elif valType=='double':
        for k in range(len(values)-1):
            f.write("%25.20e &\n" % values[k])
        f.write("%25.20e\n\n" % values[-1])

def readDOEparams(fileName = "DOE_Params.txt"):
    varValues = []
    with open(fileName) as f:
        while True:
            line = f.readline()
            splitLine = line.split(' ')
            if splitLine[0] == 'variable':
                name = splitLine[1]
                val = []
                while True:
                    line = f.readline().strip()
                    splitLine = line.split(' ')
                    if len(splitLine[0]) > 0:
                        val.append(float(splitLine[0]))
                    else:
                        break
                varValues.append([name,val])
            elif len(line) == 0:
                break        
    return varValues




if __name__ == "__main__":
    N_bond_damp_vals = 4
    N_bond_Youngs_modulus = 4

    bond_damp_vals = linspace(0,10,N_bond_damp_vals)
    bond_Youngs_modulus_vals = linspace(1e9,1e10,N_bond_Youngs_modulus)

    variable_names = [
        "run_num",
        "bond_damp_val",
        "bond_youngs_modulus"
        ]

    bvs, bYms = meshgrid(bond_damp_vals,bond_Youngs_modulus_vals)

    make_DOE_txt(variable_names,[bvs.flatten(),bYms.flatten()], "DOE_Params.txt")

    print(readDOEparams())