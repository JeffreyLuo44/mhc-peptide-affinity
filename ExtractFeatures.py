def typeOf(AA):
    switch_dict = {
        'D': "Acidic",
        'E': "Acidic",
        'H': "Basic",
        'K': "Basic",
        'R': "Basic",
        'A': "Hydrophobic",
        'C': "Hydrophobic",
        'F': "Hydrophobic",
        'I': "Hydrophobic",
        'L': "Hydrophobic",
        'M': "Hydrophobic",
        'P': "Hydrophobic",
        'V': "Hydrophobic",
        'W': "Hydrophobic",
        'G': "Neutral",
        'N': "Neutral",
        'Q': "Neutral",
        'S': "Neutral",
        'T': "Neutral",
        'Y': "Neutral"
    }
    return switch_dict.get(AA, None)

def sizeOf(AA):
    switch_dict = {
        'A': "XS",
        'G': "XS",
        'S': "XS",
        'N': "S",
        'D': "S",
        'C': "S",
        'P': "S",
        'T': "S",
        'V': "M",
        'H': "M",
        'E': "M",
        'Q': "M",
        'I': "L",
        'L': "L",
        'M': "L",
        'K': "L",
        'R': "L",
        'F': "XL",
        'W': "XL",
        'Y': "XL",
    }
    return switch_dict.get(AA, None)


def volumeOf(AA):
    switch_dict = {
        'A': 88.6,
        'G': 60.1,
        'S': 89.0,
        'N': 114.1,
        'D': 111.1,
        'C': 108.5,
        'P': 112.7,
        'T': 116.1,
        'V': 140.0,
        'H': 153.2,
        'E': 138.4,
        'Q': 143.8,
        'I': 166.7,
        'L': 166.7,
        'M': 162.9,
        'K': 168.6,
        'R': 173.4,
        'F': 189.9,
        'W': 227.8,
        'Y': 193.6,
    }
    return switch_dict.get(AA, None)


def hydropatyOf(AA):
    switch_dict = {
        'A': 1.8,
        'G': -0.4,
        'S': -0.8,
        'N': -3.5,
        'D': -3.5,
        'C': 2.5,
        'P': -1.6,
        'T': -0.7,
        'V': 4.2,
        'H': -3.2,
        'E': -3.5,
        'Q': -3.5,
        'I': 4.5,
        'L': 3.8,
        'M': 1.9,
        'K': -3.9,
        'R': -4.5,
        'F': 2.8,
        'W': -0.9,
        'Y': -1.3,
    }
    return switch_dict.get(AA, None)


positionsAccordingPaper = [0,3,5,8]
positionsInvestigate = [0,1,2,3,4,5,6,7,8]
def features(feature):
    switch_dict = {
        "posAA":True,
        "typeAA" : True,
        "sizeAA" : True,
        "numberOfResidue" : True,
        "totalVol" : True,
        "totalHydropathy1469" : True,
        "averageHydropathy1469" : True,
        "averageHydropathy" : True,

    }
    return switch_dict.get(feature, None)

AA_list = ['D','E','H','K','R','A','C','F','I','L','M','P','V','W','G','N','Q','S','T','Y']

# positions = [0,3,5]


outputData = open("outputExtractFeatures.csv", "w")
outputData.write("length,")
for position in positionsInvestigate:
    if features("posAA"):
        outputData.write("pos"+str(position+1))
        outputData.write(",")
    if features("typeAA"):
        outputData.write("type"+str(position+1))
        outputData.write(",")
    if features("sizeAA"):
        outputData.write("size"+str(position+1))
        outputData.write(",")
    
# LABEL Number of residues
if features("numberOfResidue"):
    for residue in AA_list:
        outputData.write("numOf"+residue)
        outputData.write(",")
    
    
if features("totalVol"):
    outputData.write("totalVol")
    outputData.write(",")

if features("totalHydropathy1469"):
    outputData.write("totalHydropathy1469")
    outputData.write(",")

if features("averageHydropathy1469"):
    outputData.write("averageHydropathy1469")
    outputData.write(",")

if features("averageHydropathy"):
    outputData.write("averageHydropathy")
    outputData.write(",")


outputData.write("bind")
outputData.write("\n")

inputData = open("data.txt", "r")

lines = inputData.read().splitlines()
for line in lines:

    data = line.split(" ")
    seq = data[1]
    # Length of residues
    lengthStr = str(len(data[1]))
    outputData.write(lengthStr)
    outputData.write(",")
    
    
    totalHydropathy=0
    for position in positionsInvestigate:
        try:
                # Residue in position
            if features("posAA"):
                outputData.write(seq[position])
                outputData.write(",")
                # Type of residue
            if features("typeAA"):
                outputData.write(typeOf(seq[position]))
                outputData.write(",")
                # Size of residue
            if features("sizeAA"):
                outputData.write(sizeOf(seq[position]))
            if position in positionsAccordingPaper:
                totalHydropathy+=hydropatyOf(seq[position])
        except:
            if features("posAA"):
                outputData.write("None")
                outputData.write(",")
            if features("typeAA"):
                outputData.write("None")
                outputData.write(",")
            if features("sizeAA"):
                outputData.write("None")
        outputData.write(",")
    
    # Num of Residues
    if features("numberOfResidue"):
        for residue in AA_list:
            outputData.write(str(seq.count(residue)))
            outputData.write(",")

        # Volume total
    if features("totalVol"):
        totalVol=0
        for AA in seq:
            totalVol+=volumeOf(AA)
        outputData.write(str(totalVol))
        outputData.write(",")

        # Total hydropathy 1 4 6 9
    if features("totalHydropathy1469"):
        outputData.write(str(totalHydropathy))
        outputData.write(",")
        
        # Total hydropathy 1 4 6 9/ length
    if features("averageHydropathy1469"):
        outputData.write(str(totalHydropathy/len(seq)))
        outputData.write(",")
        
        # Average hydropathy
    if features("averageHydropathy"):
        totalHydropathy=0
        for AA in seq:
            totalHydropathy+=hydropatyOf(AA)
        outputData.write(str(totalHydropathy/len(seq)))
        outputData.write(",")

    #Class bind
    if data[0]=="1":
        outputData.write("Yes")
    else :
        outputData.write("No")
    outputData.write("\n")

outputData.close()
