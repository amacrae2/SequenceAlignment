'''
Created on Oct 7, 2013

@author: alecmacrae
'''

# project 1: sequence alignment

# establishes the arguments that should be used to run the script
import sys
input_file = sys.argv[1]
output_file = sys.argv[2]


def getScoreMatrix(lettersA, lettersB, f):
    """creates the score matrix and a dictionary mapping from letter pairs to indices in the score matrix using the input file"""
    letterMap = {}
    scoreMatrix = []
    for i in xrange(len(lettersA)):
        scoreMatrix.append([])
        for j in xrange(len(lettersB)):
            line = f.readline().strip().split()
            letterMap[(line[2],line[3])] = (line[0],line[1])
            scoreMatrix[i].append(line[4])
    return (letterMap,scoreMatrix)

def getScore(a, b, scoreMatrixTuple):
    """Calculates and returns a score for a pair of characters"""
    letterMap = scoreMatrixTuple[0]
    scoreMatrix = scoreMatrixTuple[1]
    coords = letterMap.get((a,b))
    score = scoreMatrix[int(coords[0])-1][int(coords[1])-1]
    score = float(score)
    return score

def findScore(matrix, i, j, addOn):
    """calculates the score of s of a give cell in one of the matrices"""
    return round(matrix[i][j] + addOn, 5)

def getScoreMap(MScore, IxScore, IyScore, matrixName):
    """finds and returns a dictionary mapping from scores to the matrix(s) that led to that score"""
    scoreMap = {MScore:"M"}
    if matrixName != "Iy":
        if scoreMap.get(IxScore) != None:
            value = scoreMap.get(IxScore)
            scoreMap[IxScore] = [value,"Ix"]
        else:
            scoreMap[IxScore] = ["Ix"]
    if matrixName != "Ix":
        if scoreMap.get(IyScore) != None:
            value = scoreMap.get(IyScore)
            scoreMap[IyScore] = [value,"Iy"]
        else:
            scoreMap[IyScore] = ["Iy"]
    return scoreMap

def fillMatrixM(matrices, row, col, A, B, scoreMatrixTuple, local):
    """fills in the matrices M and Mp (the pointer matrix for M"""
    M = matrices[0]
    Ix = matrices[1]
    Iy = matrices[2]
    Mp = matrices[3]
    matchScore = getScore(A[col-1], B[row-1], scoreMatrixTuple)
    MScore = findScore(M, row-1, col-1, matchScore)
    IxScore = findScore(Ix, row-1, col-1, matchScore)
    IyScore = findScore(Iy, row-1, col-1, matchScore)
    scoreMap = getScoreMap(MScore, IxScore, IyScore, "M")
    score = max(MScore,IxScore,IyScore)
    fromMatrix = scoreMap.get(score)
    if local and score<0:
        score = 0
    M[row][col] = score
    Mp[row][col] = fromMatrix
    matrices[0] = M
    matrices[3] = Mp
    return matrices

def fillMatrixI(matrices, row, col, d, e, scoreMatrixTuple, matrixName, lastRow, lastCol, local):
    """fills in the matrices Ix, Iy, Ixp, and Iyp"""
    M = matrices[0]
    Ix = matrices[1]
    Iy = matrices[2]
    Ixp = matrices[4]
    Iyp = matrices[5]
    if matrixName == "Ix":
        if col == lastCol:
            d = 0
            e = 0
        MScore = findScore(M, row-1, col, -d)
        IScore = findScore(Ix, row-1, col, -e)
        scoreMap = getScoreMap(MScore, IScore, "NA", "Ix")
    elif matrixName == "Iy":
        if row == lastRow:
            d = 0
            e = 0
        MScore = findScore(M, row, col-1, -d)
        IScore = findScore(Iy, row, col-1, -e)
        scoreMap = getScoreMap(MScore,"NA", IScore, "Iy")
    score = max(MScore, IScore)
    fromMatrix = scoreMap.get(score)
    if local and score<0:
        score = 0
    if matrixName == "Ix":
        Ix[row][col] = score
        Ixp[row][col] = fromMatrix
        matrices[1] = Ix
        matrices[4] = Ixp
    elif matrixName == "Iy":
        Iy[row][col] = score
        Iyp[row][col] = fromMatrix
        matrices[2] = Iy
        matrices[5] = Iyp
    return matrices

def findStartMatrix(matrices, A, B):
    """finds the matrix(s) to start traceback from"""
    lastRow = len(B)
    lastCol = len(A)
    MScore = matrices[0][lastRow][lastCol]
    IxScore = matrices[1][lastRow][lastCol]
    IyScore = matrices[2][lastRow][lastCol]
    scoreMap = getScoreMap(MScore, IxScore, IyScore, "M")
    score = max(MScore,IxScore,IyScore)
    startMatrix = scoreMap.get(score)
    return startMatrix

def getStartingFromMatrices(currMatrix, matrices, row, col):
    """Finds the lower right entry for the matrix(s) that start traceback"""
    if currMatrix == "M":
        fromMatrices = matrices[3][row][col]
    elif currMatrix == "Ix":
        fromMatrices = matrices[4][row][col]
    else:
        fromMatrices = matrices[5][row][col]
    return fromMatrices

def flatten(matrixName):
    """forms a list of de-tupled strings from possible tuple entries"""
    if len(matrixName) == 2:
        if len(matrixName[0]) == 2:
            result = [matrixName[0][0],matrixName[0][1],matrixName[1]]
            return result
    return matrixName

def getPair(row,col,A,B,currMatrix):
    """finds the nucleotide or AA pair based on which matrix traceback is in"""
    if currMatrix == "M":
        pair = [A[col-1], B[row-1]]
    elif currMatrix == "Ix" or currMatrix == ["Ix"]:
        pair = ["_", B[row-1]]
    else:
        pair = [A[col-1], "_"]
    return pair

def conSeq(currSeq, newPair):
    """concatenates the current sequence with the new pair"""
    if currSeq == []:
        return newPair
    for i in xrange(len(currSeq)):
        for j in xrange(2):
            currSeq[i][j] = currSeq[i][j]+newPair[j]
    return currSeq

def getMatrixPointer(row, col, pointerIndex, matrices):
    """finds the pointer(s) to the matrix(s) that led to the score in the current location of a matrix"""
    return matrices[pointerIndex][row-1][col-1]

def fromMatrixM(matrixName,row,col):
    if matrixName == "M":
        return matrices[0][row-1][col-1]
    elif matrixName == "Ix":
        return matrices[0][row-1][col]
    else:
        return matrices[0][row][col-1]

def getIndex(matrixName):
    if matrixName == "M":
        return 0
    elif matrixName == "Ix":
        return 1
    else:
        return 2
    
# the row index runs from 0 to the length of the matrix both ends inclusive
# same for the column index
def solve(row, col, matrixNames, currMatrix, A, B, matrices, local):
    """recursive finction that traces back through the matrices to find the output"""
    soFar = []
    fromMatrices = flatten(matrixNames)
    
    # base case: in local alignment traceback stops when the score reaches 0
    if local and "M" in fromMatrices and fromMatrixM(currMatrix,row,col) == 0 and matrices[getIndex(currMatrix)][row][col]>0:
        soFar.append(getPair(row,col,A,B,"M"))
        return soFar
    
    for fromMatrix in fromMatrices:
        pair = getPair(row,col,A,B,currMatrix)
            
        # base case: when row and col reach one, traceback is done
        if row == 1 and col == 1:
            soFar.append(pair)
            return soFar
        
        # recursive case: adds current pair to previous pairings
        else:
            if fromMatrix == "M":
                pairs = conSeq(solve(row-1,col-1,getMatrixPointer(row,col,3,matrices),fromMatrix,A,B, matrices,local), pair)
            elif fromMatrix == "Ix" or fromMatrix == ["Ix"]:
                pairs = conSeq(solve(row-1,col,getMatrixPointer(row,col,4,matrices),fromMatrix,A,B, matrices,local), pair)
            else:
                pairs = conSeq(solve(row,col-1,getMatrixPointer(row,col,5,matrices),fromMatrix,A,B, matrices,local), pair)
        soFar.extend(pairs)
    return soFar


def fillMatrices(A,B,dx,ex,dy,ey,scoreMatrixTuple, local):
    """fills in the matrices with the proper values based on affine gap penalty"""
    # defining constants
    negInf = -99999
    nCol = len(A) + 1
    
    # Initializing matrices (adding correct entries to 1st row and 1st column)
    M = [[negInf]*nCol]
    M[0][0] = 0
    Ix = [[negInf]*nCol]
    Iy = [[negInf]+[0]*(nCol-1)]
    Mp = [["x"]*nCol]
    Ixp = [["x"]*nCol]
    Iyp = [["x"]+["M"]+[["Iy"]]*(nCol-2)]
    for j in xrange(len(B)):
        M.append([negInf]+[0]*(nCol-1))
        Ix.append([0]*nCol)
        Iy.append([negInf]+[0]*(nCol-1))
        Mp.append(["x"]+[0]*(nCol-1))
        if j == 0:
            Ixp.append(["M"]+[0]*(nCol-1))
        else:
            Ixp.append([["Ix"]]+[0]*(nCol-1))
        Iyp.append(["x"]+[0]*(nCol-1))
    matrices = [M, Ix, Iy, Mp, Ixp, Iyp]
    # Initialization of matrices complete
    
    # populates the matrices and returns them
    for i in xrange(len(B)):
        row = i+1
        for j in xrange(len(A)):
            col = j+1
            matrices = fillMatrixM(matrices,row,col,A,B,scoreMatrixTuple,local)
            matrices = fillMatrixI(matrices,row,col,dx,ex,scoreMatrixTuple,"Ix",len(B),len(A),local)
            matrices = fillMatrixI(matrices,row,col,dy,ey,scoreMatrixTuple, "Iy",len(B),len(A),local)
    return matrices

def getStartPositions(matrixM, maxRow, maxCol):
    """Gets the starting positions and score for local alignment traceback"""
    topScoreIndecies = [-99999]
    for i in xrange(maxRow+1):
        for j in xrange(maxCol+1):
            score = matrixM[i][j]
            if score > topScoreIndecies[0]:
                topScoreIndecies = [score]
                topScoreIndecies.append([i,j])
            elif score == topScoreIndecies[0]:
                topScoreIndecies.append([i,j])
    return topScoreIndecies

def traceback(matrices,A,B, local):
    """starts the traceback through the matrices to find the paths taken to achieve the highest score"""
    if local:
        startMatrices = "M"
        startPositions = getStartPositions(matrices[0],len(B),len(A))
        startPositions = startPositions[1:]
    else:
        startMatrices = findStartMatrix(matrices,A,B)
        startPositions = [[len(B),len(A)]]
    
    # finds the list of matched sequences based on each of the starting matrices and concatenates
    # a list for all of them together and returns it
    sequences = []
    for currMatrix in startMatrices:
        for startPosition in startPositions:
            fromMatrices = getStartingFromMatrices(currMatrix,matrices,startPosition[0],startPosition[1])
            rowBuffer = 0
            colBuffer = 0
            if currMatrix == "Ix":
                colBuffer = 1
            elif currMatrix == "Iy":
                rowBuffer = 1
            sequences.extend(solve(startPosition[0]+rowBuffer,startPosition[1]+colBuffer,fromMatrices,currMatrix,A,B,matrices,local))
    return sequences

def maxScore(matrices, A, B):
    """finds the max score from all 3 matrices"""
    lastRow = len(B)
    lastCol = len(A)
    MScore = matrices[0][lastRow][lastCol]
    IxScore = matrices[1][lastRow][lastCol]
    IyScore = matrices[2][lastRow][lastCol]
    score = max(MScore,IxScore,IyScore)
    return score

def removeEndGaps(sequences):
    """Removes any gap-inclusive pairs at either end of the given sequences"""
    for i in xrange(len(sequences)): 
        while sequences[i][0][0] == "_" or sequences[i][1][0] == "_":
            for j in xrange(2):
                sequences[i][j] = sequences[i][j][1:]
        while sequences[i][0][-1] == "_" or sequences[i][1][-1] == "_":
            for j in xrange(2):
                sequences[i][j] = sequences[i][j][0:-1]
    return sequences
        
# initializes a set of values from the given input file 
f = open(input_file)
A = f.readline().strip()
B = f.readline().strip()
alignmentType = f.readline().strip()
if alignmentType == "0":
    local = False
else:
    local = True
gapPenalties = f.readline().strip()
gapPenalties = gapPenalties.strip().split()
dx = float(gapPenalties[0])
ex = float(gapPenalties[1])
dy = float(gapPenalties[2])
ey = float(gapPenalties[3])
NA = f.readline().strip()
lettersA = f.readline().strip()
NB = f.readline().strip()
lettersB = f.readline().strip()
scoreMatrixTuple = getScoreMatrix(lettersA, lettersB, f)
f.close()

        
matrices = fillMatrices(A, B, dx, ex, dy, ey, scoreMatrixTuple, local)
for i in xrange(6):
    print matrices[i]
sequences = traceback(matrices,A,B,local)
strippedSequences = removeEndGaps(sequences)
if local:
    outputScore = getStartPositions(matrices[0],len(B),len(A))
    outputScore = outputScore[0]
else:
    outputScore = maxScore(matrices,A,B)
outputScore = round(outputScore,1)
f = open(output_file, 'w')
f.write(str(outputScore)+"\n")
f.write("\n")
for seqPair in strippedSequences:
    for i in xrange(2):
        f.write(seqPair[i]+"\n")
    f.write("\n")
f.close()