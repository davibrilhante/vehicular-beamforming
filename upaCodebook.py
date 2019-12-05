import numpy as np

frequency = 60e9
c = 3e8

Lambda = c/frequency

nElementsX = 4
nElementsY = 4
elementsDist = 0.5 * Lambda

numRandomVectorsX = nElementsX/2
numRandomVectorsY = nElementsY/2

numDFTPairVectorsX = nElementsX/2
numDFTPairVectorsY = nElementsY/2

numSteeredVectorsX = nElementsX
numSteeredVectorsY = nElementsY

numVectorsX = Nax + numRandomVectorsX + numDFTPairVectorsX + numSteeredVectorsX
numVectorsY = Nay + numRandomVectorsY + numDFTPairVectorsY + numSteeredVectorsY

Wx = [[0 for i in range(numVectorsX)] for j in range(nElementsX)]
Wy = [[0 for i in range(numVectorsy)] for j in range(nElementsY)]

#dftMatrixX = 
#dftMatrixY = 

####################################### Random Vectors #########################################
randomX = np.random.normal(nElementsX,numRandomVectorsX)+1j*np.random.(nElementsX,numRandomVectorsX)
randomY = np.random.normal(nElementsY,numRandomVectorsY)+1j*np.random.(nElementsY,numRandomVectorsY)

Wx[:][nElementsX+1:nElementsX+numRandomVectorsX] = randomX
Wx[:][nElementsY+1:nElementsY+numRandomVectorsY] = randomY

####################################### DFT Vectors #############################################
dftPairsX = [[0 for i in range(numDFTPairVectorsX)] for j in range(nElementsX)]
dftPairsY = [[0 for i in range(numDFTPairVectorsY)] for j in range(nElementsY)]

for i in range(numDFTPairVectorsX):
    randNum1 = np.random.randint(0,nElementsX)
    randNum2 = np.random.randint(0,nElementsX)
    dftPairsX[:][i] = 1


for i in range(numDFTPairVectorsY):
    randNum1 = np.random.randint(0,nElementsY)
    randNum2 = np.random.randint(0,nElementsY)
    dftPairsY[:][i] = 1


