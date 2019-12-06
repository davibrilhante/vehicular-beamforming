import numpy as np

def fftmtx(N, option=1):
    # Implemented by Aldebaro Klautau
    # Also Available at:
    # https://github.com/aldebaro/dsp-telecom-book-code/blob/master/MatlabOctaveFunctions/ak_fftmtx.m
    # function [Ah, A] = ak_fftmtx(N, option)
    # FFT 'DIRECT Ah' and 'INVERSE A' matrices with 3 options for
    # the normalization factors:
    # 1 (default) ->orthonormal matrices
    # 2->conventional discrete-time (DT) Fourier series
    # 3->used in Matlab/Octave, corresponds to samples of DTFT
    # Example that gives the same result as Matlab/Octave:
    # Ah=ak_fftmtx(4,3); x=[1:4]'; Ah*x, fft(x)

    W = np.exp(-1j*2*np.pi/N)#Twiddle factor W_N
    #Create the matrix with twiddle factors
    Ah = [[W**(i*j) for j in range(N)] for i in range(N)]
    A = []
    if option == 1:#Orthonormal or unitary
        Ah /= np.sqrt(N)
        A = np.conj(Ah)
    elif option == 2:#Read X(K) in Volts in case x(n) is in Volts
        A = np.conj(Ah)
        Ah /= N
    elif option == 3: #As in matlab/octave Ah = Ah
        A = np.conj(Ah)
    else:
        print('Invalid option value: %d' % option)
        return 0

    return [Ah, A]


frequency = 60e9
c = 3e8

Lambda = c/frequency

nElementsX = 4
nElementsY = 4
normDist = 0.5
elementsDist = normDist * Lambda

numRandomVectorsX = int(nElementsX/2)
numRandomVectorsY = int(nElementsY/2)

numDFTPairVectorsX = int(nElementsX/2)
numDFTPairVectorsY = int(nElementsY/2)

numSteeredVectorsX = nElementsX
numSteeredVectorsY = nElementsY

numVectorsX = nElementsX + numRandomVectorsX + numDFTPairVectorsX + numSteeredVectorsX
numVectorsY = nElementsY + numRandomVectorsY + numDFTPairVectorsY + numSteeredVectorsY

Wx = [[0 for j in range(numVectorsX)] for i in range(nElementsX)]
Wy = [[0 for j in range(numVectorsY)] for i in range(nElementsY)]

dftMatrixX = fftmtx(nElementsX, 3)
dftMatrixY = fftmtx(nElementsY, 3)

####################################### Random Vectors #########################################

#Creates a random normal sample matrix with nElementsX by numRandomVectorsX dimension
randomX = np.random.normal(size=(nElementsX,numRandomVectorsX))+\
                1j*np.random.normal(size=(nElementsX,numRandomVectorsX))
randomY = np.random.normal(size=(nElementsY,numRandomVectorsY))+\
                1j*np.random.normal(size=(nElementsY,numRandomVectorsY))


####################################### DFT Vectors #############################################
dftPairsX = [[0 for i in range(numDFTPairVectorsX)] for j in range(nElementsX)]
dftPairsY = [[0 for i in range(numDFTPairVectorsY)] for j in range(nElementsY)]

for i in range(numDFTPairVectorsX):
    randNum1 = np.random.randint(0,nElementsX-1)
    randNum2 = np.random.randint(0,nElementsX-1)
    for j in range(nElementsX):
        dftPairsX[j][i] = dftMatrixX[0][j][randNum1] + dftMatrixX[0][j][randNum2]


for i in range(numDFTPairVectorsY):
    randNum1 = np.random.randint(0,nElementsY)
    randNum2 = np.random.randint(0,nElementsY)
    for j in range(nElementsY):
        dftPairsY[j][i] = dftMatrixY[0][j][randNum1] + dftMatrixY[0][j][randNum2]

##################################### STEERED Vectors ###########################################
# Simply steer along the y direction with elevation of 90 degrees
minAzimuth = np.radians(70)
maxAzimuth = np.radians(110)
elevationPoint = np.pi/2

steeredX = [[0 for i in range(numSteeredVectorsX)] for j in range(nElementsX)]
azimuths = np.linspace(minAzimuth,maxAzimuth,numSteeredVectorsX)

for i in range(numSteeredVectorsX):
    for j in range(nElementsX):
        beta = -2*np.pi*normDist*np.sin(elevationPoint)*np.cos(azimuths[i])
        steeredX[j][i] = np.exp(1j*j*beta)

steeredY = [[0 for i in range(numSteeredVectorsY)] for j in range(nElementsY)]
azimuths = np.linspace(minAzimuth,maxAzimuth,numSteeredVectorsY)
for i in range(numSteeredVectorsY):
    for j in range(nElementsY):
        beta = -2*np.pi*normDist*np.sin(elevationPoint)*np.cos(azimuths[i])
        steeredY[j][i] = np.exp(1j*j*beta)

#Populates Wx and Wy with the random matrix, DFT Pairs and steered vector
for i in range(nElementsX):
    counter = 0
    for j in range(numRandomVectorsX):
        Wx[i][nElementsX+j] = randomX[i][j]
    counter+=j

    for j in range(numDFTPairVectorsX):
        Wx[i][nElementsX+counter+j] = dftPairsX[i][j]
    counter+=j

    for j in range(numSteeredVectorsX):
        Wx[i][nElementsX+counter+j] = steeredX[i][j]


for i in range(nElementsY):
    counter = 0
    for j in range(numRandomVectorsY):
        Wy[i][nElementsY+j] = randomY[i][j]
    counter+=j

    for j in range(numDFTPairVectorsY):
        Wy[i][nElementsY+counter+j] = dftPairsY[i][j]
    counter+=j

    for j in range(numSteeredVectorsY):
        Wy[i][nElementsY+counter+j] = steeredY[i][j]

print(Wx)
