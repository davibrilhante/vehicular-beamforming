import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

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
    Ah = []
    for i in range(N):
        Ah.append([])
        for j in range(N):
            Ah[-1].append(0)
            Ah[i][j]= W**(i*j) 
    A = []
    if option == 1:#Orthonormal or unitary
        Ah /= np.sqrt(N)
        A = np.conj(Ah)
    elif option == 2:#Read X(K) in Volts in case x(n) is in Volts
        A = np.conj(Ah)
        Ah /= N
    elif option == 3: #As in matlab/octave Ah = Ah
        A = np.conj(Ah)/N
    else:
        print('Invalid option value: %d' % option)
        return 0

    return [Ah, A]

def upaCodebookCreator(frequency = 60e9, nElementsX=4, nElementsY=4, normDist=0.5):
    '''
        Returns the codebook used in: [REF] paper aldebaro et al

        Inputs:
            freq = the frequency in hertx (default 60 GHz)
            elementsX = number of elements in x
            elementsY = number of elements in y
            dist = distance between elements normalized by lambda
    '''
    #frequency = 60e9
    c = 3e8

    Lambda = c/frequency

    #nElementsX = 4
    #nElementsY = 4
    #normDist = 0.5
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
    #print(dftMatrixX)
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
        for j in range(nElementsX):
            Wx[i][j] = dftMatrixX[0][i][j]

        for j in range(numRandomVectorsX):
            Wx[i][nElementsX+j] = randomX[i][j]
        counter+=numRandomVectorsX

        for j in range(numDFTPairVectorsX):
            Wx[i][nElementsX+counter+j] = dftPairsX[i][j]
        counter+=numDFTPairVectorsX

        for j in range(numSteeredVectorsX):
            Wx[i][nElementsX+counter+j] = steeredX[i][j]


    for i in range(nElementsY):
        counter = 0
        for j in range(nElementsY):
            Wy[i][j] = dftMatrixY[0][i][j]

        for j in range(numRandomVectorsY):
            Wy[i][nElementsY+j] = randomY[i][j]
        counter+=numRandomVectorsY

        for j in range(numDFTPairVectorsY):
            Wy[i][nElementsY+counter+j] = dftPairsY[i][j]
        counter+=numDFTPairVectorsY

        for j in range(numSteeredVectorsY):
            Wy[i][nElementsY+counter+j] = steeredY[i][j]


    ### Normalization
    Wx = np.divide(Wx, np.sqrt(sum(np.absolute(np.power(Wx,2)))))
    Wy = np.divide(Wy, np.sqrt(sum(np.absolute(np.power(Wy,2)))))

    #print(Wx)


    ##################################### Codebook Creation #########################################
    numCodebookIndices = nElementsX*nElementsY + numRandomVectorsX*numRandomVectorsY +\
            nElementsX*numDFTPairVectorsY + nElementsY*numDFTPairVectorsX +\
            numSteeredVectorsX*numSteeredVectorsY

    codebook = [[0,0] for i in range(numCodebookIndices)]
    numElement = 0

    for i in range(nElementsX):
        for j in range(nElementsY):
            codebook[numElement][0] = i
            codebook[numElement][1] = j
            numElement += 1


    #Random Combinations
    for i in range(numRandomVectorsX):
        for j in range(numRandomVectorsY):
            codebook[numElement][0] = nElementsX + i
            codebook[numElement][1] = nElementsY + j
            numElement += 1

    #Combine single and double DFT for x and y
    for i in range(nElementsX):
        for j in range(numDFTPairVectorsY):
            codebook[numElement][0] = i
            codebook[numElement][1] = nElementsY + numRandomVectorsY + j
            numElement += 1


    for i in range(nElementsY):
        for j in range(numDFTPairVectorsX):
            codebook[numElement][0] = i
            codebook[numElement][1] = nElementsX + numRandomVectorsX + j
            numElement += 1

    #Steering Vectors
    for i in range(numSteeredVectorsX):
        for j in range(numSteeredVectorsY):
            codebook[numElement][0] = nElementsX+numRandomVectorsX+numDFTPairVectorsX+i
            codebook[numElement][1] = nElementsY+numRandomVectorsY+numDFTPairVectorsY+j
            numElement += 1


    #print(codebook)
    #Use the Kronecker to represent the matrix for a pair of wx and wy as a single array
    #See  John Brady, Akbar Sayeed, Millimeter-Wave MIMO Transceivers - Chap 10
    #Section 10.5
    #http://dune.ece.wisc.edu/publications/ -> 2015

    numCodewords = len(codebook)
    W = [[0 for i in range(numCodewords)] for j in range(nElementsX*nElementsY)]


    for i in range(nElementsX*nElementsY):
        for j in range(numCodewords):
            tempX = [Wx[k][codebook[j][0]] for k in range(nElementsX)]
            tempY = [Wy[k][codebook[j][1]] for k in range(nElementsY)]
            W[i][j] = np.kron(tempY, tempX)

    return W

def upaCodebookFullyRandom(frequency=60e9, nElementsX=4, nElementsY=4, normDist=0.5, bookLength=64):
    W = []
    for i in range(bookLength):
        randomX = np.random.normal(size=(nElementsX,nElementsX))+\
                        1j*np.random.normal(size=(nElementsX,nElementsX))
        randomY = np.random.normal(size=(nElementsY,nElementsY))+\
                        1j*np.random.normal(size=(nElementsY,nElementsY))
        W.append(np.dot(randomX, randomY).reshape(-1))
    return [W]

def plotBeamPattern(W, nElementsX=4, nElementsY=4, normDist=0.5):
    '''
        Plots the beam pattern based on arrat factor as given in [REF]
        based on an input codebook W [array bidimensional N X M, where M = Nx*Ny]
    '''
    ##################################### Array Factor Plot #########################################
    angles = np.angle(W[0][0])
    modulus = np.real(W[0][0])/np.cos(angles)
    arrayFactor = []
    im = []
    #for i in range()
    for p in range(-90,270):
        phi = np.radians(p)
        temp = 0
        counter = 0
        for n in range(nElementsY):
            for m in range(nElementsX):
                temp += modulus[counter]*np.exp(1j*m*angles[counter])*np.exp(1j*n*angles[counter])*np.exp(1j*2*np.pi*normDist*(m*np.cos(phi) + n*np.sin(phi)))
                counter += 1
        arrayFactor.append(temp)


    plt.polar(np.radians(range(-90, 270)), arrayFactor)
    #plt.plot(arrayFactor)
    plt.show()

if __name__ =='__main__':
    from sys import argv
    print(len(argv))
    if len(argv) > 1:
        W = upaCodebookCreator(float(argv[1]), int(argv[2]), int(argv[3]), float(argv[4]))
        plotBeamPattern(W, int(argv[2]), int(argv[3]), float(argv[4]))
    elif len(argv) <= 1:
        W = upaCodebookCreator()
        plotBeamPattern(W)
        W = upaCodebookFullyRandom()
        #print(W[0][0])
        plotBeamPattern(W)

