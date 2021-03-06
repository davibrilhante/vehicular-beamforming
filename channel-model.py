import numpy as np
import commpy
from upaCodebook import upaCodebookCreator, upaCodebookFullyRandom
from matplotlib import pyplot as plt

np.random.seed(1)



class Interaction:
    def __init__(self, Type, x, y, z):
        self.type = Type
        self.x = x
        self.y = y
        self.z = z

    def Print(self):
        print('[type: '+self.type+', x: %f, y: %f, z: %f]' % (self.x, self.y, self.z))

class Path:
    def __init__(self, n=1, pwr=0, toa=0, atheta=0, aphi=0, dtheta=0, dphi=0):
        self.nInteractions = n
        self.receivedPower =  pwr
        self.timeOfArrival = toa
        self.arrivalTheta = atheta
        self.arrivalPhi = aphi
        self.departureTheta = dtheta
        self.departurePhi = dphi
        self.interactions = []
        for i in range(n):
            self.interactions.append(Interaction('',0,0,0))

    def setInteraction(self, index, Type, x, y, z):
        self.interactions[index] = Interaction(Type, x, y, z)

    def Print(self):
        print('interactions: %d' % self.nInteractions)
        print('received power: %f' % self.receivedPower)
        print('time of arrival: %f' % self.timeOfArrival)
        print('arrival theta: %f' % self.arrivalTheta)
        print('arrival phi: %f'% self.arrivalPhi)
        print('departure theta: %f' % self.departureTheta)
        print('departure phi: %f' % self.departurePhi)
        

class Point:
    def __init__(self, length=10, pwr=0, toa=0, delay=0):
        self.nPaths = length
        self.avgReceivedPower = pwr
        self.avgTimeOfArrival = toa
        self.delaySpread = delay
        self.paths = []
        for i in range(length):
            self.paths.append(Path())

    def setPath(self, index, n, pwr, toa, atheta, aphi, dtheta, dphi):
        self.paths[index] = Path(n, pwr, toa, atheta, aphi, dtheta, dphi)

    def Print(self):
        print('[nPaths: %d, avgReceivedPower: %f, avgTimeOfArrival: %f, delaySpread: %f]' % (self.nPaths, self.avgReceivedPower, self.avgTimeOfArrival, self.delaySpread))

def raisedCosine(t, period, samplerate, rolloff):
    if abs(t) == period/(2*rolloff):
        return (np.pi/4*period)*np.sinic(1/(2*rolloff))
    else:
        return (1/period)*np.sinc(1/period)*(np.cos(np.pi*rolloff*t/period)/(1-(2*rolloff*t/period)**2))


Nr = [4,4] #4X4 antenna array
Nt = [4,4]

d = 0.5 #distance between antenna elements in function of lambda
frequency = 60e6 #Operating frequency
bandwidth = 1760e6
c = 3e8 # light speed in m/s
Lambda = c/frequency #wavelength

dx = Lambda*d #explicit distance between antenna elements in x
dy = Lambda*d

finalAK = []
finalRand = []
theorical = []

########################### Load Raw Data ######################################
filename = "model.paths.t001_01.r002.p2m"
f = open(filename)
lines = f.readlines()
counterPoints = -1
counterPaths = 0
counterInteractions = 0
pointData = True
points = []
types=[]

for l in lines:
    l = l.strip().split()
    if l[0] == '#': #skip the format of the file
        continue
    else:
        if counterPoints == -1: #search for the number of points
            nPoints = int(l[0])
            for i in range(nPoints):
                points.append(Point())
            counterPoints = 0
            continue

        if len(l) == 2: #search for the point index and the number of its paths
            points[int(l[0])-1] = Point(int(l[1]))
            pointData = True
            counterPaths = 0
            continue

        if pointData:
            points[counterPoints].avgReceivedPower = float(l[0])
            points[counterPoints].avgTimeOfArrival = float(l[1])
            points[counterPoints].delaySpread = float(l[2])
            #points[counterPoints].Print()
            pointData = False
            continue

        if counterPaths >= 0 and counterPaths < points[counterPoints].nPaths and not pointData:
            if len(l) == 8:
                #print(counterPoints, counterPaths)
                points[counterPoints].setPath(counterPaths, int(l[1]), float(l[2]), float(l[3]), 
                        np.radians(float(l[4])), np.radians(float(l[5])), np.radians(float(l[6])), np.radians(float(l[7])))
                #points[counterPoints].paths[counterPaths].Print()
                continue

            if len(l)==1:
                types = l[0].split('-')
                counterInteractions = 0
                continue

            else:
                if counterInteractions >=1 and counterInteractions < len(types)-1:
                    points[counterPoints].paths[counterPaths].setInteraction(counterInteractions-1, types[counterInteractions], float(l[0]),float(l[1]),float(l[2]))
                    #points[counterPoints].paths[counterPaths].interactions[counterInteractions-1].Print()
                    counterInteractions += 1
                    continue
                else:
                    if types[counterInteractions] == 'Rx':
                        pass
                    else:
                        counterInteractions += 1
                        continue

            counterPaths+=1
            if counterPaths == points[counterPoints].nPaths:
                counterPoints += 1


H = []
k = 2*np.pi/Lambda
rolloff = 0.1
symbolPeriod = 1/bandwidth
sampleF = 2*(frequency+bandwidth)
samples = 10000
rcosTime, rcosResponse = commpy.filters.rcosfilter(samples, rolloff, symbolPeriod, sampleF)

totalTxAntennaElements = Nt[0]*Nt[1] 
totalRxAntennaElements = Nr[0]*Nr[1] 

for p in points:
    counter = 0
    Hmn = []
    Hrc = []
    for m in range(totalTxAntennaElements):
        Hmn.append([])
        Hrc.append([])
        for n in range(totalRxAntennaElements):
            summ = 0
            summRC = 0
            for t in p.paths:
                txSteeringVector = [[],[]]
                rxSteeringVector = [[],[]]
                #t.Print()
                ########################### Steering Vectors ###################################
                for q in range(Nt[0]): 
                    omegaY = k*dy*np.sin(t.departureTheta)*np.sin(t.departurePhi)
                    txSteeringVector[0].append(np.exp(1j*q*omegaY))
                for q in range(Nt[1]): 
                    omegaX = k*dx*np.sin(t.departureTheta)*np.cos(t.departurePhi)
                    txSteeringVector[1].append(np.exp(1j*q*omegaX))

                for q in range(Nr[0]): 
                    omegaY = k*dy*np.sin(t.arrivalTheta)*np.sin(t.arrivalPhi)
                    rxSteeringVector[0].append(np.exp(1j*q*omegaY))
                for q in range(Nr[1]): 
                    omegaX = k*dx*np.sin(t.arrivalTheta)*np.cos(t.arrivalPhi)
                    rxSteeringVector[1].append(np.exp(1j*q*omegaX))

                #Finishes the calculum of the steering vectors multplying by the antenna factor (scalar)
                #and does the kronecker product between theta and phi components
                txSteeringVector = (1/np.sqrt(Nt[0]*Nt[1]))*np.kron(txSteeringVector[0],txSteeringVector[1])
                rxSteeringVector = (1/np.sqrt(Nr[0]*Nr[1]))*np.kron(rxSteeringVector[0],rxSteeringVector[1])
                randomTheta = np.radians(np.random.uniform(0,360))

                # As the gain is not available in a complex form, only real, the phase is randomly picked
                complexGain = t.receivedPower*(np.cos(randomTheta) + 1j*np.sin(randomTheta))

                #Multiplies the steering vectors among themselves, but the tx steering vectors suffers an
                #hermitian transformation or the conjugate transpose
                summ += np.sqrt(Nt[0]*Nt[1]*Nr[0]*Nr[1])*complexGain*np.dot(rxSteeringVector, txSteeringVector.conj().T)
                
                # result: array with the results with time near the nT - tl
                # choice: one of elements in result is choosen randomly
                # gfilter: impulse response function of raised cosine filter in the time chosen by choice and result
                '''
                result = np.where(np.isclose(rcosTime, counter*symbolPeriod - t.timeOfArrival))
                if len(result) > 1:
                    choice = np.random.choice(result[0])
                    gfilter = rcosResponse[choice]
                else:
                '''
                gfilter = raisedCosine(counter*symbolPeriod - t.timeOfArrival, symbolPeriod, sampleF, rolloff) #rcosResponse[choice]
                counter += 1
                    
                summRC += gfilter*summ
            Hmn[-1].append(summ)
            Hrc[-1].append(summRC)

    H.append(Hmn)
    #H.append(Hrc)
#print('H matrix dimensions: %d points with %d tx elements and %d rx elements' % (len(H), len(H[0]), len(H[0][0])))
#print(H[0])#Printing the Matrix of point 0

#Two different codebooks with tx and rx each
codebooks = [[[],[]],[[],[]]]
bestCodewords = [[[],[]],[[],[]]]
p = [[],[]]
channel = []
dummy = [1 for i in range(totalTxAntennaElements)]

for method in range(2):
#CODEBOOKS + CHANNEL
    if method == 0:
        print('======================AK Method==============================')
        codebooks[method][0] = upaCodebookCreator(frequency, Nt[0], Nt[1], d)
        codebooks[method][1] = upaCodebookCreator(frequency, Nr[0], Nr[1], d)
    elif method == 1:
        print('=====================Random Method===========================')
        codebooks[method][0] = upaCodebookFullyRandom(frequency, Nt[0], Nt[1], d)
        codebooks[method][1] = upaCodebookFullyRandom(frequency, Nr[0], Nr[1], d)

    for h in H:
        powers = []
        for w in codebooks[method][0][0]:
            for f in codebooks[method][1][0]:
                factor1 = np.dot(w,h)
                powers.append(np.dot(factor1, f))
        maxGain = np.max(np.power(np.absolute(powers),2))
        p[method].append(maxGain)
        maxPair = np.argmax(powers)
        maxTxCodebook = int(maxPair/len(codebooks[method][1][0]))
        maxRxCodebook = maxPair % len(codebooks[method][1][0])
        bestCodewords[method][0].append(maxTxCodebook)
        bestCodewords[method][1].append(maxRxCodebook)
        print(maxGain, maxTxCodebook, maxRxCodebook)
        if method == 1:
            continue
        else:
            M = np.dot(np.dot(dummy, h),dummy)
            channel.append(M)

#angles = np.angle(p)
gains = np.absolute(p)#np.real(p/np.cos(angles)) 
theorical = np.power(np.absolute(channel),2)#np.real(channel/np.cos(np.angle(channel)))
if len(gains[0]) == 10:
    finalAK.append(gains[0])
    finalRand.append(gains[1])

#theorical = np.multiply(np.sqrt(Nt[0]*Nt[1]*Nr[0]*Nr[1]),np.absolute(theorical))
#theorical = np.absolute(theorical)
#print(len(finalAK), len(np.mean(finalAK, axis=0)))
print(finalAK[0])
#plt.plot(np.mean(finalAK,axis=0), label='AK')
plt.plot(finalAK[0], label='AK')
#plt.plot(np.mean(finalRand,axis=0), label='Random')
plt.plot(finalRand[0], label='Random')
plt.plot(theorical, label='Channel')
plt.legend()
plt.show()
