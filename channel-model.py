import numpy as np


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




Nr = 4*4 #4X4 antenna array
Nt = 4*4

d = 0.5 #distance between antenna elements in function of lambda
frequency = 60e6 #Operating frequency
c = 3e8 # light speed in m/s
Lambda = c/frequency #wavelength

dx = Lambda*d #explicit distance between antenna elements in x
dy = Lambda*d


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
                points[counterPoints].setPath(counterPaths, int(l[1]), float(l[2]), float(l[3]), float(l[4]), float(l[5]), float(l[6]), float(l[7]))
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

for p in points:
    ########################### Steering Vectors ###################################
    txSteeringVector = [[1],[1]]
    rxSteeringVector = [[1],[1]]
    pathgain = t.
    for t in p.paths:
        #t.Print()
        omegaY = k*dy*np.sin(t.departureTheta)*np.sin(t.departurePhi)
        omegaX = k*dx*np.sin(t.departureTheta)*np.cos(t.departurePhi)
        txSteeringVector[0].append(np.exp(1j*omegaY))
        txSteeringVector[1].append(np.exp(1j*omegaX))

        omegaY = k*dy*np.sin(t.arrivalTheta)*np.sin(t.arrivalPhi)
        omegaX = k*dx*np.sin(t.arrivalTheta)*np.cos(t.arrivalPhi)
        rxSteeringVector[0].append(np.exp(1j*omegaY))
        rxSteeringVector[1].append(np.exp(1j*omegaX))
    #Finishes the calculum of the steering vectors multplying by the antenna factor (scalar)
    #and does the kronecker product between theta and phi components
    txSteeringVector = (1/np.sqrt(Nt))*np.kron(txSteeringVector[0],txSteeringVector[1])
    rxSteeringVector = (1/np.sqrt(Nr))*np.kron(rxSteeringVector[0],rxSteeringVector[1])

    #Multiplies the steering vectors among themselves, but the tx steering vectors suffers an
    #hermitian transformation or the conjugate transpose
    H.append(np.dot(rxSteeringVector, txSteeringVector.conj().T))
    print(H[-1])

