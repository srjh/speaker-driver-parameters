import numpy.fft as npfft
import scipy.io.wavfile
import scipy.signal
import matplotlib.pyplot as plt
import cmath
import math
import sys
import re

#Calculate from here:
#http://sound.whsites.net/tsp.htm


#Return the next-highest power of two for the given input.
def nextPowerTwo(i):
    return math.pow(2,1+math.floor(math.log(i,2)))

if len(sys.argv) == 1:
    print ("Pass the sample filename as a parameter, e.g.\n")
    print ("python tsparam.py test.wav\n")
    print ("Optionally pass the series resistance (10 Ohms assumed):\n")
    print ("python tsparam.py test.wav Rs=9.4\n")
    print ("See diagram.pdf for circuit details.")
    exit()
    
#Check input file. Add .wav if needed.      
try:
    arate, adata = scipy.io.wavfile.read(sys.argv[1])
except IOError:
    try:
        arate, adata = scipy.io.wavfile.read(sys.argv[1]+".wav")
    except IOError:
        print ("Invalid file.")
        raise  
    
#Check resistance input, assume 10 ohms if not given
for args in sys.argv:        
    regex=re.compile('Rs=(.*)') 
    regexOut = regex.match(args)
    if regexOut:
        Rseries=float(regexOut.group(1))
    else:
        Rseries=10.0
   
#Rseries is the resistance of the series resistor.
#
#This program assumes that the right channel is the voltage divided left channel.
#
#Connection to the recording device "line in" and amplifier output are as follows:
#
#AMP OUT -----+--------> LINE In Right
#             |
#           Ztest
#             |
#             +--------> LINE In Left
#             |   
#          Rseries
#             |
#GND----------+--------> LINE In GND
 
#Display File Properties
print ("Sample Rate:" + str(arate))
print ("Sample Count:" + str(len(adata)))
print ("Time: {:.2f} s".format(len(adata)/arate))
try:
    print ("Num Channels:" + str(len(adata[0]))) #Fails if error.
except:
    print ("Two channels are needed.")
    raise
    
#Limit input to 4 million samples
if len(adata)>4194304:
    bothSignal=adata[0:4194304]
    zeropad=4194304
else:
    bothSignal=adata
    zeropad=int(nextPowerTwo(len(adata)))
    
transb = list(zip(*bothSignal))
leftSignal = transb[0] #Reduced Signal
rightSignal = transb[1] #Full Signal

del transb,bothSignal,adata

#Fourier transforms. More efficient padded with zeros to 2^n.
leftFFT = npfft.fft(leftSignal,zeropad)
rightFFT = npfft.fft(rightSignal,zeropad)
freqFFT = npfft.fftfreq(zeropad, 1.0/arate)

#Transfer function
H = leftFFT/rightFFT

#Infer Impedance on the assumption that the transfer function
#is a voltage divider with the series resistor.
l = Rseries*(1-H)/(H)
del leftSignal,rightSignal,leftFFT,rightFFT

#Resample; take only ascending frequency
h2=scipy.signal.resample(H,2*arate)
resampledImpedance=scipy.signal.resample(l,2*arate)
resampledFreq=scipy.signal.resample(freqFFT,2*arate)
resampledImpedance=resampledImpedance[0:arate]
resampledFreq=resampledFreq[0:arate]

#Initial fit found by smoothing and peak finding:
b, a = scipy.signal.butter(4, 0.01)
smoothedFFT = scipy.signal.filtfilt(b, a, l)
del l
maxr = 0
maxru = 0
minz = 0
minzu = 0
maxz = 0
maxzu = 0
minr = 1e20
minru = 0
#Go between 10 Hz and the rising zero crossing of the reactance.
for u in range(0,len(smoothedFFT)//2):
    if (freqFFT[u] > 10):
        if smoothedFFT[u].imag > maxz:
            maxz=smoothedFFT[u].imag
            maxzu=u
        if smoothedFFT[u].imag < minz:
            minz=smoothedFFT[u].imag
            minzu=u
        if smoothedFFT[u].real > maxr:
            maxr = smoothedFFT[u].real
            maxru =u
        if smoothedFFT[u].real < minr:
            minr = smoothedFFT[u].real
            minru = u
        if smoothedFFT[u].imag >0 and smoothedFFT[u-1].imag<0:
            break

if len(smoothedFFT)/2 > 6*u:           
    slope = (smoothedFFT[5*u].imag-smoothedFFT[3*u].imag)/(freqFFT[5*u]-freqFFT[3*u])
else:
    slope = (smoothedFFT[len(smoothedFFT)/3].imag-smoothedFFT[len(smoothedFFT)/4].imag)/(freqFFT[len(smoothedFFT)/3]-freqFFT[len(smoothedFFT)/4])
le = slope/(2*math.pi)
            
#qms*(f/fs-fs/f)=1
qms1 = 1/(freqFFT[maxru]/freqFFT[maxzu]-freqFFT[maxzu]/freqFFT[maxru])
qms2 = 1/(freqFFT[minzu]/freqFFT[maxru]-freqFFT[maxru]/freqFFT[minzu])

#Write Frequency response to file
output = open(sys.argv[1]+'.csv','w')
output.write ("Freq(Hz) , Real Z(Ohm), Imag Z(j*Ohm), Z magnitude (Ohm)\n")
for x in range(1,arate):
    output.write(str(abs(resampledFreq[x])) + ',' + str(resampledImpedance[x].real)+ ',' + str(resampledImpedance[x].imag) + ',' + str(abs(resampledImpedance[x])) + ',' + str(h2[x].real) + ',' + str(h2[x].imag) +'\n')
output.close()

re = minr*0.95  #Guesstimate
qms = math.sqrt(qms1*qms2)
fs = freqFFT[maxru]
rec = maxr-re
qes = re*qms/rec
qts = qes*qms/(qms+qes)

print ("\nResults as follows:")
print ("Fs: %.2f" % freqFFT[maxru])
print ("Qms: %.2f" % qms)
print ("Qes: %.2f" % qes)
print ("Qts: %.2f" % qts)
print ("Zmax: %.2f" % maxr)
print ("Re: %.2f" % re)

#The following is looped as many times as the user requires.
#Parameters are as calculated initially, but can be refined.
refine = "y"
while refine is "y":
    simImpedanceReal = []
    simImpedanceImag=[]
    freqSim = [x*0.25 for x in range(1,40000)] 
    for f in freqSim:
      zc = -1j*rec*fs/(f*qms)
      zl = 1j*rec*f/(qms*fs)
      if zl+zc==0:
        zrlc = rec
      else:
        zlc = zc*zl/(zl+zc)
        zrlc = zlc*rec/(zlc+rec)
      zImp = zrlc + re + le*1j*f*2*math.pi + le*f*math.pi
      simImpedanceReal.append(zImp.real)
      simImpedanceImag.append(zImp.imag)
    print ("\nClose plot to refine parameters")
    plt.semilogx(resampledFreq,resampledImpedance.real, label='Measured Real Z')
    plt.semilogx(resampledFreq,resampledImpedance.imag,label='Measured Imag Z')
    plt.semilogx(freqSim,simImpedanceReal,'r--',label='Simulated Real Z')
    plt.semilogx(freqSim,simImpedanceImag,'c--', label='Simulated Imag Z')
    plt.xlabel("Freq (Hz)")
    plt.ylabel("Impedance (Ohm)")
    legend = plt.legend(loc='upper right')
    plt.grid(b=True,which='both',axis='both')
    plt.axis([10,min(arate/2,10000),minz*1.1,maxr*1.1])
    plt.show()
    refine=input("Manually refine parameters (y/n)?")
    if refine is "y":
        fs = float(input("Fs:"))
        qms = float(input("Qms:"))
        qes = float(input("Qes:"))
        re = float(input("Re:"))
        rec = re*qms/qes
        qts = qes*qms/(qms+qes)
        le = float(input("Le:"))

#Finally, calculate Vas from the change in Fs if known (with mass or volume):
vasCalc = input("Calculate Vas (y/n)?")
if vasCalc == "y":
    vasMode = input("Mass or volume mode (m/v)?")
    if vasMode=="m":
        madd = float(input("Added mass (g)?"))
        fsm = float(input("New Fs (Hz): "))
        coneA = float(input("Cone Area (cubic cm): "))
        cmass = madd/((fs/fsm)**2-1)
        print ("Cone mass:" + str(cmass))
        Vas = 0.0012*345**2*10*coneA**2/(4*math.pi**2*fs**2*cmass)
        print ("Vas:" + str(Vas))
    elif vasMode=="v":
        testVol = input("Test Volume (your choice of units):")
        fsm = input("New Fs (Hz): ")
        Vas = testVol*((fsm/fs)**2-1)
        print ("Vas (your choice of units):" + str(Vas))
    else:
        print ("Invalid choice, select m or v")
