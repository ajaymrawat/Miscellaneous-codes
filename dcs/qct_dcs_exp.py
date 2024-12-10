# This code calculates DCSs using histogram binning.
# It can calculate DCSs for any number of collision energies at once.
#import numpy as np
import math
pi = math.pi
au2ang = 0.529177249**2
#sab=file with name of binned files sab.dat (Binned means real values to closer integer values)
with open('sab.dat','r') as bsname, \
open('ics-total-gb-hlihp.dat','r') as csfile:  #values of cross-sections at different col en (Ecol vs ICS)
    for line in bsname:
        cs = csfile.readline()  #reading the cross-section
#       cs = float(cs.strip())  #removing the newline \n symbol .. 
        # ... and converting string to real no.
        cs = cs.strip()
        print(cs)
        cs = cs.split()
        print(cs)
        cs = float(cs[1])
        print(cs)
        nt = [0]*181       #grid initiation
        tp = 0.0
        tdcs = 0.0
        dt = (pi)/180.0
        ntraj = 0
        with open(line.strip(),'r') as bsfile:
            for data in bsfile:
                ntraj += 1
                data = int(data.strip())
# the 180-theta was used beacuse the 180-theta was used in code (wrong convention)
                nt[data] = nt[data] + 1
        s = 0
#       print('nt[0]: ',nt[0])
#       nt[1]=nt[0]+nt[1]
        print('total trajectories in file'+' '+line.strip()+'=  '+str(ntraj))
        with open('dcs'+line[2:7]+'.dat','w') as dcsfile:
            for i in range(1,180):
                p=nt[i]/ntraj
                s=s+nt[i]
#                print(i,nt[i],s, ' p: ',p)
                tp=tp+p
                if (i == 0):
                    rad=(pi*0.01)/180
                else:
                    rad=(pi*i)/180.0
                dcs=(cs*p)/(2*dt*pi*math.sin(rad))
                tdcs=tdcs+((dcs*(math.sin(rad)*2*pi))*dt) # dcs to ics
                dcsfile.write(str(i)+' '+str(dcs)+'\n')
        print('ics read from file'+' '+str(cs))
        print('total ics calculated'+' '+str(tdcs))
        dics=abs(cs-tdcs)
        print('difference of both ics= '+str(dics))
        print('total probability (should be = 1)'+' '+str(tp))
        print(' ')
