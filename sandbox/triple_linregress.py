from pylab import *
#from pylab import nan
from neural_network.netutils import utils

D = utils.Data_Handlers()

params = {'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 10,
          'font.family': 'serif',
          'text.usetex': True }

rcParams.update(params)

def plot_stuff(fname):

    d = fromfile(fname, sep='\n')
    d = d.reshape((len(d)/2,2))

    print d.shape

    a1 = where(d[:,0] < 2)[0]
    a2 = where(d[:,0] < 3.)[0]
    
    r1 = a1[-1]+1
    r2 = a2[-1]+1

    print "r1 ", r1
    print "r2 ", r2

    dvec = [(d[:r1],d[r1,0]), (d[r1:r2],d[r2,0]), (d[r2:],0)]

    lr = []
    for x in dvec:
        s = D.linregress(x[0][:,0],x[0][:,1])
        y = s[0]*x[0][:,0] + s[1]

        lr.append((x,y,s[0]))

    figure()
    clf()
    title("DFA for 5 cell network")
    xlabel("Log of window size n")
    ylabel("Log of DFA")

    color = ['g','r','m']
    i=0
    for a in lr:
        if a[0][1] != 0:
            plot(a[0][0][:,0],a[1], color=color[i],linewidth=1.5,label="linear reg. slope="+str(a[2]))
            i+=1
        else:
           plot(a[0][0][:,0],a[1], 'm', linewidth=1.5,label="linear reg. slope="+str(a[2]))
            
    plot(d[:,0],d[:,1], 'b-', linewidth=1, label="DFA curve")

    legend(loc=4)

    show()
        
