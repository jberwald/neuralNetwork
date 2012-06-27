from pylab import *

fit = []

def run(e0, dir, **args):
    
    fargs = {'step': '1e5',
             'direction': 0,
             'range': 20}

    for k in args.keys():
        fargs[k] = args[k]

    for i in range(fargs['range']):
        try:
            fh = open(dir+'fitness.'+str(i)+'.out', 'r')
        except:
            continue
        a = fh.read(12)
        print a
        fit.append(float(a[1:]))
        fh.close()

    if fargs['step'] is '1e5':
        step = 1e-5
    elif fargs['step'] is '1e4':
        step = 1e-4
    elif fargs['step'] is '1e3':
        step = 1e-3

    xt = [e0 + (step)*i for i in range(fargs['range'])]

    print "xt ", xt
    print "fit ", fit

    fig = figure()
    title('Change along eps'+str(fargs['direction'])+' direction')
    xlabel('epsilon')
    ylabel('fitness')
    plot(xt, fit)
    show()

    fig.clf()
    close(fig)

def plot_fitvec(eps, fitvec, d, fname='fitvec.eps'):

    n = len(fitvec)

    eps = [fitvec[i][1][d] for i in range(n)]
    vals = [fitvec[i][0] for i in range(n)]

    params = {'text.usetex': True}
    rcParams.update(params)

    fig = figure()
    title("Change in fitness along epsilon "+str(d)+" axis")
    xlabel("epsilon")
    ylabel("fitness")
    plot(eps,vals)
    savefig(fname)

    fig.clf()
    close(fig)
