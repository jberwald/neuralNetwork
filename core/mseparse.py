#!/usr/bin/python


def parse_mse(fname):
    ent = [] # list of (tau, entropy) sequence, sep'd by m and r

    data = open(fname,'r')
    data.seek(0)

    for line in data:
        if 'm' in line:
            mr = []  # list of info for each m,r pair
            linesplit = line.split(',')

            l0 = linesplit[0].split('=')
            l1 = linesplit[1].split('=')
            mr.append((int(l0[1].strip()), float(l1[1].strip())))

            newline = data.next()
            entvals = []        
            while 'm' not in newline:
                if newline == '\n':
                    newline = data.next()
                    continue
                ns =newline.split('\t')
                entvals.append((int(ns[0]),float(ns[1])))
                try:
                    newline = data.next()
                except:
                    print 'end of file'
                    break
            mr.append(entvals)
            ent.append(mr)

    data.close()
    return ent
    
    
