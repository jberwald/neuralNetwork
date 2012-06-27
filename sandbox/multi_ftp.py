"""
Transfer multiple files from simplex -> laptop
"""
import os, sys, time, paramiko

# hardcoded. might add options for conley in future.
username = 'jberwald'
password = 'PIN8diki'
hostname = '153.90.246.232'
hostport = 22

def get(sc, remotefile, localfile):
    """
    Use get() method to get a file from remote host
    T         : Transport instance
    
    remotefile: file to get on remote machine 
    
    localfile : file where to put it (incl. filename)
    [both remote and local expect full path]
    """
    try:
        sc.get(remotefile, localfile)
    except:
        print "failed to get file from remote host"
        sys.exit(1)

def transfer_many(T, remotefolder, localfolder, filelist, 
                  cmd='get', verbose=True):
    """
    Remote and local folders are supersets that contain files from
    filelist. filelist is a static list of files to
    transfer. remotefolderand localfolder are changed as we loop
    through various simulation numbers.

    T   : sftp_client object

    Eg., remotefolder = 'remote/path/n10g0e1c1s0'
         localfolder = 'local/path/n10g0e0c1s0'

         loop formation of folder names spurred by simtype =
         {n,g,e,c}. Must provide number of sims (eg., g0->g9, so
         numsims=10)
    """
    for f in filelist:
        rfile = remotefolder + f
        lfile = localfolder + f

        if cmd == 'get':
            try: 
                get(T, rfile, lfile)
                if verbose:
                    print "transfer of ", f, " successful"
            except:
                T.close()
                print "failed to transfer from remote host"
                sys.exit(1)
        elif cmd == 'put':
            try: 
                put(T, rfile, lfile)
                if verbose:
                    print "transfer of ", f, " successful"
            except:
                print "failed to transfer from local host"
                sys.exit(1)

def goto_fail():
    t.close()


if __name__ == "__main__":
    
    # open ssh transport
    t = paramiko.Transport((hostname, hostport))
    t.connect(username=username, password=password)
    time.sleep(2)

    rdir = '/data/jberwald/tests/fitness_tests/network_topology/'
    ldir = '/home/data/jberwald/tests/fitness_tests/network_topology/'

    dpre = 'n10'
    dsuf = 's0/'

    flist =  ['dfaplot_v.0.0.pdf','findiff_dfa.pdf','ts_binned.v.0.0.pdf','eps','graph.pkl', 'avgvolt.sr.1.v.0.dfa']

    print "opening connections..."
    try:
        sc = t.open_sftp_client()
        time.sleep(2)
    except:
        print "failed to open sftp connection"
        sys.exit(1)

    print "success!"
    print "your transfer object is ", t

    for i in range(10):#[6,7,8,9]:
        for j in [5,6,7,8,9,10]:#[1,2,3]:
            for k in range(1,11):#[8,9,10]:
                print 'g ', i
                print 'e ', j
                print 'c ', k
                folder = dpre + 'g'+str(i)+'e'+str(j)+'c'+str(k)+dsuf
                rfile = rdir+folder
                lfile = ldir+folder
                transfer_many(sc, rfile, lfile, flist)

    try:      
        t.close()
    except:
        t.close()
        print "close() failed somehow"
        sys.exit(1)
    finally: 
        print "done"

    
    
