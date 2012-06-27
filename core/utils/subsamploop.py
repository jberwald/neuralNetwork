import subsample
import os

srate = 10

slash = "/"
prefix = '/data/jberwald/neurons/10node/gamma/sigma0.2/tfinal2e5/g1/'
suffix = 'avgvolt.trunc.100.out'

dl = os.listdir( prefix )
       
volt_files = [ prefix + x + slash + suffix
               for x in dl if x.startswith('n')]

#print volt_files
print "sampling at rate", srate
for data in volt_files:
    print "working on", data
    print "sampling data..."
    try:
        sdata = subsample.subsample_data( data, srate )
    except:
        print "NO AVGVOLT FILE: ", data
        continue

    print "running dfa..."
    subsample.run_dfa( sdata )

    print "running mse..."
    subsample.run_mse( sdata )

    print ""
