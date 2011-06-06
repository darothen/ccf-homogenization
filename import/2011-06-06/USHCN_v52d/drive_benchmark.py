#!/usr/bin/env python

import os
import subprocess
import logging
import optparse

logging.basicConfig(level=logging.INFO)

def main():
    
    ## Set-up the command line options. Not super-necessary right now,
    ## but might be useful later when broadening how the program works.
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("-u", "--ucp", dest="ucp",
                      action="store", type="string",
                      default="UCPM24aM21fC6bL1Si6cR5aH1")
    parser.add_option("-m", "--mode", dest="mode",
                      action="store", type="string",
                      default="FAST")
    parser.add_option("-c", "--case", dest="case",
                      action="store", type="string",
                      default="21_21")
    (options, args) = parser.parse_args()
    options_dict = options.__dict__
    
    print "Executing USHCN_v52d Driver Program\n"
    print "   Options:"
    for option, value in options_dict.iteritems():
        print "    %s - %s" % (option, value)
    print "\n"
    
    ## Begin setting up environmental variables
    USHCNBASE = os.getcwd()
    UCP = options.ucp
    MODE = options.mode
    CASE = options.case
    
    os.environ["USHCNBASE"] = USHCNBASE # necessary environmental variable
                                        # for the actual USHCN scripts
    print "$USHCNBASE is: %s" % USHCNBASE
                                          
    ## Create a 'bin' top-level directory to house the executable
    print "Creating $USHCNBASE/bin"
    bin_dir = os.path.join(USHCNBASE, "bin")
    if not os.path.exists(bin_dir):
        os.makedirs(bin_dir)
    else:
        old_executable = "%s.%s.MLY.%s" % (UCP, MODE, CASE)
        if os.path.exists(os.path.join(bin_dir, old_executable)):
            os.remove(os.path.join(bin_dir, old_executable))
            print "   Removed old executable: %s" % old_executable
    
    ## Call the compile script
    print "\nCompiling USHCN_v52d source code"
    print "--------------------------------"
    subprocess.call(["./src_codes/lib/scripts/compile.csh",
                     MODE,
                     CASE,
                     USHCNBASE])
    print "--------------------------------"
    
    ## Call the actual benchmark script
    print "\nRunning the USHCN_v52d program"
    print "--------------------------------"
    subprocess.call(["./src_codes/benchmark/scripts/case7_bench.sh",
                     UCP,
                     MODE,
                     CASE,
                     USHCNBASE])
    print "--------------------------------"
    
    ## All done!
    print "\n done."                              
   
if __name__ == "__main__":
    main()

    