#/usr/bin/env python
from __future__ import print_function
import os, sys, time, threading, argparse
try:
        from subprocess import getstatusoutput
except:
        from commands import getstatusoutput
global prog_path


class progress_bar_loading(threading.Thread):

        def run(self):
                global stop, kill
                sys.stdout.flush()
                i = 0
                while stop != True:
                        if (i%4) == 0: 
                                sys.stdout.write('\b/')
                        elif (i%4) == 1:
                                sys.stdout.write('\b-')
                        elif (i%4) == 2:
                                sys.stdout.write('\b\\')
                        elif (i%4) == 3:
                                sys.stdout.write('\b|')
                        sys.stdout.flush()
                        time.sleep(0.2)
                        i+=1
                if kill == True:
                        print('\b Abort!', end=' ')
                        sys.exit()
                else:
                        print('\b ', end=' ')
                
                        

def task_1():
        cmd='cd '+prog_dir+'/utils/;'
        cmd=cmd+'git clone https://github.com/soedinglab/hh-suite.git;'
        cmd=cmd+'mkdir -p hh-suite/build && cd hh-suite/build;'
        cmd=cmd+'cmake -DCMAKE_INSTALL_PREFIX=.. ..;'
        cmd=cmd+'make -j 4 && make install'
        out=getstatusoutput(cmd)
        if out[0]!=0:
                print('ERROR: hhblits not installed\n',out[1], file=sys.stderr)
                os._exit(1)
                #sys.exit(1)
        cmd='cd '+prog_dir+'/utils/hh-suite/bin;pwd;./hhblits -h'
        out=getstatusoutput(cmd)
        if out[0]!=0:
                print('\b ')
                print('ERROR: Incorrect hhblits installation\n',out[1], file=sys.stderr)
                os._exit(2)
                #sys.exit(2)
        print('\b   done!')
        return


def task_2():
        www_uc30='http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz'
        file_uc30='uniclust30_2018_08_hhsuite.tar.gz'
        cmd='cd '+prog_dir+'/data/;'
        cmd=cmd+'wget -q '+www_uc30+';'
        cmd=cmd+'ls '+file_uc30
        out=getstatusoutput(cmd)
        if out[0]!=0:
                print('ERROR: uniclust30_2018_08 not downloaded\n',out[1], file=sys.stderr)
                os._exit(3)
                #sys.exit(3)
        print('\b   done!')
        return


def task_3():
        file_uc30='uniclust30_2018_08_hhsuite.tar.gz'
        cmd='cd '+prog_dir+'/data/;'
        cmd=cmd+'tar -xzvf '+file_uc30
        out=getstatusoutput(cmd)
        if out[0]!=0:
                print('ERROR: untar uniclust30_2018_08 not completed\n',out[1], file=sys.stderr)
                os._exit(4)
                #sys.exit(4)
        else:
                cmd='cd '+prog_dir+'/data/;rm uniclust30_2018_08_hhsuite.tar.gz'
                getstatusoutput(cmd)
        print('\b   done!')
        return



def main():
        global prog_dir, stop, kill
        db=False
        parser = argparse.ArgumentParser(description='Program for the installation of DDGun.')
        parser.add_argument ("-d", "--db",action="store_true",dest="db", help="Install DB only")
        args = parser.parse_args()
        if args.db: db=True
        prog_dir = os.path.dirname(os.path.abspath(__file__))
        kill = False
        stop = False
        p = progress_bar_loading()
        p.start()
        c=1
        try:
                if not db:
                        print(str(c)+') Install hhblits', file=sys.stderr)        
                        task_1()
                        c=c+1
                print(str(c)+') Download uniclust30_2018_08 (25Gb)', file=sys.stderr)
                task_2()
                c=c+1
                print(str(c)+') Untar uniclust30_2018_08 (25Gb)', file=sys.stderr)
                task_3()
                time.sleep(1)
                stop=True
        except KeyboardInterrupt or EOFError:
                kill = True
                stop = True
        


if __name__ == '__main__':
        main()
        
