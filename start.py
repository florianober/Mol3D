#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
#~ import math
import sys
import time
import random
import numpy as np

import threading 
import queue 

path_home = './'
path_results = path_home+'results/'
qu = queue.Queue()

if 'vis' in sys.argv:
    vis = True
else:
    vis = False
    

random.seed(101)
    
class worker(threading.Thread): 
    Ergebnis = {} 
    ErgebnisLock = threading.Lock() 
 
    def __init__(self, queue_in):
        threading.Thread.__init__(self)
        self.queue = queue_in
        
    def run(self): 
        while True: 
            para_in = self.queue.get() 
            self.ErgebnisLock.acquire()
            # generate input parameter and start Mol3d
            #
            
            self.Ergebnis[para_in[0]] = 'working'
            
            input_parameter = para_in[1]
            set_input_file(input_parameter)
            self.ErgebnisLock.release()
            
            time.sleep(1.2)
            #os.system('./mol3d input/input_'+threading.current_thread().name+'.dat > %s.dat' %(para_in[0]))
            self.ErgebnisLock.acquire() 
            self.Ergebnis[para_in[0]] = 'done'
            self.ErgebnisLock.release() 
 
            self.queue.task_done()


def set_input_file(paralist):
    #~ print('this will be implemented soon')
    print(threading.current_thread().name)
    if paralist == 'dummy':
        # just use the unchanged input file -> "input/input.dat"
        pass
    elif isinstance(paralist,dict):
        
        #~ load most recent input file
        f = open(path_home+"input/input.dat")
        lines = f.readlines()
        f.close()
        
        # write a new input file based on the parameters given in 'paralist'
        # tbd. I should create an class which can handle all input parameter and provide a paralist
        
        f = open(path_home+"input/input_"+threading.current_thread().name+'.dat', 'w')
        for entry in lines:
            val = entry.split()
            if len(val) >=1:
                if 'th_map' == val[0]:
                    entry = 'th_map = {%.1f}                           inclination to observer\n' %(paralist['inc'])
                if 'proname' == val[0]:
                    entry = 'proname = {%s}                       project name\n' %(paralist['pname'])
                if 'old_proname' == val[0]:
                    entry = 'old_proname = {%s}                       project name\n' %(paralist['old_pname'])
                if 'line' == val[0] :
                    entry = 'line = {%d}                                line transition no (see gas_input file)\n' %(paralist['line'])
            f.write(entry)
        f.close()
        
def start_mol3d(vis):
    t1 = time.time()
    os.chdir(path_home)
    # no of threads
    N = 5
    print('')
    print('#        welcome         #')
    print('--------------------------')
    print(' using up to %d cores' %N)
    print(' compile mol3d')
    # make mol3d
    #~ os.system('make CO=fast FC=ifort new >/dev/null')
    if os.path.isfile('mol3d'):
        print(' executable succesfully created')
    else:
        sys.exit('could not find the executable, maybe Mol3d is not compiled (correctly)?')
    # create list with parameters for the run
    inc         = [60,45,10,0]
    tr          = [2,7]
    old_pname   = ['temp_mcgap','temp_mcnogap']
    arr = [] # more general tbd
    counter = 0
    
    for op in range(len(old_pname)):
        for i in range(len(inc)):
            for t in range(len(tr)):
                para_list = {}
                if 'no' in old_pname[op]:
                    para_list['pname'] = '2to9nogapi%.1f_CO%d' %(inc[i],tr[t])
                else:
                    para_list['pname'] = '2to9gapi%.1f_CO%d' %(inc[i],tr[t])
                
                para_list['old_pname'] = old_pname[op]
                para_list['inc']   = inc[i]
                para_list['line']  = tr[t]
                arr.append(para_list)
    
    my_threads = [worker(qu) for i in range(N)] 
    for thread in my_threads: 
        thread.setDaemon(True) 
        thread.start() 

    for g in range(len(arr)):  # put jobs into the queue

        parameter = [arr[g]['pname'],arr[g]]
        #~ parameter = [arr[g]['pname'],'dummy']
        worker.ErgebnisLock.acquire() 
        worker.Ergebnis[parameter[0]] = "queued" 
        worker.ErgebnisLock.release()

        qu.put(parameter)
        # wait a bit before submitting a new job to the queue
        # This is mainly a workaround to let the code start one instance of Mol3d at this moment
        #~ no = random.random()*2.1
        #~ time.sleep(no) 
  
    # some visualisation if wanted
    
    if vis:
        while True:
            print(' for some informations press "r", to quit press "x"')
            eingabe = input(" > ") 
            if eingabe == "x": 
                break 
        
            elif eingabe == "r": 
                print("  -------- current status --------")
                worker.ErgebnisLock.acquire()
                #print(qu.qsize())
                for z, e in worker.Ergebnis.items(): 
                    print(( "  %s: %s" % (z, e))) 
                worker.ErgebnisLock.release() 
                print("  ----------------------------------")
    else:
        print('  working....please wait')
        
    qu.join()
    print(" elapsed total time: %2.2f sec" % (time.time() - t1))
    print('')
    print(' thank you for using Mol3d')
    print('--------------------------')
    print('')
    print('byebye')
    
start_mol3d(vis)
