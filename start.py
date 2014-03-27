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
            self.Ergebnis[para_in[0]] = 'working'
            self.ErgebnisLock.release()
            
            no = random.random()*10
            print('sleeping for %2.2f seconds' %no)
            time.sleep(no)
            
            
            self.ErgebnisLock.acquire() 
            self.Ergebnis[para_in[0]] = 'done'
            self.ErgebnisLock.release() 
 
            self.queue.task_done()
        
def start_mol3d(vis):
    t1 = time.time()
    os.chdir(path_home)
    # no of threads
    N = 2
    
    arr = np.arange(20)
    my_threads = [worker(qu) for i in range(N)] 
    for thread in my_threads: 
        thread.setDaemon(True) 
        thread.start() 

                         
    for g in range(len(arr)):  # put jobs into the queue
        
        pname = 'test_%d' %g
        
        parameter = [pname,arr[g]]
        worker.ErgebnisLock.acquire() 
        worker.Ergebnis[parameter[0]] = "queued" 
        worker.ErgebnisLock.release()

        qu.put(parameter) 
        # wait a bit
    # some visualisation if wanted
    
    if vis:
        while True: 
            eingabe = input("> ") 
            if eingabe == "ende": 
                break 
        
            elif eingabe == "status": 
                print("-------- Aktueller Status --------")
                worker.ErgebnisLock.acquire()
                print(qu.qsize())
                for z, e in worker.Ergebnis.items(): 
                    print(( "%s: %s" % (z, e))) 
                worker.ErgebnisLock.release() 
                print("----------------------------------")

    qu.join()
    print("Elapsed Time: %2.2f sec" % (time.time() - t1))
    
start_mol3d(vis)
