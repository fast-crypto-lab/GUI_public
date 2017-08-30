#!/usr/bin/python

import os
import shutil
import fileinput

tar = [ "quartz959", "quartz9417", "quartz1279" , "quartz103129", "quartzref", "quartz128ref", "ecdonaldp160", "ecdonaldp192", "ecdonaldp256", "ronald1024", "ronald2048" ]
tar2 = [ "cycles" , "open_cycles" ]

accu = [ [0,0] for i in tar ]

def parse_line(line):
  items = line.split(" ")
  #any( k in s for k in keywords )
  if ( len(items) <= 7 ) : return 
  if ( any( items[5]==s for s in tar ) and any( items[6] == s for s in tar2) and items[7]=="59" ) :
    print line
    tt = [ int(i) for i in items[8:len(items)-1] ]
    idx = tar.index(items[5])
    s_idx = tar2.index(items[6])
    accu[idx][s_idx] += sum(tt)/len(tt)

#print items[5],items[6],items[7],"\t[",str(count),"]","\t(",str(accu/count),")"
       

for line in fileinput.input():
  parse_line(line)

for i,j in enumerate(accu):
  print tar[i],"\t\tcycles: ",str(j[0]/3),"\topen_cycles: ",str(j[1]/3)

