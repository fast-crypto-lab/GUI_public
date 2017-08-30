#!/usr/bin/python

import os
import shutil
import fileinput

tar = [ "quartzref", "quartz128ref", "quartz103129", "quartz9417", "quartz959", "ecdonaldp160", "ecdonaldp192", "ecdonaldp256", "ronald1024", "ronald2048" ]
tar2 = [ "cycles" , "open_cycles" ]

def parse_line(line):
  items = line.split(" ")
  count = 0
  accu = 0
  #any( k in s for k in keywords )
  if ( any( items[5]==s for s in tar ) and any( items[6]==s for s in tar2 ) and items[7]=="59" ) :
    #print line
    for i in items[8:]:
      if( i.isdigit() ):
        count = count+1;
        accu = accu + int(i)
    print items[5],items[6],items[7],"\t[",str(count),"]","\t(",str(accu/count),")"
       

for line in fileinput.input():
  parse_line(line)



