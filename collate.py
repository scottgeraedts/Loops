#!/usr/bin/python
import os
import string
import math
pi=3.141592654
Ls=[8,12]
n=1;
for Lnum in Ls:
	L=str(Lnum)
	os.chdir(L);
	Lfile=open("L"+L,"w")
	Lfile.write("t1 t2 E C C3 M chi GG RRx RRyz\n")
	ts=os.listdir(".")
	ts.sort()
	for t in ts:
		try:
			os.chdir(t)
		except OSError:
			continue
		final=[]
		for i in range(6):
			final.append(0)
		j=0
		seeds=os.listdir(".")
		for seed in seeds:
			try:
				os.chdir(seed)
			except OSError:
				continue	

			print os.getcwd()
			
			#calculate the derived quantities and put them in out
			dat=open("out","r")
			out=open("temp","w")
			lines=dat.readlines();
			for line in lines:
				if(line[0]=='#'): continue
				rawdata=line.split()

			rawdata=map(float,rawdata)
			for element in rawdata:
				out.write('%9.7f ' % element)
			
			#compute derived quantities
			out.close()
			dat.close()
			dat=open('temp','r')
			lines=dat.readlines()
			
			#compute averages
			for line in lines:
				if(line[0]=="#"): continue
				temp=line.split()
				for i in range(len(temp)):
					final[i]=final[i]+float(temp[i])
			dat.close()
			os.chdir("..")
			j=j+1
		#print the simple averaged quantities
		for i in range(len(final)):
			final[i]=final[i]/j
			Lfile.write('%7.5f  ' % final[i])

		Lfile.write('\n')
		os.chdir("..")
	Lfile.close()
	os.system("mv L"+L+" ..")
	os.chdir("..")
