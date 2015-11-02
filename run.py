#!/usr/bin/python
import os
As=[0.1*x for x in range(10)]
Ls=["8","12"]
seeds=1
offset=0

for L in Ls:
	os.system("mkdir "+L)
	os.chdir(L)
	os.system("cp ../a.out ../sub.pbs ./")
	os.system("cp ../params temp")

	#edit params file for size
	file=open("params","w")
	input=open("temp","r")
	file.write(L+'\n')
	data=input.readlines()
	for i in range(1,len(data)):
		file.write(data[i])
	
	file.close()
	input.close()
	os.system("rm temp")

	for A in As:
		os.system("mkdir "+str(A))
		os.chdir(str(A))
		os.system("cp ../a.out ../sub.pbs ./")
		os.system("cp ../params temp")

		#edit params file for t
		file=open("params","w")
		input=open("temp","r")
		data=input.readlines()
		for i in range(len(data)):
			if(i==5 or i==6): file.write(str(A)+'\n')
			else: file.write(data[i])
		file.close();
		input.close();
		os.system("rm temp")

		for i in range(seeds):
			os.system("mkdir "+str(i+offset))
			os.chdir(str(i+offset))
			os.system("cp ../a.out ../sub.pbs ./")
			os.system("cp ../params temp")

			#edit params file for seed
			file=open("params","w")
			input=open("temp","r")
			data=input.readlines()
			for j in range(len(data)):
				if(j==4): file.write(str(i+offset)+'\n')
				else: file.write(data[j])
			file.close();
			input.close();
			os.system("rm temp")
			#run job
			os.system("./a.out>rate &")
			os.chdir("..")
		os.chdir("..")
	os.chdir("..")
