#!/usr/bin/python

import os, time
from subprocess import call



count = 1

while 1:
	import pdb
	# pdb.set_trace()
	time.sleep(1)
	files_list = os.listdir("./daemon.watch")
	if len(files_list) > 0:
		print(files_list)
		job_name = files_list[0].split(".txt")[0]

		#run pipeline
		call(["./main.R", "-n", job_name, "-f", "./daemon.watch/" + files_list[0]])	

		os.remove("./daemon.watch/" + job_name)

	
