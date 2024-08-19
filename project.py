import math
import sys

class Process:
	def __init__(self):
		self.id = -1
		self.arrival_time = -1
		self.CPU_burst_times = []
		self.IO_burst_times = []
		self.turnaround_times = []
		self.wait_times = []
		self.current_burst_no = 0
		self.tau = -1
		self.tau_remaining = -1

'''
	|----|----|----|
	| A0 | A1 | A2 | ...
	|----|----|----|
	| B0 | B1 | B2 |
	|----|----|----|
	...
'''

processes_table = []
ROWS = 26
COLS = 10
x_n = 0

def srand48(seed):
	global x_n
	x_n = (seed << 16) + 0x330E
    
def drand48():
	global x_n
	a = 0x5DEECE66D
	c = 0xB
	m = pow(2,48)
	x_n_plus_1 = (a * x_n + c) % m
	x_n = x_n_plus_1
	return x_n_plus_1/m

def next_exp(lmbda, bound):
	r = drand48()
	x = -1 * math.log(r) / lmbda
	while True:
		if (math.ceil(x) > bound):
			r = drand48()
			x = -1 * math.log(r) / lmbda
		else:
			break
	return x

def flattenProcessTable():
	ret = []
	for i in range(len(processes_table)):
		for j in range(len(processes_table[i])):
			if (processes_table[i][j].id == -1):
				return ret
			ret.append(processes_table[i][j])
	return ret

def parseQueue(queue):
	if (len(queue) == 0): return "[Q empty]"
	out = "[Q"
	for p in queue:
		out += " " + p.id
	out += "]"
	return out

def defaultAllProcesses():
	for i in range(ROWS):
		for j in range(COLS):
			processes_table[i][j].current_burst_no = 0
			processes_table[i][j].tau = math.ceil(1 / lmbda)
			processes_table[i][j].tau_remaining = math.ceil(1 / lmbda)

def FCFS(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	current_burst = [-1,-1,-1,-1]		# [Process object, When will it finish switching in?, When will it complete?, When will it finish switching out?]

	ready_queue = []
	blocked_processes = []
	using_CPU = False

	print("time 0ms: Simulator started for FCFS [Q empty]")

	while (True):

		# Check if any process finished blocking on I/O
		if (blocked_processes):
			for blocked in blocked_processes:
				if (time == blocked[1]): 		# blocked process finished blocking?
					process = blocked[0]
					ready_queue.append(process)
					if (time < 10000): 
						print("time {:d}ms: Process {} completed I/O; added to ready queue {}"
							.format(time, process.id, parseQueue(ready_queue)))
					blocked_processes.remove(blocked)

		# Has process finished switching in? If yes, we can start the process
		if (time == current_burst[1]):
			if (time < 10000): 
				process = current_burst[0]
				process_burst_time = current_burst[2] - current_burst[1]
				print("time {:d}ms: Process {} started using the CPU for {:d}ms burst {}"
					.format(time, process.id, process_burst_time, parseQueue(ready_queue)))
				
		# Has process finished switching out? If yes, we can block on I/O and mark the CPU ready
		if (using_CPU and time == current_burst[3]):
			process = current_burst[0]
			if (process.current_burst_no == len(process.CPU_burst_times)):
				processes.remove(process)
			else:
				blocked = [-1, -1]			# [Process object, when will it be unblocked?]
				blocked[0] = process
				blocked[1] = time + t_cs // 2 + process.IO_burst_times[process.current_burst_no - 1]
				blocked_processes.append(b)
				current_burst = [-1,-1,-1,-1]	# Reset current CPU burst data
			using_CPU = False

		# Has process completed? If yes, we should start switching it out
		if (using_CPU):
			if (time == current_burst[2]):
				process = current_burst[0]
				process.current_burst_no += 1

				if (process.current_burst_no == len(process.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(time, process.id, parseQueue(ready_queue)))
				else:
					bursts_left = len(process.CPU_burst_times) - process.current_burst_no
					if (bursts_left == 1):
						if (time < 10000): 
							print("time {:d}ms: Process {} completed a CPU burst; 1 burst to go {}"
			 					.format(time, process.id, parseQueue(ready_queue)))
					else:
						if (time < 10000): 
							print("time {:d}ms: Process {} completed a CPU burst; {:d} bursts to go {}"
								.format(time, process.id, bursts_left, parseQueue(ready_queue)))

					b = [-1, -1]		# [Process object, when will it be unblocked?]
					b[0] = process
					b[1] = time + t_cs // 2 + process.IO_burst_times[process.current_burst_no - 1]
					blocked_processes.append(b)

					if (time < 10000): 
						print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}"
							.format(time, process.id, b[1], parseQueue(ready_queue)))

		# Check if any process has arrivied
		for process in processes:							
			if (time == process.arrival_time):
				ready_queue.append(process)
				if (time < 10000): 
					print("time {:d}ms: Process {} arrived; added to ready queue {}"
		   				.format(time, process.id, parseQueue(ready_queue)))

		# If not using CPU and there is a process in the ready queue, we should start switching it in
		if (not using_CPU and ready_queue):		
			process = ready_queue.pop(0)			
			process_burst_time = process.CPU_burst_times[process.current_burst_no]
			current_burst[0] = process												# Process object
			current_burst[1] = time + t_cs // 2										# When will it finish switching in?
			current_burst[2] = time + t_cs // 2 + process_burst_time				# When will it complete?
			current_burst[3] = time + t_cs // 2 + process_burst_time + t_cs // 2	# When will it finish switching out?
			using_CPU = True

		if (not processes):
			print("time {:d}ms: Simulator ended for FCFS [Q empty]".format(time))
			break
		
		time += 1

def SJF(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	current_burst = [-1,-1,-1,-1]		# [Process object, When will it finish switching in?, When will it complete?, When will it finish switching out?]

	ready_queue = []
	blocked_processes = []
	using_CPU = False

	print("time 0ms: Simulator started for SJF [Q empty]")

	while (True):

		# Check if any process finished blocking on I/O
		if (blocked_processes):
			for blocked in blocked_processes:
				if (time == blocked[1]): 		# blocked process finished blocking?
					process = blocked[0]
					ready_queue.append(process)
					ready_queue.sort(key=lambda x: (x.tau, x.id[0], x.id[1]))
					if (time < 10000): 
						print("time {:d}ms: Process {} (tau {:d}ms) completed I/O; added to ready queue {}"
							.format(time, process.id, process.tau, parseQueue(ready_queue)))
					blocked_processes.remove(blocked)

		# Has process finished switching in? If yes, we can start the process
		if (time == current_burst[1]):
			if (time < 10000): 
				process = current_burst[0]
				process_burst_time = current_burst[2] - current_burst[1]
				print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}"
					.format(time, process.id, process.tau, process_burst_time, parseQueue(ready_queue)))
				
		# Has process finished switching out? If yes, we can block on I/O and mark the CPU ready
		if (using_CPU and time == current_burst[3]):
			process = current_burst[0]
			if (process.current_burst_no == len(process.CPU_burst_times)):
				processes.remove(process)
			else:
				blocked = [-1, -1]			# [Process object, when will it be unblocked?]
				blocked[0] = process
				blocked[1] = time + t_cs // 2 + process.IO_burst_times[process.current_burst_no - 1]
				blocked_processes.append(b)
				current_burst = [-1,-1,-1,-1]	# Reset current CPU burst data
			using_CPU = False


		# Has process completed? If yes, we should start switching it out
		if (using_CPU):
			if (time == current_burst[2]):
				process = current_burst[0]
				process.current_burst_no += 1
				process_old_tau = process.tau

				if (process.current_burst_no == len(process.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(time, process.id, parseQueue(ready_queue)))
				else:
					bursts_left = len(process.CPU_burst_times) - process.current_burst_no
					if (bursts_left == 1):
						if (time < 10000): 
							print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; 1 burst to go {}"
			 					.format(time, process.id, process.tau, parseQueue(ready_queue)))
					else:
						if (time < 10000): 
							print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; {:d} bursts to go {}"
								.format(time, process.id, process.tau, bursts_left, parseQueue(ready_queue)))
							
					process.tau = math.ceil(process.CPU_burst_times[process.current_burst_no-1] * alpha + (1 - alpha) * process.tau)
					if (time < 10000): print("time {:d}ms: Recalculated tau for process {}: old tau {:d}ms ==> new tau {:d}ms {}"
							.format(time, process.id, process_old_tau, process.tau, parseQueue(ready_queue)))

					b = [-1, -1]		# [Process object, when will it be unblocked?]
					b[0] = process
					b[1] = time + t_cs // 2 + process.IO_burst_times[process.current_burst_no - 1]
					blocked_processes.append(b)

					if (time < 10000): 
						print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}"
							.format(time, process.id, b[1], parseQueue(ready_queue)))

		# Check if any process has arrivied
		for process in processes:							
			if (time == process.arrival_time):
				ready_queue.append(process)
				ready_queue.sort(key=lambda x: (x.tau, x.id[0], x.id[1]))
				if (time < 10000): 
					print("time {:d}ms: Process {} (tau {:d}ms) arrived; added to ready queue {}"
		   				.format(time, process.id, process.tau, parseQueue(ready_queue)))
					
		# If not using CPU and there is a process in the ready queue, we should start switching it in
		if (not using_CPU and ready_queue):		
			process = ready_queue.pop(0)			
			process_burst_time = process.CPU_burst_times[process.current_burst_no]
			current_burst[0] = process												# Process object
			current_burst[1] = time + t_cs // 2										# When will it finish switching in?
			current_burst[2] = time + t_cs // 2 + process_burst_time				# When will it complete?
			current_burst[3] = time + t_cs // 2 + process_burst_time + t_cs // 2	# When will it finish switching out?
			using_CPU = True

		if (not processes):
			print("time {:d}ms: Simulator ended for SJF [Q empty]".format(time))
			break
		
		time += 1

def SRT(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	current_burst = [-1,-1,-1,-1,-1]		
	# [ Process object,
	#   When will it finish switching in?,
	#   When will it complete?, 
	#   When will it finish switching out?
	#   Preempted? ]

	ready_queue = []
	blocked_processes = []
	preempted_processes = []
	using_CPU = False

	print("time 0ms: Simulator started for SRT [Q empty]")

	while (True):

		# Check if any process finished blocking on I/O
		if (blocked_processes):
			for blocked in blocked_processes:
				if (time == blocked[1]): 		# blocked process finished blocking?
					process = blocked[0]

					current_burst[4] = 0
					ready_queue.append(process)

					# Preemption?
					if (using_CPU and process.tau < current_burst[0].tau_remaining):
						ready_queue.sort(key=lambda x: (x.tau_remaining, x.tau, x.id[0], x.id[1]))
						c = current_burst[0]

						c.tau_remaining += t_cs // 2		# <----- ?????????????
					

						# Switch current process out	(preempted)
						# ready_queue.insert(0, c)
						preempted = [-1, -1]						# [Process object, time remaining]
						preempted[0] = c
						preempted[1] = current_burst[1] - time
						preempted_processes.append(preempted)

						# Switch blocked process in
						process_burst_time = process.CPU_burst_times[process.current_burst_no]
						current_burst[0] = process															# Process object
						current_burst[1] = time + t_cs // 2	+ t_cs // 2										# When will it finish switching in?
						current_burst[2] = time + t_cs // 2 + process_burst_time + t_cs // 2				# When will it complete?
						current_burst[3] = time + t_cs // 2 + process_burst_time + t_cs // 2 + t_cs // 2	# When will it finish switching out?
						if (time < 10000):
							print("time {:d}ms: Process {} (tau {:d}ms) completed I/O; preempting {} (predicted remaining time {:d}ms) {}"
								.format(time, process.id, process.tau, c.id, c.tau_remaining, parseQueue(ready_queue)))

						ready_queue.insert(0, c)
						ready_queue.remove(process)
						
					else:
						ready_queue.sort(key=lambda x: (x.tau_remaining, x.tau, x.id[0], x.id[1]))
						if (time < 10000):
							print("time {:d}ms: Process {} (tau {:d}ms) completed I/O; added to ready queue {}"
								.format(time, process.id, process.tau, parseQueue(ready_queue)))
								
					blocked_processes.remove(blocked)

		# Has process finished switching in? If yes, we can start the process
		if (time == current_burst[1]):
			if (time < 10000): 
				process = current_burst[0]

				if (current_burst[4] == 0):
					process_burst_time = current_burst[2] - current_burst[1]
					print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}"
						.format(time, process.id, process.tau, process_burst_time, parseQueue(ready_queue)))
				else:
					process_remaining_time = current_burst[2] - current_burst[1]
					process_burst_time = process.CPU_burst_times[process.current_burst_no]
					print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for remaining {:d}ms of {:d}ms burst {}"
		   				.format(time, process.id, process.tau, process_remaining_time, process_burst_time, parseQueue(ready_queue)))

		# Has process finished switching out? If yes, we can block on I/O and mark the CPU ready
		if (using_CPU and time == current_burst[3]):
			process = current_burst[0]
			if (process.current_burst_no == len(process.CPU_burst_times)):
				processes.remove(process)
			else:
				blocked = [-1, -1]			# [Process object, when will it be unblocked?]
				blocked[0] = process
				blocked[1] = time + t_cs // 2 + process.IO_burst_times[process.current_burst_no - 1]
				blocked_processes.append(b)
				current_burst = [-1,-1,-1,-1,-1]	# Reset current CPU burst data
			using_CPU = False

		# Has process completed? If yes, we should start switching it out
		if (using_CPU):
			if (time == current_burst[2]):
				process = current_burst[0]
				process.current_burst_no += 1
				process_old_tau = process.tau

				if (process.current_burst_no == len(process.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(time, process.id, parseQueue(ready_queue)))
				else:
					bursts_left = len(process.CPU_burst_times) - process.current_burst_no
					if (bursts_left == 1):
						if (time < 10000): 
							print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; 1 burst to go {}"
			 					.format(time, process.id, process.tau, parseQueue(ready_queue)))
					else:
						if (time < 10000): 
							print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; {:d} bursts to go {}"
								.format(time, process.id, process.tau, bursts_left, parseQueue(ready_queue)))
					
					process_new_tau = math.ceil(process.CPU_burst_times[process.current_burst_no-1] * alpha + (1 - alpha) * process.tau)
					process.tau = process_new_tau
					process.tau_remaining = process_new_tau

					if (time < 10000): 
						print("time {:d}ms: Recalculated tau for process {}: old tau {:d}ms ==> new tau {:d}ms {}"
							.format(time, process.id, process_old_tau, process.tau, parseQueue(ready_queue)))

					b = [-1, -1]		# [Process object, when will it be unblocked?]
					b[0] = process
					b[1] = time + t_cs // 2 + process.IO_burst_times[process.current_burst_no - 1]
					blocked_processes.append(b)

					if (time < 10000): 
						print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}"
							.format(time, process.id, b[1], parseQueue(ready_queue)))

		# Check if any process has arrivied
		for process in processes:							
			if (time == process.arrival_time):
				ready_queue.append(process)
				ready_queue.sort(key=lambda x: (x.tau_remaining, x.tau, x.id[0], x.id[1]))
				if (time < 10000): 
					print("time {:d}ms: Process {} (tau {:d}ms) arrived; added to ready queue {}"
		   				.format(time, process.id, process.tau, parseQueue(ready_queue)))
					
		# If not using CPU and there is a process in the ready queue, we should start switching it in
		if (not using_CPU and ready_queue):		
			resume_preempted = False
			process = ready_queue.pop(0)			

			for preempted in preempted_processes:
				if (process == preempted[0]):
					resume_preempted = True
					process_burst_time = process.CPU_burst_times[process.current_burst_no] - (process.tau - process.tau_remaining)
					current_burst[0] = process													# Process object
					current_burst[1] = time + t_cs // 2											# When will it finish switching in?
					current_burst[2] = time + t_cs // 2 + process_burst_time					# When will it complete?
					current_burst[3] = time + t_cs // 2 + process_burst_time + t_cs // 2		# When will it finish switching out?
					current_burst[4] = 1
					preempted_processes.remove(preempted)
					break

			if (not resume_preempted):
				process_burst_time = process.CPU_burst_times[process.current_burst_no]
				current_burst[0] = process												# Process object
				current_burst[1] = time + t_cs // 2										# When will it finish switching in?
				current_burst[2] = time + t_cs // 2 + process_burst_time				# When will it complete?
				current_burst[3] = time + t_cs // 2 + process_burst_time + t_cs // 2	# When will it finish switching out?
				current_burst[4] = 0
			using_CPU = True

		if (not processes):
			print("time {:d}ms: Simulator ended for SRT [Q empty]".format(time))
			break
		
		if (time >= current_burst[1] and current_burst[0] != -1):	# track remaining tau of current burst
			current_burst[0].tau_remaining -= 1

		time += 1

def RR(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	current_burst = [-1,-1,-1,-1]		# [Process object, When will it finish switching in?, When will it complete?, When will it finish switching out?]

	ready_queue = []
	blocked_processes = []
	using_CPU = False

	print("time 0ms: Simulator started for RR [Q empty]")


if (__name__ == "__main__"):

	args = sys.argv
	if (len(args) != 9):
		sys.stderr.write("ERROR: invalid number of arguments\n")
		sys.exit(1)
	
	n = int(args[1])
	n_CPU = int(args[2])
	seed = int(args[3])
	srand48(seed)
	lmbda = float(args[4])
	bound = int(args[5])

	t_cs = int(args[6])
	alpha = float(args[7])
	t_slice = int(args[8])

	if (n <= 0):
		sys.stderr.write("ERROR: invalid number of processes\n")
		sys.exit(1)

	if (n_CPU < 0 or n_CPU > n):
		sys.stderr.write("ERROR: invalid number of processes that are CPU-bound\n")
		sys.exit(1)

	if (lmbda <= 0):
		sys.stderr.write("ERROR: invalid lambda value\n")
		sys.exit(1)
	
	if (bound <= 0):
		sys.stderr.write("ERROR: invalid bound value\n")
		sys.exit(1)

	for i in range(ROWS):
		row = []
		for j in range(COLS):
			row.append(Process())
		processes_table.append(row)

	type = 0

	CPU_total_CPU_burst_time = 0.0
	CPU_total_IO_burst_time = 0.0
	CPU_num_CPU_burst = 0
	CPU_num_IO_burst = 0

	IO_total_CPU_burst_time = 0.0
	IO_total_IO_burst_time = 0.0
	IO_num_CPU_burst = 0
	IO_num_IO_burst = 0
	
	print("<<< PROJECT PART I", flush=True)
	if (n_CPU == 1):
		print("<<< -- process set (n={:d}) with 1 CPU-bound process".format(n))
	else:
		print("<<< -- process set (n={:d}) with {:d} CPU-bound processes".format(n, n_CPU))
	
	print("<<< -- seed={:d}; lambda={:.6f}; bound={:d}".format(seed, lmbda, bound))

	for i in range(n):
		if (i >= n_CPU):
			type = 1
		
		letter = 65 + i // 10
		number = i % 10

		row = int(letter - 65)
		col = int(number)

		processes_table[row][col].tau = math.ceil(1 / lmbda)
		processes_table[row][col].tau_remaining = math.ceil(1 / lmbda)

		processes_table[row][col].id = str(chr(letter)) + str(number)


		arrival_time = math.floor(next_exp(lmbda, bound))
		processes_table[row][col].arrival_time = arrival_time

		uniform_distribution = drand48()
		CPU_bursts = math.ceil(uniform_distribution * 32)

		if (type == 0):
			if (CPU_bursts == 1):
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst".format(letter, number, arrival_time))
			else:
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU bursts".format(letter, number, arrival_time, CPU_bursts))
		else:
			if (CPU_bursts == 1):
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst".format(letter, number, arrival_time))
			else:
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU bursts".format(letter, number, arrival_time, CPU_bursts))

		for j in range(CPU_bursts):
			CPU_burst_time = math.ceil(next_exp(lmbda, bound))

			IO_burst_time = -1
			if (j < CPU_bursts - 1):
				IO_burst_time = math.ceil(next_exp(lmbda, bound))

			if (type == 0):
				CPU_burst_time *= 4
				CPU_total_CPU_burst_time += CPU_burst_time
				CPU_num_CPU_burst += 1
				if (IO_burst_time != -1):
					CPU_total_IO_burst_time += IO_burst_time
					CPU_num_IO_burst += 1
					processes_table[row][col].IO_burst_times.append(IO_burst_time)
			else:
				IO_burst_time *= 8
				IO_total_CPU_burst_time += CPU_burst_time
				IO_num_CPU_burst += 1
				if (IO_burst_time != -8):
					IO_total_IO_burst_time += IO_burst_time
					IO_num_IO_burst += 1
					processes_table[row][col].IO_burst_times.append(IO_burst_time)

			processes_table[row][col].CPU_burst_times.append(CPU_burst_time)

			# if (j == CPU_bursts-1):
			# 	print("==> CPU burst {:d}ms".format(CPU_burst_time))
			# else:
			# 	print("==> CPU burst {:d}ms ==> I/O burst {:d}ms".format(CPU_burst_time, IO_burst_time))

	print()
	print("<<< PROJECT PART II")
	print("<<< -- t_cs={:d}ms; alpha={:.2f}; t_slice={:d}ms".format(t_cs, alpha, t_slice))

	FCFS(t_cs, alpha, t_slice)
	print()
	SJF(t_cs, alpha, t_slice)
	print()
	SRT(t_cs, alpha, t_slice)
	print()
	RR(t_cs, alpha, t_slice)

	CPU_avg_CPU_burst_time = 0
	if (CPU_num_CPU_burst != 0):
		CPU_avg_CPU_burst_time = CPU_total_CPU_burst_time / CPU_num_CPU_burst
	
	IO_avg_CPU_burst_time = 0
	if (IO_num_CPU_burst != 0):
		IO_avg_CPU_burst_time = IO_total_CPU_burst_time / IO_num_CPU_burst
	
	avg_CPU_burst_time = 0
	if (CPU_num_CPU_burst + IO_num_CPU_burst != 0):
		avg_CPU_burst_time = (CPU_total_CPU_burst_time + IO_total_CPU_burst_time) / (CPU_num_CPU_burst + IO_num_CPU_burst)

	CPU_avg_IO_burst_time = 0
	if (CPU_num_IO_burst != 0):
		CPU_avg_IO_burst_time = CPU_total_IO_burst_time / CPU_num_IO_burst
	
	IO_avg_IO_burst_time = 0
	if (IO_num_IO_burst != 0):
		IO_avg_IO_burst_time = IO_total_IO_burst_time / IO_num_IO_burst

	avg_IO_burst_time = 0
	if (CPU_num_IO_burst + IO_num_IO_burst != 0):
		avg_IO_burst_time = (CPU_total_IO_burst_time + IO_total_IO_burst_time) / (CPU_num_IO_burst + IO_num_IO_burst)
	
	f = open("simout.txt", "w")
	f.write("-- number of processes: {:d}\n".format(n))
	f.write("-- number of CPU-bound processes: {:d}\n".format(n_CPU))
	f.write("-- number of I/O-bound processes: {:d}\n".format(n - n_CPU))
	f.write("-- CPU-bound average CPU burst time: {:.3f} ms\n".format(math.ceil(CPU_avg_CPU_burst_time * 1000) / 1000))
	f.write("-- I/O-bound average CPU burst time: {:.3f} ms\n".format(math.ceil(IO_avg_CPU_burst_time * 1000) / 1000))
	f.write("-- overall average CPU burst time: {:.3f} ms\n".format(math.ceil(avg_CPU_burst_time * 1000) / 1000))
	f.write("-- CPU-bound average I/O burst time: {:.3f} ms\n".format(math.ceil(CPU_avg_IO_burst_time * 1000) / 1000))
	f.write("-- I/O-bound average I/O burst time: {:.3f} ms\n".format(math.ceil(IO_avg_IO_burst_time * 1000) / 1000))
	f.write("-- overall average I/O burst time: {:.3f} ms\n".format(math.ceil(avg_IO_burst_time * 1000) / 1000))
	f.close()
