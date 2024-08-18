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
lmbda_c = -1

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
			processes_table[i][j].tau = math.ceil(1 / lmbda_c)
			processes_table[i][j].tau_remaining = math.ceil(1 / lmbda_c)

def FCFS(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	queue = []
	
	using_CPU = False
	current_burst = [-1, -1, -1]		# [Process object, when will it finish?, when will it start?]

	blocked = []

	print("time 0ms: Simulator started for FCFS [Q empty]", flush=True)

	while (processes):
		
		if (time == current_burst[2]):
			p = current_burst[0]
			
			p_burst_time = current_burst[1] - current_burst[2]
			if (time < 10000): 
				print("time {:d}ms: Process {} started using the CPU for {:d}ms burst {}"
					.format(time, p.id, p_burst_time, parseQueue(queue)), flush=True)

		for p in processes:
			if (time == p.arrival_time):
				queue.append(p)
				if (time < 10000): print("time {:d}ms: Process {} arrived; added to ready queue {}".format(time, p.id, parseQueue(queue)), flush=True)

		if (blocked):
			for blocked_p in blocked:
				if (time == blocked_p[1]): 	# blocked process finished blocking?
					p = blocked_p[0]
					queue.append(p)
					if (time < 10000): print("time {:d}ms: Process {} completed I/O; added to ready queue {}".format(time, p.id, parseQueue(queue)), flush=True)
					blocked.remove(blocked_p)
		
		if (using_CPU):
			if (time-t_cs//2 == current_burst[1]):	# current CPU burst completed?
				ctime = time - t_cs//2
				using_CPU = False
				p = current_burst[0]
				current_burst = [-1, -1, -1]	# empty current burst info
				p.current_burst_no += 1

				if (p.current_burst_no == len(p.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(ctime, p.id, parseQueue(queue)), flush=True)
					#time += t_cs // 2
					processes.remove(p)
					continue
				
				bursts_left = len(p.CPU_burst_times) - p.current_burst_no
				if (bursts_left == 1):
					if (time < 10000): print("time {:d}ms: Process {} completed a CPU burst; 1 burst to go {}".format(ctime, p.id, parseQueue(queue)), flush=True)
				else:
					if (time < 10000): print("time {:d}ms: Process {} completed a CPU burst; {:d} bursts to go {}".format(ctime, p.id, bursts_left, parseQueue(queue)), flush=True)

				b = [-1, -1]		# [Process object, when will it be unblocked?]
				b[0] = p
				b[1] = time  + p.IO_burst_times[p.current_burst_no - 1]
				blocked.append(b)

				if (time < 10000): print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}".format(ctime, p.id, b[1], parseQueue(queue)), flush=True)

				# time += t_cs // 2
				# continue
				
		if (not using_CPU and queue):
			using_CPU = True
			p = queue.pop(0)
			p_burst_time = p.CPU_burst_times[p.current_burst_no]
			current_burst[0] = p
			current_burst[1] = time + t_cs // 2 + p_burst_time
			current_burst[2] = time + t_cs // 2
			#time += t_cs // 2
			#if (time < 10000): print("time {:d}ms: Process {} started using the CPU for {:d}ms burst {}".format(time, p.id, p_burst_time, parseQueue(queue)))

		time += 1

	print("time {:d}ms: Simulator ended for FCFS [Q empty]".format(time), flush=True)

def SJF(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	queue = []

	using_CPU = False
	current_burst = [-1, -1, -1]		# [Process object, when will it finish?, when will it start?]

	blocked = []

	print("time 0ms: Simulator started for SJF [Q empty]", flush=True)

	while (processes):

		if (time == current_burst[2]):
			p = current_burst[0]
			
			p_burst_time = current_burst[1] - current_burst[2]
			if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}"
					.format(time, p.id, p.tau, p_burst_time, parseQueue(queue)), flush=True)

		for p in processes:
			if (time == p.arrival_time):
				queue.append(p)
				queue.sort(key=lambda x: (x.tau, x.id))
				if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) arrived; added to ready queue {}".format(time, p.id, p.tau, parseQueue(queue)), flush=True)
		
		if (blocked):
			for blocked_p in blocked:
				if (time == blocked_p[1]): 	# blocked process finished blocking?
					p = blocked_p[0]
					queue.append(p)
					queue.sort(key=lambda x: (x.tau, x.id))
					if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed I/O; added to ready queue {}".format(time, p.id, p.tau, parseQueue(queue)), flush=True)
					blocked.remove(blocked_p)

		if (using_CPU):
			if (time-t_cs//2 == current_burst[1]):	# current CPU burst completed?
				ctime = time -t_cs//2
				using_CPU = False
				p = current_burst[0]
				current_burst = [-1, -1, -1]	# empty current burst info

				old_tau = p.tau
				p.current_burst_no += 1

				if (p.current_burst_no == len(p.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(ctime, p.id, parseQueue(queue)), flush=True)
					#time += t_cs // 2
					processes.remove(p)
					continue
				
				bursts_left = len(p.CPU_burst_times) - p.current_burst_no
				if (bursts_left == 1):
					if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; 1 burst to go {}".format(ctime, p.id, p.tau, parseQueue(queue)), flush=True)
				else:
					if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; {:d} bursts to go {}".format(ctime, p.id, p.tau, bursts_left, parseQueue(queue)), flush=True)

				b = [-1, -1]		# [Process object, when will it be unblocked?]
				b[0] = p
				b[1] = time + p.IO_burst_times[p.current_burst_no - 1]
				blocked.append(b)

				p.tau = math.ceil(p.CPU_burst_times[p.current_burst_no-1] * alpha + (1 - alpha) * p.tau)
				if (time < 10000): print("time {:d}ms: Recalculated tau for process {}: old tau {:d}ms ==> new tau {:d}ms {}".format(ctime, p.id, old_tau, p.tau, parseQueue(queue)), flush=True)

				if (time < 10000): print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}".format(ctime, p.id, b[1], parseQueue(queue)), flush=True)
				
				# time += t_cs // 2
				# continue

		if (not using_CPU and queue):
			using_CPU = True
			p = queue.pop(0)
			p_burst_time = p.CPU_burst_times[p.current_burst_no]
			current_burst[0] = p
			#time += t_cs // 2
			current_burst[1] = time + t_cs // 2 + p_burst_time
			current_burst[2] = time + t_cs // 2
			#if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}".format(time, p.id, p.tau, p_burst_time, parseQueue(queue)))
			#continue

		time += 1

	print("time {:d}ms: Simulator ended for SJF [Q empty]".format(time), flush=True)
			
def SRT(t_cs, alpha, t_slice):
	defaultAllProcesses()

	time = 0
	processes = flattenProcessTable()
	queue = []

	using_CPU = False
	current_burst = [-1, -1, -1, -1]		# [Process object, when will it finish?, when will it start?, preemption?]

	preempted = []
	blocked = []

	print("time 0ms: Simulator started for SRT [Q empty]")

	while (processes):

		if (time == current_burst[2]):
			p = current_burst[0]
			
			if (current_burst[3] == 0):		# not preempted
				p_burst_time = current_burst[1] - current_burst[2]
				if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}"
						.format(time, p.id, p.tau, p_burst_time, parseQueue(queue)), flush=True)
			else: 			# preempted
				remaining = current_burst[1] - time
				p_burst_time = p.CPU_burst_times[p.current_burst_no]
				if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for remaining {:d}ms of {:d}ms burst {}"
		   				.format(time, p.id, p.tau, remaining, p_burst_time, parseQueue(queue)), flush=True)

		for p in processes:
			if (time == p.arrival_time):
				queue.append(p)
				queue.sort(key=lambda x: (x.tau_remaining, x.id[0], x.id[1]))

				if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) arrived; added to ready queue {}".format(time, p.id, p.tau, parseQueue(queue)), flush=True)

		if (using_CPU):
			if (time-t_cs//2== current_burst[1]):	# current CPU burst completed?
				ctime = time - t_cs//2
				using_CPU = False
				p = current_burst[0]
				current_burst = [-1, -1, -1, -1]	# empty current burst info

				old_tau = p.tau
				p.current_burst_no += 1

				if (p.current_burst_no == len(p.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(ctime, p.id, parseQueue(queue)), flush=True)
					#time += t_cs // 2
					processes.remove(p)
					continue
				
				bursts_left = len(p.CPU_burst_times) - p.current_burst_no
				if (bursts_left == 1):
					if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; 1 burst to go {}".format(ctime, p.id, p.tau, parseQueue(queue)), flush=True)
				else:
					if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed a CPU burst; {:d} bursts to go {}".format(ctime, p.id, p.tau, bursts_left, parseQueue(queue)), flush=True)

				p.tau = math.ceil(p.CPU_burst_times[p.current_burst_no-1] * alpha + (1 - alpha) * p.tau)
				p.tau_remaining = p.tau
				if (time < 10000): print("time {:d}ms: Recalculated tau for process {}: old tau {:d}ms ==> new tau {:d}ms {}".format(ctime, p.id, old_tau, p.tau, parseQueue(queue)), flush=True)

				b = [-1, -1]		# [Process object, when will it be unblocked?]
				b[0] = p
				b[1] = time + p.IO_burst_times[p.current_burst_no - 1]
				blocked.append(b)

				if (time < 10000): print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}".format(ctime, p.id, b[1], parseQueue(queue)), flush=True)
				
				# time += t_cs // 2

				# continue

		
		if (blocked):
			for blocked_p in blocked:
				if (time == blocked_p[1]): 			# blocked process finished blocking?
					p = blocked_p[0]
					c = current_burst[0]
					queue.append(p)
					queue.sort(key=lambda x: (x.tau_remaining, x.id[0], x.id[1]))

					if (using_CPU and p.tau < c.tau_remaining):			# preemption
						if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed I/O; preempting {} (predicted remaining time {:d}ms) {}"
								.format(time, p.id, p.tau, c.id, c.tau_remaining, parseQueue(queue)), flush=True)
					
						pre = [-1, -1]						# [Process object, time remaining]
						pre[0] = c
						pre[1] = current_burst[1] - time
						preempted.append(pre)
						queue.insert(0, c)
						queue.remove(p)
						queue.sort(key=lambda x: (x.tau_remaining, x.id[0], x.id[1]))
						current_burst[0] = p
						p_burst_time = p.CPU_burst_times[p.current_burst_no]
						#time += t_cs
						current_burst[1] = time + t_cs + p_burst_time
						current_burst[2] = time + t_cs
						current_burst[3] = 0

						#print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}".format(time, p.id, p.tau, p_burst_time, parseQueue(queue)))

					else:
						
						if (time < 10000): print("time {:d}ms: Process {} (tau {:d}ms) completed I/O; added to ready queue {}".format(time, p.id, p.tau, parseQueue(queue)), flush=True)
					
					blocked.remove(blocked_p)

		if (not using_CPU and queue):
			using_CPU = True
			resume_preempted = False
			p = queue.pop(0)
			for pre in preempted:
				if (pre[0] == p): 
					resume_preempted = True	
					p_burst_time = p.CPU_burst_times[p.current_burst_no]
					current_burst[0] = p
					#time += t_cs // 2
					current_burst[1] = time + t_cs // 2 + pre[1]
					current_burst[2] = time + t_cs // 2
					current_burst[3] = 1

					# print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for remaining {:d}ms of {:d}ms burst {}"
		   			# 		.format(time, p.id, p.tau, pre[1], p_burst_time, parseQueue(queue)))
					
					preempted.remove(pre)
					break

			if (resume_preempted): continue
				
			else:	# ok to switch in process to CPU
				p_burst_time = p.CPU_burst_times[p.current_burst_no]
				current_burst[0] = p
				#time += t_cs // 2
				current_burst[1] = time + t_cs // 2 + p_burst_time
				current_burst[2] = time + t_cs // 2
				current_burst[3] = 0
				#print("time {:d}ms: Process {} (tau {:d}ms) started using the CPU for {:d}ms burst {}".format(time, p.id, p.tau, p_burst_time, parseQueue(queue)))
				#continue

		if (time >= current_burst[2] and current_burst[0] != -1):	# track remaining tau of current burst
			current_burst[0].tau_remaining -= 1

		time += 1

	print("time {:d}ms: Simulator ended for SRT [Q empty]".format(time), flush=True)

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
	lmbda_c = lmbda
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
		print("<<< -- process set (n={:d}) with 1 CPU-bound process".format(n), flush=True)
	else:
		print("<<< -- process set (n={:d}) with {:d} CPU-bound processes".format(n, n_CPU), flush=True)
	
	print("<<< -- seed={:d}; lambda={:.6f}; bound={:d}".format(seed, lmbda, bound), flush=True)

	for i in range(n):
		if (i >= n_CPU):
			type = 1
		
		letter = 65 + i // 10
		number = i % 10

		row = int(letter - 65)
		col = int(number)

		processes_table[row][col].tau = math.ceil(1 / lmbda_c)
		processes_table[row][col].tau_remaining = math.ceil(1 / lmbda_c)

		processes_table[row][col].id = str(chr(letter)) + str(number)


		arrival_time = math.floor(next_exp(lmbda, bound))
		processes_table[row][col].arrival_time = arrival_time

		uniform_distribution = drand48()
		CPU_bursts = math.ceil(uniform_distribution * 32)

		if (type == 0):
			if (CPU_bursts == 1):
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst".format(letter, number, arrival_time), flush=True)
			else:
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU bursts".format(letter, number, arrival_time, CPU_bursts), flush=True)
		else:
			if (CPU_bursts == 1):
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst".format(letter, number, arrival_time), flush=True)
			else:
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU bursts".format(letter, number, arrival_time, CPU_bursts), flush=True)

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

	print(flush=True)
	print("<<< PROJECT PART II", flush=True)
	print("<<< -- t_cs={:d}ms; alpha={:.2f}; t_slice={:d}ms".format(t_cs, alpha, t_slice), flush=True)

	FCFS(t_cs, alpha, t_slice)
	print(flush=True)
	SJF(t_cs, alpha, t_slice)
	print(flush=True)
	SRT(t_cs, alpha, t_slice)

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
