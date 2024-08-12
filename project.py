import sys
import math

class Process:
	def __init__(self):
		self.id = -1
		self.arrival_time = -1
		self.CPU_burst_times = []
		self.IO_burst_times = []
		self.turnaround_times = []
		self.wait_times = []
		self.current_burst_no = 0

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
COLS = 9
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

def FCFS(t_cs, alpha, t_slice):
	time = 0
	processes = flattenProcessTable()
	queue = []
	
	using_CPU = False
	current_burst = [-1, -1]		# Process object, when will it finish?

	blocked = []
	blocking_on_IO = False

	print("time 0ms: Simulator started for FCFS [Q empty]")
	while (processes):
		for p in processes:
			if (time == p.arrival_time):
				queue.append(p)
				print("time {:d}ms: Process {} arrived; added to ready queue {}".format(time, p.id, parseQueue(queue)))

		if (blocked):
			for blocked_p in blocked:
				if (time == blocked_p[1]): 	# blocked process finished blocking?
					p = blocked_p[0]
					queue.append(p)
					print("time {:d}ms: Process {} completed I/O; added to ready queue {}".format(time, p.id, parseQueue(queue)))
					blocked.remove(blocked_p)
		
		if (using_CPU):
			if (time == current_burst[1]):	# current CPU burst completed?
				using_CPU = False
				p = current_burst[0]
				current_burst = [-1, -1]	# empty current burst info
				p.current_burst_no += 1

				if (p.current_burst_no == len(p.CPU_burst_times)):
					print("time {:d}ms: Process {} terminated {}".format(time, p.id, parseQueue(queue)))
					time += t_cs // 2
					processes.remove(p)
					continue
				
				bursts_left = len(p.CPU_burst_times) - p.current_burst_no
				if (bursts_left == 1):
					print("time {:d}ms: Process {} completed a CPU burst; 1 burst to go {}".format(time, p.id, parseQueue(queue)))
				else:
					print("time {:d}ms: Process {} completed a CPU burst; {:d} bursts to go {}".format(time, p.id, bursts_left, parseQueue(queue)))

				b = [-1, -1]		# Process object, when will it be unblocked?
				b[0] = p
				b[1] = time + t_cs // 2 + p.IO_burst_times[p.current_burst_no - 1]
				blocked.append(b)

				print("time {:d}ms: Process {} switching out of CPU; blocking on I/O until time {:d}ms {}".format(time, p.id, b[1], parseQueue(queue)))

				time += t_cs // 2
				continue
				
		if (not using_CPU and queue):
			using_CPU = True
			p = queue.pop(0)
			p_burst_time = p.CPU_burst_times[p.current_burst_no]
			current_burst[0] = p
			time += t_cs // 2
			current_burst[1] = time + p_burst_time
			print("time {:d}ms: Process {} started using the CPU for {:d}ms burst {}".format(time, p.id, p_burst_time, parseQueue(queue)))
			continue

		time += 1

	print("time {:d}ms: Simulator ended for FCFS [Q empty]".format(time))
			

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
	
	print("<<< PROJECT PART I")
	if (n_CPU == 1):
		print("<<< -- process set (n={:d}) with 1 CPU-bound process".format(n))
	else:
		print("<<< -- process set (n={:d}) with %{:d} CPU-bound processes".format(n, n_CPU))
	
	print("<<< -- seed={:d}; lambda={:.6f}; bound={:d}".format(seed, lmbda, bound))

	for i in range(n):
		if (i >= n_CPU):
			type = 1
		
		letter = 65 + i // 10
		number = i % 10

		row = int(letter - 65)
		col = int(number)

		processes_table[row][col].id = str(chr(letter)) + str(number)

		arrival_time = math.floor(next_exp(lmbda, bound))
		processes_table[row][col].arrival_time = arrival_time

		uniform_distribution = drand48()
		CPU_bursts = math.ceil(uniform_distribution * 32)

		if (type == 0):
			if (CPU_bursts == 1):
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst".format(letter, number, arrival_time))
			else:
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU burst".format(letter, number, arrival_time, CPU_bursts))
		else:
			if (CPU_bursts == 1):
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst".format(letter, number, arrival_time))
			else:
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU burst".format(letter, number, arrival_time, CPU_bursts))

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
