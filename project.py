import sys
import math

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

		arrival_time = math.floor(next_exp(lmbda, bound))

		uniform_distribution = drand48()
		CPU_bursts = math.ceil(uniform_distribution * 32)

		if (type == 0):
			if (CPU_bursts == 1):
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst:".format(letter, number, arrival_time))
			else:
				print("CPU-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU burst:".format(letter, number, arrival_time, CPU_bursts))
		else:
			if (CPU_bursts == 1):
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; 1 CPU burst:".format(letter, number, arrival_time))
			else:
				print("I/O-bound process {:c}{:d}: arrival time {:d}ms; {:d} CPU burst:".format(letter, number, arrival_time, CPU_bursts))

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
			else:
				IO_burst_time *= 8
				IO_total_CPU_burst_time += CPU_burst_time
				IO_num_CPU_burst += 1
				if (IO_burst_time != -8):
					IO_total_IO_burst_time += IO_burst_time
					IO_num_IO_burst += 1
	
	print()
	print("<<< PROJECT PART II")
	print("<<< -- t_cs={:d}ms; alpha={:.2f}; t_slice={:d}ms\n".format(t_cs, alpha, t_slice))

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
