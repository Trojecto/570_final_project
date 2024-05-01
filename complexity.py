import sys
from resource import * import time
import psutil


def process_memory():
	process = psutil.Process() memory_info = process.memory_info()
	memory_consumed = int(memory_info.rss/1024) 
	return memory_consumed

def time_wrapper(): 
	start_time = time.time() 
	call_algorithm()
	end_time = time.time()
	time_taken = (end_time - start_time)*1000 
	return time_taken