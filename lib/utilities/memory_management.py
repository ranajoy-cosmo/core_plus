#! /usr/bin/python2

import os
import time
import random
import psutil


def pause(min_delay, max_delay):
    chosen_delay = random.uniform(min_delay, max_delay)
    time.sleep(chosen_delay)


def wait_for_free_memory(limit=6e9, delay=0.25, max_delay=None, start_with_pause=False, verbose=False):
    if max_delay is None:
        max_delay = delay
    if start_with_pause:
        pause(delay, max_delay)
    while True:
        usage = psutil.virtual_memory()
        if usage.available > limit:
            print "memory pass"
            break
        if verbose:
            os.system("free -m")
            print "waiting for %i bytes to be free" % limit
            print "    current", usage
        pause(delay, max_delay)


def check_memory(text="", machine=None, num_proc=1):
    """
    Prints the memory usage at this moment.
    text will be printed to the added message.  Useful to know where in the code this was printed.
    machine is the name of the comupter.  If defined in mem_cap_MB, the % of memory used will be printed as well.
    """
    #computer memory capacities
    mem_cap_MB = {"edison": 2500.0}

    #memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)

    #print summary to screen
    log_text = "Memory usage {0:s}: {1:.1f} MB".format(text, mem)
    if not machine is None and machine in mem_cap_MB.keys():
        log_text += ", {0:.1f}%".format(100.0*mem/mem_cap_MB[machine]/num_proc)
    print(log_text)
    return


if __name__ == "__main__":
    wait_for_free_memory(limit=1e12, delay=0.2, verbose=True)
