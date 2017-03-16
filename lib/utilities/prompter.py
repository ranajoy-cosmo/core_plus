import sys
from termcolor import colored

def prompt(text, outstream=sys.stdout, color=None):
    if color:
        text = colored(text, color) 
    outstream.write(text)
    outstream.flush()
