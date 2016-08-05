import sys
from termcolor import colored

def prompt(text, flush=True):
    print text
    if flush:
        sys.stdout.flush()

def prompt_warning(text):
    print colored('WARNING :', 'red'), text

