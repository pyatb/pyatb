import os
from pyatb.timer import timer
from pyatb.parallel import COMM, SIZE, RANK

'''
This file is used to initialize the pyatb package, 
including the creation of a workspace folder, and some global variables.
'''

# workspace folder
INPUT_PATH = os.getcwd()
OUTPUT_PATH = os.path.join(INPUT_PATH, 'Out')
RUNNING_LOG = None

if RANK == 0:
    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)

    RUNNING_LOG = os.path.join(OUTPUT_PATH, 'running.log')
else:
    RUNNING_LOG = os.path.join(OUTPUT_PATH, 'running-' + str(RANK) + '.log')

timer = timer(RUNNING_LOG, RANK, SIZE)
timer.program_start()

from pyatb import io
from pyatb.init_tb import init_tb

# if 'OMP_NUM_THREADS' not in os.environ:
#     os.environ['OMP_NUM_THREADS'] = '1'
