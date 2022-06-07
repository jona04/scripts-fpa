'''
Simple script to run all the tests.

Copyright: Paulo J. S. Silva 2012
'''

import os
from subprocess import call

# Run sythetic tests
os.chdir('synthetic_benchmark')
call(['python', '-u', 'run_tests.py'])
os.chdir('..')

# Run SVM tests
# os.chdir('svm_data/adult')
# call(['python', '-u', 'get_adult.py'])
# os.chdir('../mnist')
# call(['python', '-u', 'get_mnist.py'])
# os.chdir('../../svm')
# call(['python', '-u', 'run_tests.py'])
# os.chdir('..')
