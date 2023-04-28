'''
Script that runs all the synthetic tests for the paper.
'''
from subprocess import call, check_output
import os
import scipy as sp
import numpy as np

# Set of tests types and sizes
results_dir = 'test_results/'
tests = ['uncorrelated', "weakly_correlated", "correlated", "flow"]
sizes = [500000, 1000000, 10000000, 50000000]


if __name__ == '__main__':
    # # Rebuild the benchmark.
    # call(['make', 'clean'])
    # call(['make'])

    # # Run tests.
    # for t in tests:
    #     # Delete file with old test results.
    #     try:
    #         os.remove(t)
    #     except:
    #         pass

    #     # Run the test for each size and append the results.
    #     for size in sizes:
    #         call(['./synthetic', t, str(size)])
    #         data = np.loadtxt('times.dat')
    #         n = data.shape[0]
    #         print("teste 1")
    #         try:
    #             print("teste 2")
    #             data = np.hstack((size*np.ones((n, 1)), data))
    #             if size != sizes[0]:
    #                 old_data = np.loadtxt(results_dir + t)
    #                 data = np.vstack((old_data, data))
    #             np.savetxt(results_dir + t, data)
    #         except:
    #             print("teste 3")
    #             pass
    #         finally:
    #             print("teste 4")
    #             if size != sizes[0]:
    #                 old_data = np.loadtxt(results_dir + t)
    #                 data = np.vstack((old_data, data))
    #             np.savetxt(results_dir + t, data)

    # # Clean up times.dat file.
    # os.remove('times.dat')

    # Generate report in LaTeX format.
    os.chdir(results_dir)
    latex = check_output(['python', '../summarize.py', 'latex'])
    out = open('synthetic_results.tex', 'w', encoding='utf-8')
    out.write(latex.decode('utf-8'))
    out.close()
    os.chdir('..')
