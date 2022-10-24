'''
Script that runs all the synthetic tests for the paper.
'''
from subprocess import call, check_output
import os
import scipy as sp

# Set of tests types and sizes
results_dir = 'test_results/'
tests = ['uncorrelated', "weakly_correlated", "correlated", "flow"]
sizes = [50000, 100000, 500000, 1000000, 1500000, 2000000]
# sizes = [50000, 100000]

if __name__ == '__main__':
    # Rebuild the benchmark.
    call(['make', 'clean'])
    call(['make'])

    # Run tests.
    for t in tests:
        # Delete file with old test results.
        try:
            os.remove(t)
        except:
            pass

        # Run the test for each size and append the results.
        for size in sizes:
            call(['./synthetic', t, str(size)])
            data = sp.loadtxt('times.dat')
            n = data.shape[0]
            data = sp.hstack((size*sp.ones((n, 1)), data))
            if size != sizes[0]:
                old_data = sp.loadtxt(results_dir + t)
                data = sp.vstack((old_data, data))
            sp.savetxt(results_dir + t, data)
    # Clean up times.dat file.
    os.remove('times.dat')

    # Generate report in LaTeX format.
    os.chdir(results_dir)
    latex = check_output(['python', '../summarize.py', 'latex'])
    out = open('synthetic_results.tex', 'w', encoding='utf-8')
    out.write(latex.decode('utf-8'))
    out.close()
    os.chdir('..')
