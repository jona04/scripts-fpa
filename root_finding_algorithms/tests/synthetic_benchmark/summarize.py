'''
Read table tests and summarize results on screen or as a latex file.
'''
import sys
import numpy as np
from run_tests import tests, sizes

# Name of the solvers in the order they were run by synthetic_benchmark.
solvers = ['FPA', 'Secant', 'Regula Falsi', 'Bissection']

# Maps the columns to different statistics for each solver.
its, time, dist, res = {}, {}, {}, {}
for s in range(len(solvers)):
    its[solvers[s]] = 1 + 2*s
    time[solvers[s]] = 2 + 2*s
    dist[solvers[s]] = 3 + 2*s
    res[solvers[s]] = 4 + 2*s

def verify(test):
    '''Verify if all the solvers performed well in the test.

    It verifies if the relative residual is small enough (1.0e-8) and
    optimal value when the problem is viewed as a D-projection is
    approximatelyh equal for all solvers.
    '''
    EPS = 1.0e-4
    # if np.any(test[its.values()] < 0):
        # return False
    # if np.any(test[res.values()] > EPS):
    #     return False
    # min_dist = min([test[dist[i]] for i in solvers])
    # if np.any(
        # (test[dist.values()] - min_dist)/max(1.0, min_dist) > np.sqrt(EPS)):
        # return False
    return True
          
def screen_report(test, data):
    '''Print in stdout a good representation of the test data as a table.
    '''
    
    # Compute statistics.
    print test
    print '-' * len(test)
    print
    for solver in solvers:
        print solver
        print '-' * len(solver)
        for size in sizes:
            size_data = data[data[:,0] == size,:]
            its_data = size_data[:,its[solver]]
            time_data = size_data[:,time[solver]]
            print '%7d \t %.1f \t %d \t %d \t %.1f \t %.1f \t %.1f' % (
                int(size), its_data.mean(), its_data.max(), its_data.min(), 
                1000*time_data.mean(), 1000*time_data.max(), 
                1000*time_data.min() )
        print
    print

# Formatting strings for Latex ouput
LATEXSTART=r'''
\documentclass{article}

% Some mathematical macros.
\usepackage{amsfonts}
\usepackage{amssymb,amsopn,amsgen,amsbsy,amsmath}

\begin{document}
'''
LATEXEND=r'''
\end{document}
'''
TABLESTART=r'''
\begin{table}[t]
\begin{center}
\begin{scriptsize}
\begin{tabular}{|r||ccc|rrr||rrr|rrr||}
\hline
\!\!Dimension\!\!\! & 
\multicolumn{3}{c|}{Iterations} &
\multicolumn{3}{c||}{Time (msec)} & 
\multicolumn{3}{c|}{Iterations} &
\multicolumn{3}{c||}{Time (msec)} \\
\multicolumn{1}{|c||}{$n$}&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c|}{\!\!\!min\!\!\!}
&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c||}{\!\!\!min\!\!\!}
&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c|}{\!\!\!min\!\!\!}
&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c||}{\!\!\!min\!\!\!} \\ \hline
'''
TABLEEND=r'''
\end{tabular}
\vspace{-2ex}
\caption{%s test}
\end{scriptsize}
\end{center}
\end{table}
'''
TABLELINE = r'&\!\!\! %.1f \!\!\!&\!\!\! %d  \!\!\!&\!\!\! %d  \!\!\!&\!\!\! %.1f  \!\!\!&\!\!\! %.1f  \!\!\!&\!\!\! %.1f'
ENDLINE = r'\\ \hline'
NEWTONSECANTHEAD = r'& \multicolumn{6}{c||}{\bf \hspace{-8ex}FPA} & \multicolumn{6}{c||}{\bf \hspace{-8ex}Secant} \\ \hline'
FIXSEARCHHEAD = r'& \multicolumn{6}{c||}{\bf Regula Falsi} & \multicolumn{6}{c||}{\bf Bissection} \\ \hline'

def print_latex_table(size, solvers, data):
    '''Print a latex representation of the test data as a table. 
    '''
    print r'%7d \!\!\!' % size
    size_data = data[data[:,0] == size,:]
    for solver in solvers:
        its_data = size_data[:,its[solver]]
        time_data = size_data[:,time[solver]]
        print TABLELINE % (
            its_data.mean(), its_data.max(), its_data.min(), 
            1000*time_data.mean(), 1000*time_data.max(), 1000*time_data.min()
            )
    print ENDLINE
    
def latex_report(test, data):
    '''Print the full report of the test data in latex format.
    '''
    # Compute statistics.
    print NEWTONSECANTHEAD
    for size in sizes:
        print_latex_table(size, ['FPA', 'Secant'], data)

    print FIXSEARCHHEAD
    for size in sizes:
        print_latex_table(size, ['Regula Falsi',
                                 'Bissection'], data)
    print

# Main program
if __name__ == '__main__':
    if sys.argv[1] == 'latex': print LATEXSTART
    
    for test in tests:
        data = np.loadtxt(test)

        # Verify if all solvers did really solve the problem.
        line_num = 1
        for line in data:
            if not verify(line):
                print 'Possible failure at line ', line_num,
                print 'of', test
            line_num += 1

        # Print results in the right format.
        if sys.argv[1] == 'latex':
            print TABLESTART
            latex_report(test, data)
            print TABLEEND % test.replace('_', ' ')
        else:
            screen_report(test, data)

    if sys.argv[1] == 'latex': print LATEXEND
