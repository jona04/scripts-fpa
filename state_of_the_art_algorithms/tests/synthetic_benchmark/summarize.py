'''
Read table tests and summarize results on screen or as a latex file.
'''
import sys
import numpy as np
from run_tests import tests, sizes

# Name of the solvers in the order they were run by synthetic_benchmark.
solvers = ['FPA2','FPA','Newton','Secant','Variable fixing', 'Median search']
# solvers = ['FPA2','FPA','Newton','Secant']

# Maps the columns to different statistics for each solver.
its, time, dist, res = {}, {}, {}, {}
for s in range(len(solvers)):
    its[solvers[s]] = 1 + 4*s
    time[solvers[s]] = 2 + 4*s
    dist[solvers[s]] = 3 + 4*s
    res[solvers[s]] = 4 + 4*s

def verify(test):
    '''Verify if all the solvers performed well in the test.

    It verifies if the relative residual is small enough (1.0e-8) and
    optimal value when the problem is viewed as a D-projection is
    approximately equal for all solvers.
    '''
    EPS = 1.0e-8
    test = np.array(test, dtype=int)
    if np.any(test[list(its.values())] < 0):
        return False
    if np.any(test[list(res.values())] > EPS):
        return False
    min_dist = min([test[dist[i]] for i in solvers])
    if np.any(
        (test[list(dist.values())] - min_dist)/max(1.0, min_dist) > np.sqrt(EPS)):
        return False
    return True
          
def screen_report(test, data):
    '''Print in stdout a good representation of the test data as a table.
    '''
    
    # Compute statistics.
    print (test)
    print ('-' * len(test))
    print ()
    for solver in solvers:
        print ('-' * len(solver))
        for size in sizes:
            size_data = data[data[:,0] == size,:]
            its_data = size_data[:,its[solver]]
            time_data = size_data[:,time[solver]]
            print ('%7d \t %.1f \t %d \t %d \t %.1f \t %.1f \t %.1f' % (
                int(size), its_data.mean(), its_data.max(), its_data.min(), 
                1*time_data.mean(), 1*time_data.max(), 
                1*time_data.min() ))
        print ()
    print ()

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
\begin{tabular}{|r||ccc|rrr|r||rrr|rrr|r|}
\hline
\!\!Dimension\!\!\! & 
\multicolumn{3}{c|}{Iterations} &
\multicolumn{3}{c|}{Time (msec)} & 
\multicolumn{1}{c||}{Error} & 
\multicolumn{3}{c|}{Iterations} &
\multicolumn{3}{c|}{Time (msec)} &
\multicolumn{1}{c|}{Error} \\
\multicolumn{1}{|c||}{$n$}&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c|}{\!\!\!min\!\!\!}
&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c|}{\!\!\!min\!\!\!}
&\multicolumn{1}{c||}{---}
&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c|}{\!\!\!min\!\!\!}
&\multicolumn{1}{c}{\!\!\!avg\!\!\!}&\multicolumn{1}{c}{\!\!\!\!max\!\!\!}&\multicolumn{1}{c|}{\!\!\!min\!\!\!}
&\multicolumn{1}{c|}{---} \\ \hline
'''
TABLEEND=r'''
\end{tabular}
\vspace{-2ex}
\caption{%s test}
\end{scriptsize}
\end{center}
\end{table}
'''
TABLELINE = r'&\!\!\! %.1f \!\!\!&\!\!\! %d  \!\!\!&\!\!\! %d  \!\!\!&\!\!\! %.3f  \!\!\!&\!\!\! %.3f  \!\!\!&\!\!\! %.3f & %d'
ENDLINE = r'\\ \hline'
FPAAFPA = r'& \multicolumn{7}{c||}{\bf \hspace{-8ex}FPA2} & \multicolumn{7}{c|}{\bf \hspace{-8ex}FPA} \\ \hline'
NEWTONSECANT = r'& \multicolumn{7}{c||}{\bf Newton} & \multicolumn{7}{c|}{\bf Secant} \\ \hline'
FIXSEARCH = r'& \multicolumn{7}{c||}{\bf Variable fixing} & \multicolumn{7}{c|}{\bf Median search} \\ \hline'

def print_latex_table(size, solvers, data):
    '''Print a latex representation of the test data as a table. 
    '''
    print (r'%7d \!\!\!' % size)
    size_data = data[data[:,0] == size,:]
    for solver in solvers:
        its_data = size_data[:,its[solver]]
        time_data = size_data[:,time[solver]]
        error_iterations = np.where(its_data > 99)
        its_data = np.delete(its_data, error_iterations)
        time_data = np.delete(time_data, error_iterations)
        print (TABLELINE % (
            its_data.mean(), its_data.max(), its_data.min(), 
            1*time_data.mean(), 1*time_data.max(), 1*time_data.min(),len(error_iterations[0])
            ))
    print (ENDLINE)
    
def latex_report(test, data):
    '''Print the full report of the test data in latex format.
    '''
    # Compute statistics.

    print (FPAAFPA)
    for size in sizes:
        print_latex_table(size, ['FPA2', 'FPA'], data)

    print (NEWTONSECANT)
    for size in sizes:
        print_latex_table(size, ['Newton', 'Secant'], data)

    print (FIXSEARCH)
    for size in sizes:
        print_latex_table(size, ['Variable fixing', 'Median search'], data)

# Main program
if __name__ == '__main__':
    if sys.argv[1] == 'latex': print (LATEXSTART)
    
    for test in tests:
        data = np.loadtxt(test)

        # Verify if all solvers did really solve the problem.
        line_num = 1
        for line in data:
            if not verify(line):
                print ('Possible failure at line ', line_num,)
                print ('of', test)
            line_num += 1

        # Print results in the right format.
        if sys.argv[1] == 'latex':
            print (TABLESTART)
            latex_report(test, data)
            print (TABLEEND % test.replace('_', ' '))
        else:
            screen_report(test, data)

    if sys.argv[1] == 'latex': print (LATEXEND)
