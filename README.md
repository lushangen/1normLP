# 1normLP

This is a solver for the 1-norm regression problem, as outlined in the pdf Project(2).
The main solver is OneNormLP3033954135.m, with two helper methods to check degeneracy conditions.
test.m contains two test cases: one is a randomly generated 6x9, and the other is an always degenerate case.
Run test.m to see. 
'''
linprog(b, [], [], A', zeros(size(A,2),1), -ones(size(A,1),1), ones(size(A,1),1));
''' runs Matlab's built in solver to compare results.
