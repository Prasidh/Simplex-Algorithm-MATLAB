% Name: Prasidh Sriram
% RUID: 178004060
% Assignment: Final Project

varInput = 'Enter number of decision variables: ';
numVars = input(varInput);

constrInput = 'Enter number of constraints: ';
numConstraints = input(constrInput);

% Initialize A,b,c with random values on the interval [1 10]
A = randi(10,numConstraints,numVars);
b = randi(10,numConstraints,1);
c = randi(10,numVars,1);

% For manual LP input:
% \/\/\/\/\/\/\/\/

%A=[3,1,1,0;1,2,0,1];
%b=[8,9];
%c=[-4,-5,0,0];

tic
format short;
my_simplex(A,b,c,true)

% Save time elapsed
timeElapsed = toc;

fprintf('Elapsed time is %d seconds\n',timeElapsed);