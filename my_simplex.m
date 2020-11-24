% Name: Prasidh Sriram
% RUID: 178004060
% Assignment: Final Project

function tab = my_simplex(A,b,c,debug)

% Verifies that vector dimensions match
validateInput(A,b,c);

[A,b] = makePhaseOne(A,b,debug);
tab   = simplex(A,b,c,debug,'phase two')
printSolution(tab);
end

% Generate initial tableaux
function [A,b] = makePhaseOne(A,b,debug)
[m,n]              = size(A);
tab                = zeros(m+1,n+m+1);
tab(1:m,1:n)       = A;
tab(end,n+1:end-1) = 1;
tab(1:m,end)       = b(:);
tab(1:m,n+1:n+m)   = eye(m);

% Debug flag determines whether to print each intermediate step
if debug
    fprintf('>>>>Current tableau [phase one]\n');
    disp(tab);
end

for i = 1:m %now make all entries in bottom row zero
    tab(end,:) = tab(end,:)-tab(i,:);
end

tab = simplex(tab(1:m,1:n+m),tab(1:m,end),tab(end,1:n+m),...
                                             debug,'phase one');

A = tab(1:m,1:n);
b = tab(1:m,end);

end

function tab = simplex(A,b,c,debug,phase_name)
[m,n]        = size(A);
tab          = zeros(m+1,n+1);
tab(1:m,1:n) = A;
tab(m+1,1:n) = c(:);
tab(1:m,end) = b(:);

%Boolean determines whether optimal soln has been found
solnNotFound = true;

while solnNotFound
    % Determines whether to print each intermediate step
    if debug
        fprintf('***********************\n');
        fprintf('Current tableau [%s] \n',phase_name);
        disp(tab);
    end
    
    if any(tab(end,1:n)<0)% check if there is negative coeff in bottom row
        [~,J] = min(tab(end,1:n)); % if yes, find the most negative
        % now check if corresponding column is unbounded
        if all(tab(1:m,J)<=0) 
          error('problem unbounded. All entries <= 0 in column %d',J);
          
        % Perform row operations to make all column entries 0
        % Except pivot
        else 
            pivot_row = 0;
            min_found = inf;
            for i = 1:m
                if tab(i,J)>0
                    tmp = tab(i,end)/tab(i,J);
                    if tmp < min_found
                        min_found = tmp;
                        pivot_row = i;
                    end
                end
            end
            if debug
                fprintf('pivot row is %d\n',pivot_row);
            end
            
            %normalize
            tab(pivot_row,:) = tab(pivot_row,:)/tab(pivot_row,J); 
            
            % Perform elementary row ops
            for i=1:m+1
                if i ~= pivot_row
                    tab(i,:)=tab(i,:)-sign(tab(i,J))*...
                                 abs(tab(i,J))*tab(pivot_row,:);
                end
            end
        end
        
        % Elementary row ops done, now print new solution
        if debug 
            fprintf('current basic feasible solution is\n');
            disp(get_current_x());
        end
    else
        solnNotFound=false;       
    end
end

    % Finds current basis vector
    function current_x = get_current_x()
        current_x = zeros(n,1);
        for j=1:n
            if length(find(tab(:,j)==0))==m
                idx= tab(:,j)==1;
                current_x(j)=tab(idx,end);
            end
        end
    end
end

% Makes sure vector/matrix dimensions match
function validateInput(A,b,c)
if ~ismatrix(A)
    error('A must be matrix');
end

if ~isvector(b)
    error('b must be vector');
end
if ~isvector(c)
    error('c must be vector');
end

[m,n]=size(A);
if rank(A) <m
    error('Rank A must be equal to number of rows in A');
end

if length(b) ~= m
    error('b must have same size as number of rows of A');
end
if length(c) ~= n
    error('c must have same size as number of columns of A');
end
end

function printSolution(tab)
[nRow,nCol] = size(tab);

A = tab(1:nRow-1,1:nCol-1);
b = tab(1:nRow-1,nCol);

q = A ~= 0;
q = find(sum(q,1)==1); %find all columns with one non-zero entry;

solnVector = zeros(nCol-1,1); 

for n=1:length(q)
    j = find(A(1:nRow-1,q(n))==1);
    if isempty(j)
       solnVector(q(n)) = 0;
    else
       solnVector(q(n)) = b(j); 
    end;
end;

% Print solution
solnVector

end