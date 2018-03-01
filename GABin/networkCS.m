
function params = networkCS()


% Cost (in $/m).
C = [1637 1796 2122 2908 2940 4139 4139 4139];

% Commercial diameters (in mm).
D = [100 150 200 250 300 400 400 400];

% Demand flow rate at internal nodes (in m^3/h).
Q = [11500 11500 11500 11500 0 11500 8750 8750 0 8750 8750 12500]'; 

% Pipe lengths (in m).
L = [6300 8400 7700 9100 1400 6300 12600 3850 4900 1400 4900 4200 8400 ...
    1750 7700 5600 14700 17500 10500 15400 9100]';

pmin = 2.5;  % minimum pressure at any internal node (in bar)
pmax = 17.5; % maximum pressure at any node (in bar)

r  = 2;              % number of reservoir nodes
ps = pmax*ones(r,1); % pressure at source nodes (in bar)

ni = length(Q);  % number of internal nodes
np = length(L);  % number of pipes


% The numbering is such that the pipes and interior nodes furthest away 
% from reservior nodes have the smallest numbers (however, the figure 3 of
% El-Mahdy et al. (2010) is considered in this test).

% The signs of the pipe discharges are defined as follows: Let the ith pipe
% connect nodes number jil and ji2 with jil < ji2. The discharge qi is 
% positive when the flow direction is from node number ji2 to node number 
% jil.

% Node incidence matrix (A):
% aij = +1 if j = ji1
% aij = -1 if j = ji2
% aij =  0 otherwise
% i = 1,...,m; j = 1,...,n, and ji1 < ji2

% Incidence matrix for the test network of El-Mahdy et al. (2010).
A = [
1 0 0 0 0 0 0 0 0 0 0 0 -1 0
1 -1 0 0 0 0 0 0 0 0 0 0 0 0
0 1 -1 0 0 0 0 0 0 0 0 0 0 0
0 0 1 -1 0 0 0 0 0 0 0 0 0 0
0 0 0 1 -1 0 0 0 0 0 0 0 0 0
0 0 0 0 1 -1 0 0 0 0 0 0 0 0
0 0 0 0 0 1 -1 0 0 0 0 0 0 0
0 0 0 0 0 0 1 -1 0 0 0 0 0 0
0 0 0 0 0 0 0 1 -1 0 0 0 0 0
0 0 0 0 0 0 0 0 1 -1 0 0 0 0
0 0 0 0 0 0 0 0 0 1 -1 0 0 0
0 0 0 0 0 0 0 0 0 0 1 -1 0 0
0 0 0 0 0 0 0 0 1 0 0 0 -1 0
0 0 0 0 1 0 0 0 0 0 0 0 0 -1
0 0 0 0 0 0 0 1 0 -1 0 0 0 0
0 0 0 0 0 0 1 0 0 0 -1 0 0 0
0 0 0 0 0 1 0 0 0 0 0 -1 0 0
0 0 0 1 0 0 0 0 0 0 0 -1 0 0
0 1 0 -1 0 0 0 0 0 0 0 0 0 0
0 1 0 0 0 0 0 0 0 0 -1 0 0 0
0 0 0 0 0 0 0 0 0 1 0 0 -1 0
];

Ai = A(:,1:ni);      % interior node incidence matrix (Ai)
As = A(:,ni+1:end);  % source node incidence matrix (As)


params.Ai = Ai;
params.As = As;
params.C = C;
params.D = D;
params.L = L;
params.Q = Q;
params.r = r;
params.ni = ni;
params.np = np;
params.ps = ps;
params.pmin = pmin;
params.pmax = pmax;

