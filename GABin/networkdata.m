
function data = networkdata(params)

% Disable the following warning
warning('off','MATLAB:singularMatrix')

n = params.ni; % number of internal nodes
m = params.np; % number of pipes

% Computation of qc and Y according to equations (20) and (21-a,b).
[~,idx,idy] = licols(params.Ai'); % Extract a linearly independent set of 
                                  % columns of a given matrix
I = eye(m);
P = I(:,[idx idy]);   % nonsingular permutation matriz P 
T = params.Ai' * P;   % Ai'*P = [F G]
F = T(:,1:n);         % F must be a nonsingular matrix (ensured)
G = T(:,n+1:m);

invF  = inv(F);
condF = norm(F,2) * norm(invF,2);
if condF > 1e10       % F must be a nonsingular matrix
    error('F must be a nonsingular matrix.')    
end

qc = P * [invF*params.Q; zeros(m-n,1)]; % qc is any solution to Ai'*q = Q
Y  = P * [-invF*G; eye(m-n)];    % Y is a m x (m-1) matrix s.t. Ai'*Y = 0

data.qc = qc;
data.Y  = Y;


