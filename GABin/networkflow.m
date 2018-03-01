
function [q, pi, minp, cost] = networkflow(Pi,params,data)

%%%
%   The implementation of this function follows the paper:
%   "Methods for Analysing Pipe Networks" by
%   Hans Bruun Nielsen
%%%


%%% Known quantities:
% d1,...,dm   --> pipe diameters (in mm)
% l1,...,lm   --> pipe lengths (in m)
% Q1,...,Qn   --> demand flow rate at internal nodes (in m^3/h)
% ps1,...,psr --> pressure at source nodes (in bar)
%%%


% A candidate solution
index = Pi;

% Diameters of a candidate solution.
d = params.D(index)';

% Cost of a candidate solution.
c = params.C(index)*params.L;


%%% Unknow quantities:
% pi1,...,pin --> pressure at internal nodes (in bar)
% q1,...,qm   --> flow at pipes (in m^3/h)
%%$


n = params.ni; %length(Q);  % number of internal nodes
m = params.np; %length(L);  % number of pipes

pi = nan(n,1);  % pressure at internal nodes (in bar)
q  = nan(m,1);  % flow at pipes (in m^3/h)


% Disable the following warnings
warning('off','optim:fsolve:NonSquareSystem')


% Computation of qc and Y according to equations (20) and (21-a,b).
qc = data.qc;
Y  = data.Y;


% Parameters employed in the application.
alpha = 1.8570;                             % 1.8  <= alpha <= 2.0
E = 0.9;                                    % efficiency factor
K = 18.43 * (params.L./(d.^(4.854) * E^2)); % flow resistance


% Influence of the type of fluid: water (0) or gas (1).
typefluid = 1;
if typefluid == 1
    hs = params.ps.^2;
else
    hs = params.ps;
end


% Solve the non-linear system of (m-n) equations in the unknown "u". It
% follows the equation (14a) (or (14b)). It is important to note that since
% hn+1 = ... = hn+r, Y'*As*hs returns a null vector. However, this 
% implementation works for both hn+1 = ... = hn+r or hn+1 ~= ... ~= hn+r.

% Create a random initial vector u0.
mQ = max(params.Q);
lb = -mQ;
ub = +mQ;
u0 = (repmat(ub,m-n,1)-repmat(lb,m-n,1)).*rand(m-n,1) + repmat(lb,m-n,1);

% Calculate the vector u by using (14b).
options = optimset('Algorithm','trust-region-dogleg','Display','off');
[u,fval,exitflag,output] = fsolve(@(u) Y'*( B(u,K,qc,Y,alpha)*(qc+Y*u) + ...
    params.As*hs ), u0, options); 


% Once u has been calculated, the flow (in m^3/h) at pipes q (a m x 1 
% vector) is given according to the linear system shown in equation (13a), 
% that satisfies (Ai')*(q) - (Q) = 0.
q = qc + Y*u;


% Some consideration: 
% Hi = bi*qi, Hi = hji2 - hji1, i = 1,...,m (see equation (3))
% Hi = -sum_{j=1}^{n+r} aij*hj              (see equation (8))
% H  = -Ai*h - As*hs
% Finally, h is found by means of equation (9).

% Create a random initial vector h0.
lb = params.pmin^2;
ub = params.pmax^2;
h0 = (repmat(ub,n,1)-repmat(lb,n,1)).*rand(n,1) + repmat(lb,n,1);

% Solve the non-linear system of equation (9). The solution satisfies 
% Y'*H = 0, in which H = -Ai*h - As*hs. 
[h,fval,exitflag,output] = fsolve(@(h) (B(u,K,qc,Y,alpha))*(q) + ...
    (params.Ai)*(h) + (params.As)*(hs), h0, options);


% At last, since we are dealing with a gas fluid (see equation (4b)), the 
% pressure at internal nodes (in bar) is given in the (n x 1) vector pi.
if typefluid == 1
    pi = sign(h).*sqrt(abs(h));
else
    pi = h;
end

minp = min(pi);
cost = c;

end





% The flow resistence is given by the widely used power expression of
% equation (5a), in which q is calculated according to (13a).
function y = B(u,K,qc,Y,alpha)

y = diag(K.*(abs(qc+Y*u)).^(alpha-1));

end


