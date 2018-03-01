
function Zfobj = evaluate(x,params,data)

N = size(x,2);
nb = params.nb;

for i = 1:N
    
    % Convert each diameter value from binary to decimal
    R = base2todec(x(:,i),nb);
    
    % Calculate the network flow for the solution i
    [~,pi,~,cost] = networkflow(R,params,data);
    
    % Determine the number of nodes that violates the pressure bounds
    w = length(pi(pi<params.pmin)) + length(pi(pi>params.pmax));
    
    % State the fitness of the solution i
    Zfobj(i) = cost + w*(params.C(end)-params.C(1))*sum(params.L);
end

