
% Uniform mutation technique (binary codification)

function M = mutation(S,pm)

[n,N] = size(S);
M = S;

for i = 1:N    
    % Identify the mutate positions of the solution i
    r = find(rand(n,1) <= pm);
    
    % Apply the mutation
    M(r,i) = not(M(r,i));    
end
      

