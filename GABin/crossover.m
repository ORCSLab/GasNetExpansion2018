
% Single point crossover technique (binary codification)

function C = crossover(P,pc)

[n,N] = size(P);
C = nan(n,N);

for i = 1:2:N-1
    if(rand <= pc)
        % Define the cut position
        ic = round(1 + ((n-1)-1).*rand); 
                
        % Generate the first offspring
        C(1:ic,i)   = P(1:ic,i);
        C(ic+1:n,i) = P(ic+1:n,i+1);
 
        % Generate the second offspring 
        C(1:ic,i+1)   = P(1:ic,i+1);
        C(ic+1:n,i+1) = P(ic+1:n,i);
    else
        % There is not a recombination of the solutions
        C(:,i)   = P(:,i);
        C(:,i+1) = P(:,i+1); 
    end  
end

