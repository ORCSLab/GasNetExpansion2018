
% Elitist strategy (ensure the preservation of the best solution found)

function [U,jU,fmin] = elitism(P,jP,U,jU)

% Identify the best solution and stores it into B
[fmin,jmin] = min(jP);
B  = P(:,jmin);
jB = jP(jmin);

% Check if the best solutions already belongs to U
flag = ismember(B',U','rows');

% If flag is false, the best solution replaces the worst solution of U
if (flag == 0)    
    [~,jmax] = max(jU);
    U(:,jmax) = B;
    jU(jmax)  = jB;
end
