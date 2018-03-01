
% Roulette wheel selection

function [S,jS] = selection(P,jP,N)

fitness = 1./jP;
sumfit  = sum(fitness);
ro      = fitness./sumfit;

for i = 1:N 
	u = rand;
	j = 1;
	s = 0;
    
    while (s < u) 	
        s = s + ro(1,j);
	    j = j + 1;
    end
    
    S(:,i) = P(:,j-1);
    jS(i)  = jP(j-1); 
end


