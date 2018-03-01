%
% Genetic Algorithm for Natural Gas Pipe Networks Optimization
% Binary Codification
%

function [P,jP] = ga(networkname,N)
%
% networkname: name of the network simulated
% N          : population size
%


%
% Run example
%		 [configurations, costs] = ga('casestudy',4);
%		 [configurations, costs] = ga('berlin52a',4);
%		 [configurations, costs] = ga('eil51b',4);
%		 [configurations, costs] = ga('st70a',4);
%		 [configurations, costs] = ga('rd100b',4);
%


% % Open a pool of Matlab for parallel computation
% if matlabpool('size') == 0
%     matlabpool open
% end


% Method's parameters
g = 1;          % Generation counter
emax = 10000;   % Maximum number of function evaluations
pc = 1.00;      % Crossover probability
pm = 0.05;      % Mutation probability
nb = 3;         % Number of bits per each pipe codification

% Network's parameters
cd 'networks'
load(networkname)
cd ..
params.nb = nb;
data = networkdata(params);


% Generation of the initial population
P = round(rand(nb*params.np,N)); 

% Evaluation of the candidate solutions
jP = evaluate(P,params,data);
nfe = N;

% Store the best fitness found
fmin = min(jP);
fvec(g) = fmin;

% Cost of the solutions found by the engineers A, B and C
a(g) = 300276200;
b(g) = 324824500;
c(g) = 301744450;
med(g) = (a(g)+b(g)+c(g))./3;


% Stop criterion
while (nfe < emax)
    
    % Binary crossover 
    C = crossover(P,pc); 

    % Binary mutation
    M  = mutation(C,pm);     
    jM = evaluate(M,params,data);
    nfe = nfe + N;
    
    % Roulette wheel selection
    [U,jU] = selection([P M],[jP jM],N);
    
    % Elitism strategy
    [P,jP,fmin] = elitism([P M],[jP jM],U,jU);
    
    % Update the generation counter 
    g = g + 1;
    
    % Store the best fitness found
    fvec(g) = fmin;
    
    % Cost of the solutions found by the engineers A, B and C
    a(g) = 300276200;
    b(g) = 324824500;
    c(g) = 301744450;    
    med(g) = (a(g)+b(g)+c(g))./3;
    
    %g
end
P = base2todec(P,nb);

% Plote the evolution of the method
figure
tmax = length(fvec);
if strcmp(networkname,'casestudy')
    plot(1:tmax,fvec,'r-',1:tmax,a,'b-',1:tmax,b,'k-',1:tmax,med,'c-')
    legend('GA Solution','Designer Min. Cost','Designer Max. Cost','Designer Average Cost')
else
    plot(1:tmax,fvec,'r-')    
end
title('Binary GA')
xlabel('No. of Generation')
ylabel('Cost')

