%GASNET - Define a gas network and find cost-efficient solutions.
%   This class performs a mono-objective, cost-minimizing, search for
%   solutions that represent valid configurations for gas network problems.
%   The decision variables are the diameters of pipes. Parameters needed
%   for the optimization are defined in the properties field. The solution 
%   to the non-linear flow problem for gas networks is implemented inside 
%   the class[2] (i.e., no external solver is necessary). Optimizing a
%   particular network with this algorithm - aiGNet - results in a set 
%   of solutions with diverse properties. Parameters listed as GetAccess =
%   private may be changed to finetune the search.
%
%   GASNET public properties:
%      costs - vector with price of each diameter available on the market
%      demands - demand of each network node
%      diameters - vector with available diameters
%      fluidType - 1 for gas, 0 for water
%      minPressures - minimum acceptable pressure in each network node
%      pipes - connectivity matrix of the network
%      pipeLengths - length of each pipe
%      sourceInds - indices of source nodes
%      sourcePressures - pressures given by source nodes
%
%   GASNET set-access-private properties: used in solving the flow problem
%
%   GASNET get-access-private properties: used in aiGNet's optimization
%
%   GASNET public methods:
%      evalnet - check gasnet object and run aiGNet
%      inputtest - verify consistency of properties to solve optimization
%      networkflow - solve the flow problem for the network
%      removeinf - remove infeasible solutions for the nominal case
%
%   GASNET private methods:
%      calculatecost - calculate price of configurations
%      creepmut - run the mutation operator on solutions
%      clone - create clones of current solutions
%      diversify - include new candidate solutions in the population
%      evaluate - determine feasibility and cost of candidate solutions
%      getnetworkdata - determine relevant data to solve flow problem
%      licols - find linearly independent set of columns
%      mutate - mutate clones of current population
%      aiGNet - immune algorithm to search for cost-efficient solutions
%      rankpipes - determine relative importance of each pipe
%      select - select individuals to survive
%      suppress - eliminate 'redundant' solutions
%
%   Example:
%      load Benchmarks;
%      [configurations, costs] = berlin52a.evalnet;
%      [configurations, costs] = eil51b.evalnet;
%      st70a.plotEnable = true; % lets algorithm plot execution details
%      [configurations, costs] = st70a.evalnet;
%      [configurations, costs] = rd100b.evalnet;
%      load CaseStudy;
%      [configurations, costs] = cs_gasnet.evalnet;
%
%   References:
%      [1] Ramos, Eduardo S., Batista, Lucas S.. Natural Gas Pipeline Network 
%          Expansion Under Load-Evolution Uncertainty Using an Evolutionary
%          Algorithm, 2018
%      [2] Nielsen, Hans B. "Methods for Analyzing Pipe Networks", Journal
%          of Hydraulic Engineering, 1989
%      [3] El-Mahdy, Ahmed and Metwalli. "Computer aided optimization of
%          natural gas pipe networks using genetic algorithm", Applied Soft
%          Computing, 2010
%
%   Authors: Eduardo S. Ramos    - edusantiago.ramos@gmail.com
%            Lucas de S. Batista - lusoba@cpdee.ufmg.br

%#ok<*MCSUP>
%#ok<*AGROW>

classdef gasnet
   properties
      fluidType = 1  % type of fluid (0:water, 1:gas)
      
      costs = []     % (nd x 1) cost for each diameter available($/m)
      diameters = [] % (nd x 1) available diameters (mm)
      
      demands = []      % (ni x 1) demand of each internal node (m3/h)
      minPressures = [] % (ni x 1) minimum pressure allowed in nodes (bar)
      
      sourceInds = []      % (ns x 1) indices of source nodes
      sourcePressures = [] % (ns x 1) pressure in sources (bar)
      
      pipes = []       % (np x 2) matrix with connections in the network
      pipeLengths = [] % (np x 1) pipe lengths (in m)
   end
   
   properties (SetAccess = private)
      Ai = [] % (np x ni) connectivity of internal nodes [1]
      As = [] % (np x ns) connectivity of source nodes [1]
      pipeWeights = [] % relative importance of each pipe in the topology
      netData = []     % data relevant for solving the flow problem
      evalParams = []  % parameters for evaluation 
   end
   
   properties (GetAccess = private)
      maxnfe = 10000          % max no. of function evalutations
      initPop = 8             % initial population
      nClones = 4             % no. clones for each cell (must be even)
      epsilon = 1e-3          % diversity criterion
      suppressionLim = 1      % suppression threshold
      newCellsPercentage = 1  % rate of randomly generated new cells
      minMutStrength = 0.02   % minimum mutation strength
      maxMutStrength = 0.10   % maximum mutation strength
      plotEnable = false      % plot execution parameters
   end
   
   methods 
      function g = set.maxnfe(g,value)
         if ~(isnumeric(value) && numel(value)==1 && value>0)
            error('maxnfe: must be single number');
         end
         if round(value)~=value
            warning('maxnfe: non-integer number provided. Rounding up');
            value = ceil(value);
         end
         g.maxnfe = value;
      end
      function g = set.initPop(g,value)
         if ~(isnumeric(value) && numel(value)==1 && value>0)
            error('initPop: must be single number');
         end
         if round(value)~=value
            warning('\ninitPop: non-integer number provided. Rounding up');
            value = ceil(value);
         end
         if mod(value,2)
            error('initPop: must be even number');
         end
         g.initPop = value;
      end
      function g = set.fluidType(g,value)
         if ~(isnumeric(value) && numel(value)==1 && value>=0)
            error('fluidType: must be single number');
         end
         if value~=0 && value~=1
            error('fluidType: must be 0 (water) or 1 (gas)');
         end
         g.fluidType = value;
      end
      
      function g = set.costs(g,value)
         if isempty(value), g.costs = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('costs: must have one column');
         end
         if ~(isnumeric(value) && all(value>=0))
            error('costs: elements must be numeric greater than 0');
         end
         if (size(g.diameters,1) > 0 && size(g.diameters,1)~=size(value,1)) 
            warning(['\ncosts: number of rows is inconsistent '...
                     'with number of diameters defined (%d).'], ...
                      size(g.diameters,1));
         end
         g.costs = value;
      end
      function g = set.diameters(g,value)
         if isempty(value), g.diameters = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('diameters: must have one column');
         end
         if ~(isnumeric(value) && all(value>=0))
            error('diameters: elements must be numeric greater than 0');
         end
         if (size(g.costs,1) > 0 && size(g.costs,1) ~= size(value,1)) 
            warning(['\ndiameters: number of rows is inconsistent '...
                     'with number of costs defined (%d).'], ...
                      size(g.costs,1));
         end
         g.diameters = value;
      end
      
      function g = set.demands(g,value)
         if isempty(value), g.demands = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('demands: must have one column');
         end
         if ~(isnumeric(value) && all(value>=0))
            error('demands: elements must be numeric greater than 0');
         end
         if (size(g.minPressures,1) > 0 && ...
               size(g.minPressures,1) ~= size(value,1)) 
            warning(['\ndemands: number of rows is inconsistent '...
                     'with number of minPressures defined (%d).'], ...
                      size(g.minPressures,1));
         end
         g.demands = value;
      end
      function g = set.minPressures(g,value)
         if isempty(value), g.minPressures = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('minPressures: must have one column');
         end
         if ~(isnumeric(value) && all(value>=0))
            error('minPressures: elements must be numeric greater than 0');
         end
         if (size(g.demands,1) > 0 && size(g.demands,1) ~= size(value,1)) 
            warning(['\nminPressures: number of rows is inconsistent '...
                     'with number of demands defined (%d).'], ...
                      size(g.demands,1));
         end
         g.minPressures = value;
      end
      
      function g = set.sourceInds(g,value)
         if isempty(value), g.sourceInds = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('sourceInds: must have one column');
         end
         if ~(isnumeric(value) && all(value>=0) && ...
               all(round(value)==value))
            error('sourceInds: elements must be numeric greater than 0');
         end
         if (size(g.sourcePressures,1) > 0 && ...
               size(g.sourcePressures,1) ~= size(value,1)) 
            warning(['\nsourceInds: number of rows is inconsistent '...
                     'with number of sourcePressures defined (%d).'], ...
                      size(g.sourcePressures,1));
         end
         g.sourceInds = value;
      end
      function g = set.sourcePressures(g,value)
         if isempty(value), g.sourcePressures = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('sourcePressures: must have one column');
         end
         if ~(isnumeric(value) && all(value>=0))
            error('sourcePressures: must be numeric greater than 0');
         end
         if (size(g.sourceInds,1) > 0 && ...
               size(g.sourceInds,1) ~= size(value,1)) 
            warning(['\nsourcePressures: number of rows is inconsistent '...
                     'with number of sourceInds defined (%d).'], ...
                      size(g.sourceInds,1));
         end
         g.sourcePressures = value;
      end
      
      function g = set.pipes(g,value)
         if isempty(value), g.pipes = value; return; end
         if size(value,2)~=2
            error('pipes: must have two columns');
         end
         if ~(isnumeric(value) && all(all(value>=0)) && ...
               all(all(round(value)==value)))
            error('pipes: elements must be numeric greater than 0');
         end
         g.pipes = value;
      end
      function g = set.pipeLengths(g,value)
         if isempty(value), g.pipeLengths = value; return; end
         if size(value,1)==1, value = value'; end
         if size(value,2)~=1
            error('pipeLengths: must have two columns');
         end
         if ~(isnumeric(value) && all(all(value>=0)))
            error('pipeLengths: elements must be numeric greater than 0');
         end
         g.pipeLengths = value;
      end
      
      function g = inputtest(g)
         %INPUTTEST Verify consistency of properties to solve optimization
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %
         %   Output arguments:
         %      g: 'gasnet' object
         
         n = size(g.costs,1);
         if size(g.diameters,1)~=n
            error('Inconsistent variables: costs and diameters');
         end
         n = size(g.demands,1);
         if size(g.minPressures,1)~=n
            error('Inconsistent variables: demands and minPressures');
         end
         n = size(g.sourceInds,1);
         if size(g.sourcePressures,1)~=n
            error('Inconsistent variables:sourceInds and sourcePressures');
         end
         n = size(g.pipes,1);
         if size(g.pipeLengths,1)~=n
            error('Inconsistent variables: pipes and pipeLengths');
         end
         
         % Cardinalities
         ni = size(g.demands,1);  % number of internal nodes
         ns = size(g.sourceInds,1);  % number of source nodes
         np = size(g.pipes,1);  % number of pipes

         % Connectivity matrices
         g.Ai = sparse([(1:np)' (1:np)'], ...
                    [min(g.pipes,[],2) max(g.pipes,[],2)], ...
                    [ones(np,1) repmat(-1,np,1)], ...
                    np, ni + ns);
         g.As = g.Ai(:,g.sourceInds);
         g.Ai(:,g.sourceInds) = [];
         
         % Network data
         g.netData = g.getnetworkdata();
         
         % Parameters of evaluation
         g.evalParams.ni = ni;  % number of internal nodes
         g.evalParams.np = np;  % number of pipes
         g.evalParams.fluidType = g.fluidType;
         g.evalParams.D = g.diameters;
         g.evalParams.C = g.costs;
         g.evalParams.Q = g.demands;
         g.evalParams.L = g.pipeLengths;
         g.evalParams.sourcePressures = g.sourcePressures;
         g.evalParams.Ai = g.Ai;
         g.evalParams.As = g.As;
         g.evalParams.pmin = g.minPressures;
      end
      
      function [pop,popFits] = evalnet(g)
         %EVALNET Check object's parameters and run aiGNet to optimize
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %
         %   Output arguments:
         %      pop: matrix with final solutions
         %      popFits: fitness - cost - of each solution
         
         g = g.inputtest;
         [pop,popFits] = g.aiGNet;
      end
      
      function [pop,popFits] = removeinf(g,pop,popFits,params,networkData)
         %REMOVEINF Remove infeasible solutions for the nominal demand case
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      pop: a set of candidate solutions
         %      popFits: fitness of candidate solutions
         %      params: parameters to solve flow problem
         %      networkData: data containing properties of the network
         %                   SEE getnetworkdata function
         %
         %   Output arguments:
         %      pop: a set of candidate solutions
         %      popFits: fitness of candidate solutions
         
         if nargin==3,
            params = g.evalParams;
            networkData = g.netData;
         end
         
         % Variables
         mp = g.minPressures; % minimum pressures allowed for each node
         nSolutions = size(pop,2);
         
         toDelete = [];
         for i = 1:nSolutions
             % Calculate the network flow for solution i
             [~,nodePressures,~] = ...
                gasnet.networkflow(pop(:,i),params,networkData);

             % Determine number of nodes that violate pressure bounds
             % [3] Eq. (21)
             nInvPressures = sum(nodePressures < mp);
             
             % Select infeasible for deletion
             if nInvPressures>0, toDelete = [toDelete i]; end
         end
         pop(:,toDelete) = [];
         popFits(toDelete) = [];
      end
      
   end         
         
   methods (Access = private)
      
      function [pop,popFits] = aiGNet(g)
         %OPT_AINET Immune algorithm to search for cost-efficient solutions
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %
         %   Output arguments:
         %      pop: a set of candidate solutions
         %      popFits: fitness of candidate solutions
         
         % Cardinality
         np = size(g.pipes,1); % number of pipes
         
         % Execution variables
         popSize = g.initPop; % initial number of cells
         memCellsSize = [];   % memory cells size
         gen = 1;             % generation counter
         
         % Params to evaluate network
         params = g.evalParams;
         
         % Determine network's parameters
         networkData = g.netData; 

         % Determine importance level of each pipe
         pipeRanks = g.rankpipes();
         
         % Calculate pipeWeights
         maxrank = max(pipeRanks) + 1;
         minrank = min(pipeRanks);
         g.pipeWeights = (maxrank-pipeRanks)./(maxrank-minrank);
         g.pipeWeights = g.pipeWeights(:);
                  
         % estEval,'r-') Generate initial population
         lb = 1;
         ub = numel(g.diameters);
         pop =(lb + (ub-lb).*rand(np,popSize));
         pop = [floor(pop(:,1:popSize/2)) ceil(pop(:,popSize/2+1:end))];
         
         % Evaluation of the candidate solutions
         popFits  = g.evaluate(pop,params,networkData);
         
         % Store the best and mean fitness values found
         bestEval(gen) = min(popFits); 
         meanEval(gen) = mean(popFits);
         gen = gen + 1;
         
         % Update function evaluations counter
         nfe = popSize;
         
         % Stop criterion
         while (nfe < g.maxnfe)
            
            % Fitness normalization into [0,1]
            normalizedFits = (max(popFits)-popFits) / ...
                                 (max(popFits)-min(popFits));

            % Clone generation (nClones for each cell)
            clonedPop = g.clone(pop);

            % Clone mutation (creep mutation) and evaluation
            mutatedPop = g.mutate(clonedPop,normalizedFits);
            mutatedPopFits= g.evaluate(mutatedPop,params,networkData);
            
            % Update number of function evaluations
            nfe = nfe + popSize*g.nClones;

            % Selection
            [pop,popFits]= g.select(pop,popFits,mutatedPop,mutatedPopFits);

            % Store the best and mean fitness values found
            bestEval(gen) = min(popFits); 
            meanEval(gen) = mean(popFits);
             
            % While mean fitness of population improves > epsilon
            if abs(meanEval(gen-1) - meanEval(gen)) > g.epsilon
               % Update the generation counter 
               gen = gen + 1;
               fprintf('\naiGNet: nfe = %d ; max = %d\n',nfe,g.maxnfe);
               continue;
            end    

            % Suppress similar individuals
            [pop,popFits] = g.suppress(pop,popFits);
            nPop = numel(popFits);
          
            % Update memCellsSize
            memCellsSize = [memCellsSize nPop];

            % Generate diversity (through new individuals)
            newPop = g.diversify(np,popSize,lb,ub);
            newPopFits = g.evaluate(newPop,params,networkData);
            nfe = nfe + numel(newPopFits);

            % Update population
            pop  = [pop newPop];
            popFits = [popFits newPopFits];
            popSize  = nPop + numel(newPopFits);   

            % Update the generation counter 
            gen = gen + 1;
            fprintf('\naiGNet: nfe = %d ; max = %d\n',nfe,g.maxnfe);
         end         
         % Apply the last suppression
         [pop,popFits] = g.suppress(pop,popFits);
         [pop,popFits] = g.removeinf(pop,popFits,params,networkData);
         memCellsSize = [memCellsSize numel(popFits)];
         
         if g.plotEnable
            % Plot the evolution of the method
            figure
            tmax = length(bestEval);
            plot(1:tmax,bestEval,'r-')
            title('aiGNet Algorithm')
            xlabel('No. of Generation')
            ylabel('Cost')

            % Plot the evolution of the method
            figure
            t = 1:length(bestEval);
            h = plot(t,bestEval,'k-o',t,meanEval,'r-o',...
                                             'LineWidth',2,'MarkerSize',2);
            legend(h,'Best fitness','Mean fitness')
            title('aiGNet Algorithm')
            xlabel('No. of Generation')
            ylabel('Cost')

            % Plot memory-cells-size evolution
            figure
            t = 1:length(memCellsSize);
            plot(t,memCellsSize,'k-o','LineWidth',2,'MarkerSize',2);
            title('aiGNet Algorithm')
            xlabel('Affinity function execution')
            ylabel('Memory Cells Size')

            % Plot solutions x monetary costs
            figure
            [~,idx] = sort(popFits(1,:));
            Y = popFits(:,idx);
            plot(Y,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',2)
            xlabel('Solution'), ylabel('Monetary cost (nominal condition)')
         end
                  
      end
         
      function data = getnetworkdata(g)
         %GETNETWORKDATA Determine data used in solving the flow problem
         %               for the network
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %
         %   Output arguments:
         %      data: relevant data (SEE networkflow function,[1],[2])
         
         % Disable the following warning
         warning('off','MATLAB:singularMatrix')

         % Cardinalities
         ni = size(g.demands,1); % number of internal nodes
         np = size(g.pipes,1); % number of pipes

         % Extract a linearly independent set of columns of the node-incidence matrix
         [~,idx,idy] = gasnet.licols(g.Ai'); 
         
         % Determine matrices F and G [2]                              
         I = eye(np);
         P = I(:,[idx idy]);  % nonsingular permutation matriz P 
         T = g.Ai' * P;       % Ai'*P = [F G]
         F = T(:,1:ni);       % F must be a nonsingular matrix (ensured)
         G = T(:,ni+1:np);

         invF  = inv(F);
         condF = norm(F,2) * norm(invF,2);
         if condF > 1e10       % F must be nonsingular
             error('F must be a nonsingular matrix.')    
         end

         data.P = P;
         data.F = invF;
         data.G = G;

         % Compute qc and Y according to equations (20) and (21-a,b) [2]
         % qc is any solution to Ai'*q = Q
         qc = P * [invF*g.demands; zeros(np-ni,1)];  %#ok<MINV>
         
         % Y is an m x (m-1) matrix s.t. Ai'*Y = 0
         Y  = P * [-invF*G; eye(np-ni)];   

         data.qc = qc;
         data.Y  = Y;
      end
         
      function pipeRanks = rankpipes(g)
         %RANKPIPES Determine relative importance of each pipe
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %
         %   Output arguments:
         %      pipeRanks: ranking indices
         
         % The closer a pipe is from a source, the smaller its rank.

         % Cardinalities
         ni = size(g.demands,1);  % number of internal nodes
         ns = size(g.sourceInds,1);  % number of source nodes
         np = size(g.pipes,1);  % number of pipes
      
         firstSourceNode = ni+1;
         lastSourceNode = ni + ns;
         pipeRanks = [];

         for sourceNode = firstSourceNode:lastSourceNode
             A = [g.Ai g.As];
             branchRank = zeros(np,1);
             rank = 1;
             nodesToVisit = sourceNode;
             while ~isempty(nodesToVisit)
                 [pipesConnectedToNodes,~] = find(A(:,nodesToVisit)~=0);
                 branchRank(pipesConnectedToNodes,1) = rank;
                 A(:,nodesToVisit) = 0;
                 nextNodes  = sum(A(pipesConnectedToNodes,:),1);
                 A(pipesConnectedToNodes,:) = 0;
                 nodesToVisit = find(nextNodes(1,:)~=0);
                 rank = rank + 1;
             end
             pipeRanks = [pipeRanks branchRank];
         end
         pipeRanks = min(pipeRanks,[],2);
      end
      
      function costVec = evaluate(g,pop,params,networkData)
         %EVALUATE Evaluate fitness function of the population
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      pop: set of candidate solutions
         %      params: parameters to help solve the flow problem
         %              SEE networkflow function,[1],[2]
         %      networkData: relevant data to solve flow problem
         %                   SEE networkflow function,[1],[2]
         %
         %   Output arguments:
         %      costVec: fitness of each solution
         
         % Variables
         mp = g.minPressures; % minimum pressures allowed for each node
         deltaC = max(g.costs) - min(g.costs); % penalty parameter
         networkLength = sum(g.pipeLengths);
         nSolutions = size(pop,2);
         
         costVec = zeros(1,nSolutions);
         parfor i = 1:nSolutions
             % Calculate the network flow for solution i
             [~,nodePressures,~] = ...
                gasnet.networkflow(pop(:,i),params,networkData);
             
             % Calculate cost of configuration
             cost = gasnet.calculatecost(pop(:,i),params);

             % Determine number of nodes that violate pressure bounds
             % [3] Eq. (21)
             nInvPressures = sum(nodePressures < mp);

             % Determine fitness of solution i
             % [3] Eq. (21)
             costVec(i) = cost + nInvPressures*deltaC*networkLength;
         end
      end
      
      function clonedPop = clone(g,pop)
         %CLONEDPOP Create clones for each candidate solution
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      pop: set of candidate solutions
         %
         %   Output arguments:
         %      clonedPop: population with clones added
         
         [nPipes,popSize] = size(pop);
         clonedPop = repmat(pop,g.nClones,1);
         clonedPop = reshape(clonedPop,nPipes,g.nClones*popSize);
      end
      
      function mutatedPop = mutate(g,clonedPop,normalizedFits)
         %MUTATEDPOP Mutate clones of current candidate solutions
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      clonedPop: population with clones added
         %      normalizedFits: fitness values of the solutions, normalized
         %
         %   Output arguments:
         %      mutatedPop: population with clones mutated
         
         % Mutation strengths
         alphas = (g.maxMutStrength-g.minMutStrength) .* normalizedFits;
         alphas = g.maxMutStrength - alphas;
         
         % Radii
         np = size(g.pipes,1);
         radii = alphas.*np;
         radii = repmat(radii,g.nClones,1);
         radii = reshape(radii,1,size(clonedPop,2));
         
         % Apply mutation
         mutatedPop = g.creepmut(clonedPop,radii);
      end
      
      function mutClones = creepmut(g,clonedPop,radii)  
         %CREEPMUT Apply mutation operator on set of solutions
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      clonedPop: population with clones added
         %      radii: measure of mutation strength for each individual
         %             SEE mutate function
         %
         %   Output arguments:
         %      mutClones: population with clones mutated
         
         nSol = size(clonedPop,2); % number of solutions
         np = size(clonedPop,1);   % number of pipes
         
         % Bounds
         lb = 1;
         ub = numel(g.diameters);
                 
         toMutate = (1:nSol);
         deltaPop = zeros(np,nSol);
         tol = 0.1.*radii;
         nMutations = zeros(1,nSol);
         failedAttempts = zeros(1,nSol);
         bestSup = inf(size(radii)); % best mutation with distance > radii
         supSols = zeros(size(clonedPop)); % solutions from 'bestSup'
         bestInf = -inf(size(radii)); % best mutation with distance <=radii
         
         while ~isempty(toMutate)
            n = numel(toMutate);
            randPipes = randi(np,1,n);
            randPipes = (toMutate-1).*np + randPipes;%adjust matrix indices
            randShifts = sign(rand(1,n)-0.5); % determine mutation (+-1)
            
            % Indices of pipes before and after the mutations
            oldPipes = clonedPop(randPipes);
            newPipes = oldPipes + randShifts;
            
            % Truncate
            randShifts(newPipes < lb) = 0; randShifts(newPipes > ub) = 0;
            newPipes = oldPipes + randShifts;
            
            % Update difference between clones and mutated clones
            deltaPop(randPipes) = deltaPop(randPipes) - randShifts;
            
            % Distances from original individuals to mutated
            distVec = arrayfun(@(x) ...
               sqrt(sum((g.pipeWeights.*deltaPop(:,x)).^2)),toMutate);
            
            % Difference between distances and desired radii
            diff = distVec - radii(toMutate);
            
            % Select mutations that are closest to (and less than) radii
            infInds = diff < 0 & diff > bestInf;
            bestInf(infInds) = diff(infInds);
            clonedPop(randPipes(infInds)) = newPipes(infInds);
            
            % Select mutations that are closest to (and greater than) radii 
            supInds = diff > 0 & diff < bestSup & diff < abs(bestInf);
            bestSup(supInds) = diff(supInds);
            supSols(:,toMutate(supInds)) = clonedPop(:,toMutate(supInds));
            supSols(randPipes(supInds)) = newPipes(supInds);
            
            % All mutations
            inds = supInds | infInds;
            
            % Update parameters
            nMutations(inds) = nMutations(inds) + 1;
            failedAttempts(inds) = 0;
            failedAttempts(~inds) = failedAttempts(~inds) + 1;
            
            % Undo difference generated by last failed mutations
            deltaPop(randPipes(~inds))= deltaPop(randPipes(~inds))...
                                           + randShifts(~inds);
            
            % Mutations that are satisfactorily close to radii
            inTol = abs(diff) <= tol(toMutate);
                        
            % Individuals no longer required to mutate
            rem = find((failedAttempts > 3 & nMutations > 0) | inTol);
            
            % Select mutations > radii that are better than those < radii
            s = bestSup(rem) < abs(bestInf(rem));
            clonedPop(:,toMutate(rem(s))) = supSols(:,toMutate(rem(s)));
            
            % Remove individuals
            toMutate(rem) = []; 
            nMutations(rem) = []; 
            failedAttempts(rem) = [];
            bestInf(rem) = [];
            bestSup(rem) = [];
                                    
         end     
         % Set output
         mutClones = clonedPop;
         
      end
      
      function [S,jS] = select(g,pop,popFits,mutPop,mutPopFits)
         %SELECT Select superior individuals from each group of clones
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      pop: a set of candidate solutions
         %      popFits: fitness of candidate solutions
         %      mutPop: a set of candidate solutions
         %      mutPopFits: fitness of candidate solutions
         %
         %   Output arguments:
         %      S: new population
         %      jS: new population fitness
         
         S = zeros(size(pop));
         
         allFits = [popFits; vec2mat(mutPopFits,g.nClones)'];
         [jS,indMin] = min(allFits,[],1);
         fromPop = indMin==1;
         S(:,fromPop) = pop(:,fromPop);
         fromMut = (0:g.nClones:numel(mutPopFits)-1);
         fromMut = fromMut(~fromPop) + indMin(~fromPop) - 1;
         S(:,~fromPop) = mutPop(:,fromMut);
      end

      function [pop,popFits] = suppress(g,pop,popFits)
         %SUPPRESS Remove 'redundant' solutions
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      pop: a set of candidate solutions
         %      popFits: fitness of candidate solutions
         %
         %   Output arguments:
         %      pop: a set of candidate solutions
         %      popFits: fitness of candidate solutions
         
         distVec = pdist(pop'); distVec(distVec==0) = eps;
         distMat = squareform(distVec);

         [rowInds,colInds] = ...
            find(distMat>0 & distMat<g.suppressionLim);

         index = [];
         for k = 1:length(rowInds)
             if colInds(k) > rowInds(k) % consider only triu(dmat)  
                 % "<" for max. and ">" for minimization
                 if popFits(rowInds(k)) > popFits(colInds(k))    
                     index = [index rowInds(k)];
                 else
                     index = [index colInds(k)];
                 end
             end
         end
         pop(:,index) = [];
         popFits(index)  = [];
      end
      
      function newPop = diversify(g,np,popSize,lb,ub)
         %DIVERSIFY Include new, randomly generated, individuals
         %
         %   Input arguments:
         %      g: 'gasnet' object
         %      np: number of pipes
         %      popSize: current number of individuals
         %      lb: minimum diameter index (1)
         %      ub: maximum diameter index
         %
         %   Output arguments:
         %      newPop: new set of candidate solutions
         
         newPopSize = round(g.newCellsPercentage*popSize);
         if (mod(newPopSize,2) ~= 0), newPopSize = newPopSize + 1; end

         P = (lb + (ub-lb).*rand(np,newPopSize));
         newPop=[floor(P(:,1:newPopSize/2)) ceil(P(:,newPopSize/2+1:end))];
      end
      
   end %methods(private)
   
   methods (Static = true)      
      function [gasFlow,nodePressures,minPressure] = ...
                              networkflow(solution,params,networkData)
         %NETWORKFLOW Solve flow problem for the network
         %
         %   Input arguments:
         %      solution: solution to be evaluated
         %      params: parameters to solve the flow problem
         %      networkData: relevant data to solve flow problem
         %
         %   Output arguments:
         %      gasFlow: flow in the network's pipes
         %      nodePressures: pressure in each node
         %      minPressure: minimum node pressure throughout the network
         %
         %   SEE ALSO references
                  
         % Disable the following warnings
         warning('off','optim:fsolve:NonSquareSystem')
         
         % Diameters of the candidate solution.
         diam = params.D(solution);
         diam = diam(:);

         % % Cost of a candidate solution.
         % c = params.C(index)*params.L;

         ni = params.ni; %length(Q);  % number of internal nodes
         np = params.np; %length(L);  % number of pipes

         qc = networkData.qc;
         Y  = networkData.Y;

         % Parameters employed in the application.
         alpha = 1.8570;                             % 1.8  <= alpha <= 2.0
         E = 0.9;                                    % efficiency factor
         K = 18.43 * (params.L./(diam.^(4.854) * E^2)); % flow resistance
         % [Edit] flow resistance for high pressure

         % Influence of the type of fluid: water (0) or gas (1).
         if params.fluidType == 1
             hs = params.sourcePressures.^2;
         else
             hs = params.sourcePressures;
         end

         % Solve the non-linear system of (m-n) equations in the unknown 
         % "u". It follows the equation (14a) (or (14b)). It is important 
         % to note that since hn+1 = ... = hn+r, Y'*As*hs returns a null 
         % vector. However, this implementation works for both 
         % hn+1 = ... = hn+r or hn+1 ~= ... ~= hn+r.

         % Create a random initial vector u0.
         mQ = max(params.Q);
         lb = -mQ;
         ub = +mQ;
         u0 = (repmat(ub,np-ni,1)-repmat(lb,np-ni,1)).*rand(np-ni,1) + ...
               repmat(lb,np-ni,1);


         % Calculate the vector u by using (14b).
         % The flow resistence is given by the widely used power expression
         %  of equation (5a), in which q is calculated according to (13a).
         options = optimset('Algorithm','trust-region-dogleg',...
                            'Display','off');
         [u,~,~,~] = fsolve(@(u) ...
            Y'*( (diag(K.*(abs(qc+Y*u)).^(alpha-1)))*(qc+Y*u) + ...
            params.As*hs ), u0, options); 
         
         % Once u has been calculated, the flow (in m^3/h) at pipes q 
         % (a m x 1 vector) is given according to the linear system shown 
         % in equation (13a), that satisfies (Ai')*(q) - (Q) = 0.
         gasFlow = qc + Y*u;

         % Some consideration: 
         % Hi = bi*qi, Hi = hji2 - hji1, i = 1,...,m (see equation (3))
         % Hi = -sum_{j=1}^{n+r} aij*hj              (see equation (8))
         % H  = -Ai*h - As*hs
         % Finally, h is found by means of equation (9).

         % Create a random initial vector h0.
         lb = params.pmin.^2;
         ub = (min(params.sourcePressures))^2;
         if numel(params.pmin)==1
            h0 = (repmat(ub,ni,1) - repmat(lb,ni,1)).*rand(ni,1) + ...
                  repmat(lb,ni,1);
         else
            h0 = (repmat(ub,ni,1) - lb) .* rand(ni,1) + lb;
         end

         % Solve the non-linear system of equation (9). The solution 
         % satisfies Y'*H = 0, in which H = -Ai*h - As*hs. 
         % The flow resistence is given by the widely used power expression
         %  of equation (5a), in which q is calculated according to (13a).
         [h,~,~,~] = fsolve(@(h) ...
            ((diag(K.*(abs(qc+Y*u)).^(alpha-1)))*(gasFlow) + ...
            (params.Ai)*(h) + (params.As)*(hs)), h0, options);

         % At last, since we are dealing with a gas fluid (see equation 
         % (4b)), the pressure at internal nodes (in bar) is given in the 
         % (n x 1) vector pi.
         if params.fluidType == 1             
             nodePressures = sign(h).*sqrt(abs(h));
         else
             nodePressures = h;
         end

         minPressure = min(nodePressures);
      end
      
   end %methods(static)
   
   methods (Access = private, Static = true)      
      function cost = calculatecost(solution,params)
         %CALCULATECOST Determine INSTALLATION PRICE of a solution
         %
         %   Input arguments:
         %      solution: a configuration, i.e., vector of pipe indices
         %      params: parameters of the network
         %
         %   Output arguments:
         %      cost: installation price of network
         
         costVec = params.C(:)';
         lengthVec = params.L(:);
         cost = costVec(solution)*lengthVec;
      end
      
      function [Xsub,idx,idy]=licols(X,tol)
         %LICOLS - Reference: https://www.mathworks.com/matlabcentral/
         %           answers/49984-how-to-remove-dependent-rows-in-a-matrix
         %
         % Extract a linearly independent set of columns of a given matrix
         %
         %    [Xsub,idx,idy]=licols(X)
         %
         % in:
         %
         %  X: The given input matrix
         %  tol: A rank estimation tolerance. Default=1e-10
         %
         % out:
         %
         % Xsub: The extracted columns of X
         % idx:  The indices (into X) of the extracted columns
         % idy:  The indices (into X) of the remaining columns

         if ~nnz(X) % X has no non-zeros and hence no independent columns
             Xsub=[]; idx=[];
             return
         end

         if nargin<2, tol=1e-10; end

         [~, R, E] = qr(X,0); 

         if ~isvector(R)
             diagr = abs(diag(R));
         else
             diagr = R(1);   
         end

         % Rank estimation
         r = find(diagr >= tol*diagr(1), 1, 'last'); % rank estimation
         idx  = sort(E(1:r));
         Xsub = X(:,idx);  
         idy  = E(r+1:end);
      end
      
   end %methods(private,static)
end %classdef


