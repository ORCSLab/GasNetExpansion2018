%ANALYSIS - Multi-objective robustness analysis of a gas network
%   This class performs the multi-objective robustness analysis procedure
%   for gas network described in [1]. It takes as input a gas network from
%   the class 'gasnet' and additional parameters that guide the
%   mono-objective search for cheap solutions/ multi-objective analysis. In
%   particular, suppose 'refYears' is defined as [0,5,10] (as in the
%   default set-up); this means that 3 mono-objective optimization searches
%   will be run, using as demand the nominal demand of 'gasNetwork'
%   multiplied by the 'meanGrowth' relative to those years, i.e., for 
%   'refYears'=[0,5,10], values in meanGrowth([1,6,11]) will be used.
%
%   ANALYSIS public properties:
%      gasNetwork - 'gasnet' object detailing properties of the network
%      meanGrowth - yearly mean growth of demand for nYears
%      nScen - number of scenarios to be generated for each year
%      nYears - time horizon of the analysis
%      refYears - years whose demands will be used in the mono-obj. search
%      scenarios - all scenarios to be considered in the analysis
%      sdGrowth - std deviation of the demand growth given by meanGrowth
%
%   ANALYSIS set-access-private properties:
%      cheapest - structure with group of cheapest solutions:
%         pop  - (#cheapest_solutions x #pipes) matrix
%         fobj - (#cheapest_solutions x #merit_functions(4)) matrix
%         sens - (#nodes x #pipes x #cheapest_solutions) matrix
%      nondominated - structure with group of nondominated solutions:
%         pop  - (#nondominated_solutions x #pipes) matrix
%         fobj - (#nondominated_solutions x #merit_functions(4)) matrix
%         sens - (#nodes x #pipes x #nondominated_solutions) matrix
%      solutions - struct. with all solutions found via mono-obj. searches
%         pop  - (#solutions x #pipes) matrix
%         fobj - (#solutions x #merit_functions(4)) matrix
%
%   ANALYSIS public methods:
%      analyzesolutions - execute multi-objective robustness analysis
%      generatescenarios - generate scenarios (demand growth) to be used
%      generatesolutions - mono-objective - cost-minimizing - search 
%      plotcheapnd - plot solutions, highlighting cheapest and nondominated
%      plotmeritfunctions - plot two merit functions (f1 x f2)
%      plotsensitivitymatrix - plot sensitivity matrix
%      runanalysis - run entire analysis procedure: generation of solutions
%                    and multi-objective analysis
%
%   ANALYSIS private methods:
%      evalcriteria - evaluate 'feasibility' and 'fault cost' of solutions
%      evalsensitivity - evaluate 'sensitivity' of solutions
%      finddominated - find indices of dominated solutions
%
%   Example:
%      load CaseStudy;
%      cs_analysis = analysis;
%      cs_analysis.gasNetwork = cs_gasnet;
%      cs_analysis = cs_analysis.runanalysis;
%      cs_analysis.plotcheapnd;
%      cs_analysis.plotmeritfunctions(cs_analysis.nondominated,2,1);
%      cs_analysis.plotsensitivitymatrix(cs_analysis.nondominated,1);
%
%   References:
%      [1] Ramos, Eduardo S., Batista, Lucas S.. Natural Gas Pipeline Network 
%          Expansion Under Load-Evolution Uncertainty Using an Evolutionary
%          Algorithm, 2018
%
%   Authors: Eduardo S. Ramos    - edusantiago.ramos@gmail.com
%            Lucas S. Batista    - lusoba@cpdee.ufmg.br

%#ok<*AGROW>

classdef analysis
   properties
      meanGrowth = [1 repmat(1.025,1,10)];   % (1x(nYears+1)) vector
      sdGrowth = [1.002 repmat(1.012,1,10)]; % (1x(nYears+1)) vector
      scenarios = []; % (nScen*(nYears+1) x 1) vector
      nScen = 100;    % numeric
      nYears = 10;    % natural number
      
      refYears = [0 5 10]; % row vector of natural numbers
      
      gasNetwork = []; % object of class 'gasnet'
   end
   properties (SetAccess=private)
      solutions = [];    % structure
      nondominated = []; % structure
      cheapest = [];     % structure
   end %properties
   
   methods
      function obj = set.gasNetwork(obj,network)
         obj.gasNetwork = network.inputtest;
      end
      
      function obj = runanalysis(obj)
         %RUNANALYSIS Run entire analysis procedure.
         %
         %	 Prototype:
         %      obj = obj.runanalysis;
         %
         %   Input arguments:
         %      obj: 'analysis' object
         %
         %   Output arguments:
         %      obj: 'analysis' object
         
         if isempty(obj.gasNetwork)
            error('Please select gas network for analysis');
         end
         if isempty(obj.scenarios)
            obj.scenarios = obj.generatescenarios;
         end
         obj.solutions = obj.generatesolutions;
         obj = obj.analyzesolutions;
      end
      
      function scens = generatescenarios(obj)
         %GENERATESCENARIOS Generate scenarios (demand growth) to be used.
         %
         %	 Prototype:
         %      scens = obj.generatescenarios;
         %
         %   Input arguments:
         %      obj: 'analysis' object
         %
         %   Output arguments:
         %      scens: vector with scenarios - nominal-demand multipliers
         
         m = 1;
         s = 0;
         n = obj.nScen;
         
         % LATIN HYPERCUBE
         scens = zeros((obj.nYears+1)*n,1);
         for i=1:obj.nYears+1
            m = m*obj.meanGrowth(i);
            s = (1+s)*obj.sdGrowth(i) - 1;
            newScens = m + s.*lhsnorm(0,1,obj.nScen);
            scens(n*(i-1)+1:n*i) =  newScens';
            
            % Uncomment to plot pdfs:
%             minx = min(newScens); 
%             maxx = max(newScens); 
%             step = (maxx-minx)/10000;
%             pdf  = normpdf(minx:step:maxx,m,s);
%             plot(minx:step:maxx, pdf, 'b-', 'LineWidth', 2); hold on
         end
      end
      
      function allSol = generatesolutions(obj)
         %GENERATESOLUTIONS Mono-obj. search for cost-efficient solutions
         %
         %	 Prototype:
         %      allSol = obj.generatesolutions;
         %
         %   Input arguments:
         %      obj: 'analysis' object
         %
         %   Output arguments:
         %      allSol: structure with solutions and their costs
         %
         %   Example:
         %      obj.solutions = obj.generatesolutions;
         
         nExec = length(obj.refYears);
         tempNets = repmat(obj.gasNetwork,nExec,1);
         m = cumprod(obj.meanGrowth); m = m(obj.refYears+1);
         
         % Adjust network for each year
         for i=1:nExec 
            tempNets(i).demands = tempNets(i).demands.*m(i);
         end
         
         % Search for solutions in each network; 
         % indSol: individual solutions
         parfor i=1:nExec
            [indSol(i).pop,indSol(i).fobj(1,:)] = tempNets(i).evalnet;
         end
         
         % Join solutions
         allSol.pop = []; allSol.fobj = [];
         for i=1:nExec
            allSol.pop = [allSol.pop, indSol(i).pop];
            allSol.fobj = [allSol.fobj, indSol(i).fobj];
         end
         
         % Remove infeasible solutions
         ref = obj.gasNetwork;
         ref = ref.inputtest;
         [allSol.pop,allSol.fobj(1,:)] = ref.removeinf...
               (allSol.pop,allSol.fobj(1,:),ref.evalParams,ref.netData);
         
         % New convention: each row is a solution
         allSol.pop = allSol.pop'; allSol.fobj = allSol.fobj';
            
      end
      
      function obj = analyzesolutions(obj)
         %ANALYZESOLUTIONS Multi-objective robustness analysis of solutions
         %
         %	 Prototype:
         %      obj = obj.analyzesolutions;
         %
         %   Input arguments:
         %      obj: 'analysis' object
         %
         %   Output arguments:
         %      obj: 'analysis' object
         
         if isempty(obj.gasNetwork.evalParams)
            obj.gasNetwork = obj.gasNetwork.inputtest();
         end
         
         fprintf('\nAnalyzing solutions...\n');
         f = zeros(size(obj.solutions.pop,1),3);
         sols = obj.solutions.pop;

         % Evaluate objective functions
         [f(:,1),f(:,2)] = obj.evalcriteria(sols);
         [sensMatrix,f(:,3)] = obj.evalsensitivity(sols);

         % Save merit function values
         obj.solutions.fobj(:,2:4) = f;
                  
         % Save cheapest solutions
         nCheapSol = 5; % number of cheapest solutions to save
         obj.cheapest = [];
         [~,ord] = sort(obj.solutions.fobj(:,1));
         ord = ord(1:nCheapSol);
         obj.cheapest.pop = sols(ord,:);
         obj.cheapest.fobj = obj.solutions.fobj(ord,:);
         obj.cheapest.sens = sensMatrix(:,:,ord);
                  
         % Remove dominated solutions and save nondominated to variable
         obj.nondominated = [];
         dom = obj.finddominated(obj.solutions.fobj);
         obj.nondominated.pop = obj.solutions.pop(~dom,:);
         obj.nondominated.fobj= obj.solutions.fobj(~dom,:);
         obj.nondominated.sens = sensMatrix(:,:,~dom);

      end
      
      function plotcheapnd(obj)
         %PLOTCHEAPND Plot solutions highlighting cheapest and nondominated
         %
         %	 Prototype:
         %      obj.plotcheapnd;
         %
         %   Input arguments:
         %      obj: 'analysis' object
         
         % Plot nondominated
         nSols = size(obj.solutions.pop,1);
         sols = 1:nSols;
         dom = obj.finddominated(obj.solutions.fobj);
         figure;
         plot(sols(dom), obj.solutions.fobj(dom,1),'ok'); 
         hold on
         plot(sols(~dom),obj.solutions.fobj(~dom,1),'ok',...
               'MarkerFaceColor','k');
         xlabel('Solution Number');
         ylabel('Cost ($)');
         title('All solutions with nondominated ones highlighted');
         
         % Plot cheapest
         [~,cheapInds] = sort(obj.solutions.fobj(:,1));
         cheapInds = cheapInds(1:5);
         figure;
         plot(sols, obj.solutions.fobj(:,1),'ok'); 
         hold on
         plot(cheapInds,obj.solutions.fobj(cheapInds,1),'ok',...
               'MarkerFaceColor','k');
         xlabel('Solution Number');
         ylabel('Cost ($)');
         title('All solutions with cheapest ones highlighted');
      end
            
      function plotmeritfunctions(~,group,critx,crity)
         %PLOTMERITFUNCTIONS Plot comparison between two merit functions
         %
         %	 Prototype:
         %      obj.plotmeritfunctions(group,critx,crity);
         %
         %   Input arguments:
         %      group: group of solutions
         %      critx: merit function number of x variable
         %      crity: merit function number of y variable
         %
         %   Examples:
         %      obj.plotmeritfunctions(obj.cheapest,2,1);
         %      obj.plotmeritfunctions(obj.nondominated,3,4);
         
         if ~(ismember(critx,(1:4)) && ismember(crity,(1:4)))
            error('Criteria numbers must be in [1 4]');
         end
         
         labels = {'Cost ($)','Feasibility (%)','Fault Cost ($)',...
                   'Sensitivity Index'};
         multiplier = [1 100 1 1];
         num4text = num2str((1:size(group.fobj,1))');
         figure;
         plot(group.fobj(:,critx)*multiplier(critx),...
               group.fobj(:,crity)*multiplier(crity),'*r', 'MarkerSize',9);
         xlabel(labels{critx});
         ylabel(labels{crity});
         text(group.fobj(:,critx)*multiplier(critx),...
              group.fobj(:,crity)*multiplier(crity),num4text,...
             'HorizontalAlignment','left','VerticalAlignment','bottom');
      end
      
      function plotsensitivitymatrix(~,group,solNo)
         %PLOTSENSITIVITYMATRIX Plot sensitivity matrix of a solution
         %
         %	 Prototype:
         %      obj.plotmeritfunctions(group,critx,crity);
         %
         %   Input arguments:
         %      group: group of solutions
         %      solNo: index of solution
         %
         %   Examples:
         %      obj.plotsensitivitymatrix(obj.cheapest,1);
         %      obj.plotsensitivitymatrix(obj.nondominated,10);
         
         % Sensitivity indices < 1 are disregarded 
         sensMatrix = group.sens(:,:,solNo);
         [i,j] = find(log10(sensMatrix)>0);
         val = sensMatrix(i+(j-1)*size(sensMatrix,1));
         
         % Plot figure
         figure;
         scatter(j,i,64,log10(val),'filled','MarkerEdgeColor','k');
         map = colormap(gray);
         colormap(flipud(map));
         cb = colorbar; 
         set(cb,'Limits',[log10(min(val(val>1))) log10(max(val))]);
         set(cb,'TicksMode','manual');
         set(cb,'Ticks',linspace(cb.Limits(1),cb.Limits(2),10));
         try % versions R2014b on
            cb.Label.String = 'Sensitivity Index - log_{10} scale';
         catch
         end
         set(gca,'XTick',1:size(sensMatrix,2));
         set(gca,'YTick',1:size(sensMatrix,1));
         xlim([0 size(sensMatrix,2)+1])
         ylim([0 size(sensMatrix,1)+1])
         xlabel('Pipe Number');
         ylabel('Node Number');
         title('Sensitivity Matrix - log_{10} scale');
      end
      
   end %methods
   
   methods (Access = private)
      function [feasRate,mFaultCost] = evalcriteria(obj,sols) 
         %EVALCRITERIA Evaluate feasibility and fault cost of solutions
         %
         %	 Prototype:
         %      [feasRate,mFaultCost] = evalcriteria(obj,sols) 
         %
         %   Input arguments:
         %      obj: 'analysis' object
         %      sols: matrix with solutions
         %
         %   Output arguments:
         %      feasRate: feasibility rate of solutions
         %      mFaultCost: mean fault cost of solutions
         
         params = obj.gasNetwork.evalParams; % parameters for evaluation
         netData = obj.gasNetwork.netData;   % network data
         scens = obj.scenarios;              % scenarios
         
         nIntNodes = params.ni; % number of internal nodes
         nPipes = params.np;    % number of pipes

         nSolutions = size(sols,1);   % number of solutions
         nScenarios = length(scens);  % number of scenarios

         failureRate=ones(1,params.np).*1e-4;% mean failure rate of pipes
         faultDuration = ones(1,params.np);  % mean fault duration of pipes
         gas_tx = 0.10;                      % gas tax

         feasRate = zeros(nSolutions,1);     % infeasibility rates
         mFaultCost = zeros(nSolutions,1);   % mean fault costs
         for i = 1:nSolutions
            fprintf('\nCriteria: Solution %d/%d\n',i,nSolutions);
            currentSol = sols(i,:)';          % solution to be analyzed
            inflag = zeros(1,nScenarios);    
            faultCost = zeros(1,nScenarios); 
            parfor j = 1:nScenarios   
               % Parameters
               par = params;
               par.Q = scens(j)*par.Q;
               
               % Network data
               net = netData;
               net.qc = net.P * [net.F*par.Q; zeros(nPipes-nIntNodes,1)];
                
               % Calculate the network flow for the solution i
               [q,pi,~] = gasnet.networkflow(currentSol,par,net);

               % Determine the number of nodes that violate min pressure
               w = length(pi(pi<par.pmin));

               % Infeasibility flag
               inflag(j) = min(1,w);        

               % Fault cost of solution i
               % [Edit] Eq.(11) Carrano
               faultCost(j) = sum(failureRate .* faultDuration .* ...
                                     par.L(:)' .* abs(q(:)') * gas_tx);
               faultCost(j) = not(inflag(j))*faultCost(j);
            end

            % Feasibility rate 
            nInfeasible = sum(inflag);
            nFeasible = nScenarios - nInfeasible;
            feasRate(i) = nFeasible/nScenarios;

            % Mean fault cost of feasible solutions
            mFaultCost(i) = sum(faultCost)/nFeasible;
         end
         fprintf('\n');
      end 
      
      function [sensMatrix,sensRate] = evalsensitivity(obj,sols)
         %EVALSENSITIVITY Evaluate sensitivity index of solutions
         %
         %	 Prototype:
         %      [sensMatrix,sensRate] = evalsensitivity(obj,sols)
         %
         %   Input arguments:
         %      obj: 'analysis' object
         %      sols: matrix with solutions
         %
         %   Output arguments:
         %      sensMatrix: element sensMatrix(n,i,j): sensitivity of node
         %                  'n' of solution 'j' when pipe 'i' is removed
         %      sensRate: sensitivity index of solutions
         
         params = obj.gasNetwork.evalParams; % parameters for evaluation
         
         nSolutions = size(sols,1);  % number of solutions
         nPipes = params.np;         % number of pipes
         nInternalNodes = params.ni; % number of internal nodes

         expectedGrowth = cumprod(obj.meanGrowth); % growth in demand
         ref = obj.gasNetwork; % reference network
         
         sensMatrix = zeros(nInternalNodes,nPipes,nSolutions);
         for year = 0:obj.nYears % for each year in our time horizon
            fprintf('\nSensitivity: Year %d/%d\n',year,obj.nYears);
            
            % Mean demand for that year
            ref.demands = obj.gasNetwork.demands*expectedGrowth(year+1);
            for i = 1:nPipes % for each pipe 
               fprintf('\nSensitivity: Pipe %d/%d\n',i,nPipes);

               % Get original network, remove pipe 'i' and get parameters
               ref2 = ref;
               ref2.pipes(i,:) = [];    
               ref2.pipeLengths(i) = []; 
               ref2 = ref2.inputtest;
               params = ref2.evalParams;
               net = ref2.netData;
               newSols = sols'; newSols(i,:) = []; % solutions without pipe

               parfor j = 1:nSolutions 
                  % Calculate the network flow for the modified solution j
                  [~,pi,~] = gasnet.networkflow(newSols(:,j),params,net);

                  inf = pi<params.pmin;
                  pi(~inf) = 0;
                  pi(inf) = ((params.pmin(inf)-pi(inf)).^2)/...
                                                        max(sum(~inf),0.5);

                  % Sum all years
                  sensMatrix(:,i,j) = sensMatrix(:,i,j) + pi;
               end
            end    
         end
         sensMatrix = sensMatrix./(obj.nYears+1); % average of the values
         
         sensRate = sum(sensMatrix,1);          % sum difference per node
         sensRate = squeeze(sum(sensRate,2));   % sum difference per pipe
         sensRate = sensRate(:)';
      end 
      
      function inds = finddominated(~,popFobj)
         %FINDDOMINATED Find indices of dominated solutions
         %
         %	 Prototype:
         %      inds = finddominated(~,popFobj)
         %
         %   Input arguments:
         %      popFobj: matrix with merit functions of solutions
         %
         %   Output arguments:
         %      inds: indices of dominated solutions
         
         [nSols, nCrit] = size(popFobj);
         popFobj(:,2)= 1-popFobj(:,2); % feasibility rate: smaller is worse
         
         % I: comparison matrix
         popT = popFobj';
         I = repmat(popFobj,1,nSols) >= repmat(popT(:)',nSols,1);
         
         % I(i,j,k)=1 if sol. k is worse than j with respect to criterion i
         I = reshape(I',nCrit,nSols,nSols);
         
         % I(1,j,k)=1 if sol. k is dominated by j
         I = prod(double(I),1);
         
         % Dominated solutions
         inds = sum(reshape(I(:,:),nSols,nSols),1) > 1; 
         inds = inds';
      end
      
   end %methods(private)
   
end %classdef analysis