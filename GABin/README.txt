Instructions for executing the Binary Genetic Algorithm (GAbin) to find cost-efficient configurations for gas networks:

1. Create a structure and save it as a MAT-file with the network's info and place it inside the 'networks' folder.
	xy: 
	jnode:
	reservoir:
	inode:
	params:

2. Execute the genetic algorithm with function 'ga.m'
	>> [configurations, costs] = ga(<name_of_network_struct>, <no_of_individuals>)


MAT-file examples:
	SEE 'networks' folder
Execution example:
	[configurations, costs] = ga('berlin52a',4)


Mind the parallel processing loops (parfor). For some versions of matlab, a matlabpool must be opened:
	if matlabpool('size')==0
		matlabpool open
	end
