Instructions for executing the aiGNet algorithm to find cost-efficient configurations for gas networks:

1. Instantiate an object of class 'gasnet' with information on the network.
	Obs.: Check the public properties of the class definition to see how to do so.

2. Execute the immune algorithm to generate viable configurations 
	>> [configurations, costs] = <gasnet_object>.evalnet


gasnet-object examples:
	SEE 'Benchmarks' MAT-file
Execution example:
	load Benchmarks;
	[configurations, costs] = berlin52a.evalnet;


Mind the parallel processing loops (parfor). For some versions of matlab, a matlabpool must be opened:
	if matlabpool('size')==0
		matlabpool open
	end
	(...)
	delete(poolobj)


