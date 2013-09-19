#include <iostream>
#include <stdio.h>

int main()
{
	cudaDeviceProp properties;
	int count;

	cudaGetDeviceCount(&count);

	for(int i = 0; i < count; i++)
	{
		cudaGetDeviceProperties(&properties,i);
		std::cout << "------------ GENERAL INFORMATION ------------" << std::endl << std::endl;
		std::cout << "Name: " << properties.name << std::endl;
		std::cout << "Compute Capabilities: " << properties.major << "." <<  properties.minor << std::endl;
		std::cout << "Clock Rate: " << properties.clockRate << std::endl;
		std::cout << "Number of Multiprocessors: " << properties.multiProcessorCount << std::endl;
		std::cout << "Max Threads per Multiprocessor: " << properties.maxThreadsPerMultiProcessor << std::endl;
		std::cout << "Asynchronous Engine Count: ";
		if (properties.deviceOverlap == 2)	
			std::cout << "Both direction at a time" << std::endl;
		else if (properties.deviceOverlap == 1)
			std::cout << "One direction at a time" << std::endl;
		else if (properties.deviceOverlap == 0)
			std::cout << "Disabled" << std::endl;
		std::cout << "Kernel Execition TimeOut: ";
		if (properties.kernelExecTimeoutEnabled)	
			std::cout << "Enabled" << std::endl;
		else
			std::cout << "Disabled" << std::endl;
		std::cout << "Total Global Memory: " << properties.totalGlobalMem << std::endl;
		std::cout << "Shared Memory per Block: " << properties.sharedMemPerBlock << std::endl;
		std::cout << "Total Constant Memory: " << properties.totalConstMem << std::endl;
		std::cout << "Registers per Block: " << properties.regsPerBlock << std::endl;
		std::cout << "Max Threads per Block: " << properties.maxThreadsPerBlock << std::endl;
		std::cout << "Max Threads Dimensions: " << properties.maxThreadsDim[0] << " " << properties.maxThreadsDim[1] << " " << properties.maxThreadsDim[2] << std::endl;
		std::cout << "Max Grid Size: " << properties.maxGridSize[0] << " " << properties.maxGridSize[1] << " " << properties.maxGridSize[2] << std::endl;
	}
	return 0;
}
