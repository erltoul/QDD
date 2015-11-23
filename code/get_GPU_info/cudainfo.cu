/*
This file is a part of PW-TELEMAN project.
PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.

PW-Teleman is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PW-Teleman is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.
*/

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
