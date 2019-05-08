
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>

#include "stdafx.h"
#include "make_simulation.h"
#include "snapshot.h"

bool allocate_and_copy_sim_to_device(simVars& sim, double * d_kernel, double * d_tau_l, double * d_tau_n, double * d_tau_p, int * d_pos)
{
	int linsize = sim.getSize() * sim.getSize();

	cudaError_t cudaStatus;
	if ((cudaStatus = cudaMalloc(&d_kernel,sizeof(double) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d_tau_l,sizeof(double) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d_tau_n,sizeof(double) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d_tau_p,sizeof(double) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d_pos,sizeof(int) * linsize/2)) != cudaSuccess)
	{
		cerr << "cudaMalloc failed, " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d_kernel);
		cudaFree(d_tau_l);
		cudaFree(d_tau_n);
		cudaFree(d_tau_p);
		cudaFree(d_pos);
		return false;
	}

	if ((cudaStatus = cudaMemcpy(d_kernel,&sim.sf[0],sizeof(double) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess ||
		(cudaStatus = cudaMemcpy(d_tau_l,&sim.tau_l[0],sizeof(double) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess ||
		(cudaStatus = cudaMemcpy(d_tau_n,&sim.tau_n[0],sizeof(double) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess ||
		(cudaStatus = cudaMemcpy(d_tau_p,&sim.tau_p[0],sizeof(double) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		cerr << "cudaMemcpyHostToDevice failed, " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d_kernel);
		cudaFree(d_tau_l);
		cudaFree(d_tau_n);
		cudaFree(d_tau_p);
		cudaFree(d_pos);
		return false;
	}
}
__global__ void findMin(int linsize, double tau_ext, const double * kernel, const double * tau_l, const double * tau_n, const double * tau_p, int * pos)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id > linsize)
		return;

	pos[id] = id;

	for (int i=linsize/2; i>0; i/=2)
	{
		if (id < i)
		{
			double tau_act_A, tau_act_A_l, tau_act_B, tau_act_B_l;
			
			tau_act_A = tau_p[pos[id]] - tau_ext - tau_l[pos[id]];
			tau_act_A_l = tau_n[pos[id]] + tau_ext + tau_l[pos[id]];
			if (tau_act_A_l < 0 && tau_act_A_l < tau_act_A)
				tau_act_A = tau_act_A_l;

			tau_act_B = tau_p[pos[id+i]] - tau_ext - tau_l[pos[id+i]];
			tau_act_B_l = tau_n[pos[id+i]] + tau_ext + tau_l[pos[id+i]];
			if (tau_act_B_l < 0 && tau_act_B_l < tau_act_B)
				tau_act_B = tau_act_B_l;

			if (tau_act_B < tau_act_A)
				pos[id] = pos[id+i];

			__syncthreads();
		}
	}
}
bool getPosition(double * d_kernel, double * d_tau_l, double * d_tau_n, double * d_tau_p, int * d_pos, int& pos);