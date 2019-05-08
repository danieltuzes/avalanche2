

#include "stdafx.h"
#ifdef IS_GPU
#include "snapshot.h"
#include "make_simulation.h" //version information is stored at make_simulation.h
#define numOfThrPerBlock 512 //maximum value of 512 for cuda compute capabiltiy 1.3 or below, 1024 for 2.0 or above

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>



//defines the variables used on the device
class GPU_vars
{
public:
	GPU_vars()
	{
		kernel = NULL;
		tau_l = NULL;
		tau_n = NULL;
		tau_p = NULL;
		pos = NULL;
	}

	deftype * kernel;
	deftype * tau_l;
	deftype * tau_n;
	deftype * tau_p;
	int * pos;
	place * minPoint;
	deftype * deform;
};


//allocate the required memory on the device and copy the required variables to the device
bool allocateAndCopyToDevice(simVars& sim, GPU_vars& d)
{
	int size = sim.getSize();
	int linsize = size * size;

	cudaError_t cudaStatus;
	if ((cudaStatus = cudaMalloc(&d.kernel,sizeof(deftype) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d.tau_l,sizeof(deftype) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d.tau_n,sizeof(deftype) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d.tau_p,sizeof(deftype) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d.pos,sizeof(int) * linsize)) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d.minPoint,sizeof(place))) != cudaSuccess ||
		(cudaStatus = cudaMalloc(&d.deform,sizeof(deftype))) != cudaSuccess)
	{
		cerr << "cudaMalloc failed, " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	if ((cudaStatus = cudaMemcpy(d.kernel,&sim.sf[0],sizeof(deftype) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess ||
		(cudaStatus = cudaMemcpy(d.tau_l,&sim.tau_l[0],sizeof(deftype) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess ||
		(cudaStatus = cudaMemcpy(d.tau_n,&sim.tau_n[0],sizeof(deftype) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess ||
		(cudaStatus = cudaMemcpy(d.tau_p,&sim.tau_p[0],sizeof(deftype) * linsize,cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		cerr << "cudaMemcpyHostToDevice failed, " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	return true;
}


//copy back the tau_l space
//tau_l cannot be in a synchronised state with host device in an effective way
bool tau_lCopyToHost(simVars& sim, GPU_vars& d)
{
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(&sim.tau_l[0],d.tau_l,sizeof(deftype) * sim.getSize() * sim.getSize(),cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaMemcpy(&sim.tau_l,d.tau_l,sizeof(int) * sim.getSize() * sim.getSize(),cudaMemcpyDeviceToHost) failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	return true;
}


//copy back the tau_n and tau_p space
//these variables can be in a synchronised state with the host device
bool tau_npCopytoHost(simVars& sim, GPU_vars& d)
{
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(&sim.tau_n[0],d.tau_n,sizeof(deftype) * sim.getSize() * sim.getSize(),cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		cerr << "cudaMemcpy(&sim.tau_n,d.tau_n,sizeof(int) * sim.getSize() * sim.getSize(),cudaMemcpyDeviceToHost) failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaMemcpy(&sim.tau_p[0],d.tau_p,sizeof(deftype) * sim.getSize() * sim.getSize(),cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		cerr << "cudaMemcpy(&sim.tau_p,d.tau_p,sizeof(int) * sim.getSize() * sim.getSize(),cudaMemcpyDeviceToHost) failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	return true;
}

__global__ void resetPos(int * pos)
{
	const int id = threadIdx.x + blockIdx.x * blockDim.x;
	pos[id] = id;
}

/// <summary>Find minimum values in array of size of threadPerblock</summary>
/// <param name="linsize">The number of elemets that has to compared.
/// It is not sure that every elemets will be compared to every otther one,
/// but minimum values will be selected from every threadPerblock number of elements.</param>
/// <param name="threadPower">The lb of the minimum distance - 1 between two cells that has to compared. It is equal with lb(sim.getSize() * sim.getSize()  / linsize) - 1.</param>
__global__ void findMinIterate(int threadPower, deftype tau_ext, deftype * tau_l, deftype * tau_n, deftype * tau_p, int * pos)
{
	const int id = (threadIdx.x + (blockIdx.x * blockDim.x << 1)) << threadPower;

	//find the minimum place
	for (int i = blockDim.x; i > 0; i>>=1)
	{
		if (threadIdx.x < i)
		{
			deftype tau_act_A, tau_act_A_l, tau_act_B, tau_act_B_l;
			
			tau_act_A   = tau_p[pos[id]] - tau_ext - tau_l[pos[id]];
			tau_act_A_l = tau_n[pos[id]] + tau_ext + tau_l[pos[id]];
			if (tau_act_A_l < 0 && tau_act_A_l < tau_act_A)
				tau_act_A = tau_act_A_l;

			tau_act_B   = tau_p[pos[id + (i<<threadPower)]] - tau_ext - tau_l[pos[id + (i<<threadPower)]];
			tau_act_B_l = tau_n[pos[id + (i<<threadPower)]] + tau_ext + tau_l[pos[id + (i<<threadPower)]];
			if (tau_act_B_l < 0 && tau_act_B_l < tau_act_B)
				tau_act_B = tau_act_B_l;

			if (tau_act_B < tau_act_A)
				pos[id] = pos[id + (i<<threadPower)];

		}
		__syncthreads();
	}
}

__global__ void setMinPointResetPos(int size, int sizePower, deftype tau_ext, deftype * tau_l, deftype * tau_n, deftype * tau_p, int * pos, place * minPoint)
{
	const int id = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (id == 0)
	{
		minPoint->tau_l = tau_l[pos[0]];
		minPoint->x = pos[0] >> sizePower;
		minPoint->y = pos[0] & (size-1);
	}
	pos[id] = id;
}

__global__ void addKernelRefreshYield(int size, int power, deftype * kernel, deftype * tau_l, deftype * tau_n, deftype * tau_p, place * minPoint, direction dir, deftype newYield, deftype deform)
{
	const int id = threadIdx.x + blockIdx.x * blockDim.x;
	
	//x and y values are needed to correctly referr to kernel values
	const int prevX = minPoint->x;
	const int prevY = minPoint->y;
	int prevPos = (prevX << power) + prevY;

	const int idx = id >> power;
	const int idy = id & (size-1);

	//refresh the Yield point at the previous deformation
	if (id == prevPos)
	{
		if (dir == direction::left)
			tau_n[prevPos] = newYield;
		else
			tau_p[prevPos] = newYield;
	}

	tau_l[id] += kernel[(((idx-prevX+size) & (size - 1)) << power) + ((idy-prevY+size) & (size - 1))] * deform;
}




//find the minimum points, where deformation incidence will appear
//and also its direction and value
bool findMin(simVars& sim, place& minPoint, GPU_vars& d)
{
	resetPos<<<sim.getSize() * sim.getSize() / numOfThrPerBlock,numOfThrPerBlock>>>(d.pos);
	cudaError_t cudaStatus;
	
#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "resetPos launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching resetPos!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion
	int threadPower = 0;
	const int iterationPower = getPower(numOfThrPerBlock) + 1;
	const int linSize = sim.getSize() * sim.getSize();
	for (int restSize =  linSize; restSize > 1; restSize /= numOfThrPerBlock * 2)
	{
		int numOfBlock = restSize / (2 * numOfThrPerBlock);
		if (numOfBlock == 1 || numOfBlock == 0)
		{
			findMinIterate<<<1,restSize/2>>>(threadPower, sim.tau_ext, d.tau_l, d.tau_n, d.tau_p, d.pos);
#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "findMinIterate<<<1,numOfThrPerBlock>>> launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching findMinIterate<<<1,numOfThrPerBlock>>>!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion
			break;
		}
		findMinIterate<<<numOfBlock,numOfThrPerBlock>>>(threadPower, sim.tau_ext, d.tau_l, d.tau_n, d.tau_p, d.pos);
#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "findMinIterate<<<numOfBlock,numOfThrPerBlock>>> launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching findMinIterate<<<numOfBlock,numOfThrPerBlock>>>!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion
		threadPower += iterationPower;
	}

	setMinPointResetPos<<<sim.getSize() * sim.getSize() / numOfThrPerBlock,numOfThrPerBlock>>>(sim.getSize(), sim.getPower(), sim.tau_ext, d.tau_l, d.tau_n, d.tau_p, d.pos, d.minPoint);

#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "setMinPointResetPos launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching setMinPointResetPos!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion

	cudaStatus = cudaMemcpy(&minPoint,d.minPoint,sizeof(place),cudaMemcpyDeviceToHost);

#pragma region check cudaStatus
#ifdef _DEBUG
	if (cudaStatus != cudaSuccess)
	{
		cerr << "cudaMemcpy(&minPoint,d.minPoint,sizeof(place),cudaMemcpyDeviceToHost) failed, " << cudaGetErrorString(cudaStatus) << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion

	return true;
}


//at minPoint, calclulates the needed deformation based on tau_l and
//             refresh the tau_l field with the previous deformation placed at minPoint
//than finds the minimum point and
//refresh the yield point of the deformed cell
bool refreshFindMin(const simVars& sim, place& minPoint, direction nowFlow, deftype deform, GPU_vars& d)
{
	cudaError_t cudaStatus;
	
	deftype newYield = (nowFlow == direction::right) ? sim.tau_p(minPoint.x,minPoint.y) : sim.tau_n(minPoint.x,minPoint.y);
	addKernelRefreshYield<<<sim.getSize() * sim.getSize() / numOfThrPerBlock,numOfThrPerBlock>>>(
		sim.getSize(),
		sim.getPower(),
		d.kernel, d.tau_l, d.tau_n, d.tau_p, d.minPoint,
		nowFlow,
		newYield,
		deform);

#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "addKernelRefreshYield launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching addKernelRefreshYield!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion

	int threadPower = 0;
	const int iterationPower = getPower(numOfThrPerBlock) + 1;
	const int linSize = sim.getSize() * sim.getSize();
	for (int restSize =  linSize; restSize > 1; restSize /= numOfThrPerBlock * 2)
	{
		int numOfBlock = restSize / (2 * numOfThrPerBlock);
		if (numOfBlock == 1 || numOfBlock == 0)
		{
			findMinIterate<<<1,restSize/2>>>(threadPower, sim.tau_ext, d.tau_l, d.tau_n, d.tau_p, d.pos);
#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "findMinIterate<<<1,numOfThrPerBlock>>> launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching findMinIterate<<<1,numOfThrPerBlock>>>!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion
			break;
		}
		findMinIterate<<<numOfBlock,numOfThrPerBlock>>>(threadPower, sim.tau_ext, d.tau_l, d.tau_n, d.tau_p, d.pos);
#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "findMinIterate<<<numOfBlock,numOfThrPerBlock>>> launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching findMinIterate<<<numOfBlock,numOfThrPerBlock>>>!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion
		threadPower += iterationPower;
	}

	setMinPointResetPos<<<sim.getSize() * sim.getSize() / numOfThrPerBlock,numOfThrPerBlock>>>(sim.getSize(), sim.getPower(), sim.tau_ext, d.tau_l, d.tau_n, d.tau_p, d.pos, d.minPoint);

#pragma region check cudaStatus
#ifdef _DEBUG
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cerr << "setMinPointResetPos launch failed: " << cudaGetErrorString(cudaStatus) << endl;;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaDeviceSynchronize returned error code " << cudaGetErrorString(cudaStatus) << " after launching setMinPointResetPos!" << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion

	cudaStatus = cudaMemcpy(&minPoint,d.minPoint,sizeof(place),cudaMemcpyDeviceToHost);

#pragma region check cudaStatus
#ifdef _DEBUG
	if (cudaStatus != cudaSuccess)
	{
		cerr << "cudaMemcpy(&minPoint,d.minPoint,sizeof(place),cudaMemcpyDeviceToHost) failed, " << cudaGetErrorString(cudaStatus) << endl;
		cudaFree(d.kernel);
		cudaFree(d.tau_l);
		cudaFree(d.tau_n);
		cudaFree(d.tau_p);
		cudaFree(d.pos);
		cudaFree(d.minPoint);
		cudaFree(d.deform);
		return false;
	}
#endif
#pragma endregion
	
	return true;

}


bool simVars::makeDeformSerial_GPU(const simPars& pars, int nominalDef, bool printToTau_g)
{

	time_t lastWrite = time(NULL); 

	place minPoint;

	GPU_vars d;

#pragma region device_setup

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaError_t cudaStatus;
    if (cudaSetDevice(0) != cudaSuccess) {
        cerr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?" << endl;;
        return false;
    }
	// cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    if (cudaDeviceReset() != cudaSuccess) {
        cerr << "cudaDeviceReset failed!" << endl;;
        return false;
    }
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		cerr << "cudaDeviceSynchronize failed, " << cudaGetErrorString(cudaStatus) << endl;
		return false;
	}
	if (cudaThreadSynchronize() != cudaSuccess)
	{
		cerr << "cudaThreadSynchronize failed" << endl;
		return false;
	}
	
	if (!allocateAndCopyToDevice(*this,d))
	{
		cerr << "Error at allocateAndCopyToDevice" << endl;
		return false;
	}
#pragma endregion
	
	if (!findMin(*this,minPoint,d))
	{
		cerr << "Error at findMin" << endl;
		return false;
	}

	printInfDescr(pars);

	ofstream tau_g(fnameWp("tau_g",pars),ios_base::app);
	if (!tau_g)
	{
		cerr << "Cannot app to " << fnameWp("tau_g",pars) << endl;
		return false;
	}
	tau_g.precision(17); //http://en.wikipedia.org/wiki/Double-precision_floating-point_format#IEEE_754_double-precision_binary_floating-point_format:_binary64
	while (getDeform() < nominalDef)
	{
		if (DEBUG_TEST)
			printInf(minPoint,*this,pars);

		if (willYield(*this,minPoint) == none) //elõfordulhatna amúgy, h a betöltéskor épp egy lavina van, és akkor nem kell megnövelni; tehát ha egy lavina elején vagyunk
		{
			tau_ext =  tau_p(minPoint.x,minPoint.y) - minPoint.tau_l + EPSILON;
			av_size = 0;
			av_size_n = 0;
		}

		do
		{
			if (static_cast<int>(time(NULL) - lastWrite) > pars.tbs()) //Gamma and tau_l are in consistent state, not as like at the end of this loop
			{
				tau_lCopyToHost(*this,d); //copy back the data from the device
				tau_npCopytoHost(*this,d);
				snapshot bak(*this,static_cast<int>(getDeform()));
				bak.make(pars);
				delMarkedSnapshots(pars);
				bak.mark(pars);
				lastWrite = time_t(0);
			}
			
			direction nowFlow = willYield(*this,minPoint);
			deftype deform = calcDeform(minPoint,tau_ext,sf[0],pars.DG(),nowFlow);
			Gamma[minPoint.x * size + minPoint.y] += deform;
			if (nowFlow == direction::right)
			{
				av_size += deform;
				sG += deform;
			}
			else // nowFlow == direction::left must hold
			{
				av_size_n -= deform; //av_size_n stay positive
				sGn -= deform; //sGn stay positive
			}

			if (pars.fsrg())
			{
				if (nowFlow == direction::right)
				{
					tau_p(minPoint.x, minPoint.y) = rand_gen();
					if (pars.lfs().MyType == leftFlowStress::same)
						tau_n(minPoint.x, minPoint.y) = tau_p(minPoint.x, minPoint.y);
				}
				else // nowFlow == direction::left
				{
					tau_n(minPoint.x, minPoint.y) = rand_gen();
					if (pars.lfs().MyType == leftFlowStress::same)
						tau_p(minPoint.x, minPoint.y) = tau_n(minPoint.x, minPoint.y);
				}
			}
				
			if (!refreshFindMin(*this,minPoint,nowFlow,deform,d))
			{
				cerr << "Error at refreshFindMin" << endl;
				return false;
			}

			if (DEBUG_TEST)
				printInf(minPoint,*this,pars);

		} while (willYield(*this,minPoint) != direction::none && getDeform() < nominalDef);

		if (willYield(*this,minPoint) == direction::none && printToTau_g)
			printTau_g(tau_g);
	}
	
	tau_lCopyToHost(*this,d);
	tau_npCopytoHost(*this,d);
	if (!snapshot(*this,static_cast<int>(getDeform())).make(pars) || !delMarkedSnapshots(pars))
		return false;
	
	cudaFree(d.kernel);
	cudaFree(d.tau_l);
	cudaFree(d.tau_n);
	cudaFree(d.tau_p);
	cudaFree(d.pos);
	cudaFree(d.minPoint);
	cudaFree(d.deform);

	return true;
}

#endif