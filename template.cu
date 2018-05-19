/*

 
Use "REAL" as floating point type instead of double or float
 
Compile with optional flag:
    -DDOUBLE to use double instead of  float
     requires also -arch sm_21
 
 
single precision code:
 
> nvcc template.cu fileutils.cpp stringutils.cpp graphicstools.cpp -lcufft -o demo_single
 
double precision code:
 
> nvcc template.cu fileutils.cpp stringutils.cpp graphicstools.cpp -lcufft -DDOUBLE -arch sm_21 -o demo_single

 
*/

#include <stdio.h>
#include <stdlib.h> /* for rand() */
#include <unistd.h> /* for getpid() */
#include <time.h> /* for time() */
#include <math.h>
#include <assert.h>

#include <cuda.h>
#include <cufft.h>

#include "fileutils.h"
#include "stringutils.h"
#include "graphicstools.h"

// ******************************************************

#define PI	3.1415926535897932384626433832795
#define TWOPI 6.28318530717958647692528676655901

// construct REAL "type," depending on desired precision
// set the maximum number of threads

#ifdef DOUBLE
 #define REAL double
 #define MAXT 256
#else
 #define REAL float
 #define MAXT 512
#endif

typedef struct {
	REAL re;
	REAL im;
} COMPLEX;

// ******************************************************
//calculate the k-index in order to determine the correct k-vector for a given x,y,z-index

#define k_INDEX(i,L) ((i)<=((L)/2)?(i):((i)-(L)))

// ******************************************************

//initialize a real GPU array with a constant
__global__ void G_setrealconst(int N,REAL *a,REAL val) {
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    if(idx<N) a[idx]=val;
};


//multiply two complex GPU arrays (A,B) and store result in A
__global__ void G_mulcarray(int N,COMPLEX *A,COMPLEX *B)
{
    int i=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    REAL re,im,re2,im2;
    if(i<N)
    {
        re=A[i].re;im=A[i].im;
        re2=B[i].re;im2=B[i].im;
        A[i].re=re*re2-im*im2;
        A[i].im=im*re2+re*im2;
    }
};


// ******************************************************


//execute the FFT on the GPU, zin and zout can be the same array for "in-place" FFT (a little slower)
//set "fwd" to false for inverse FFT
void G_FFT(COMPLEX *zin,COMPLEX *zout,cufftHandle &fftPlan,bool fwd=true)
{
#ifdef DOUBLE
    if(fwd) cufftExecZ2Z(fftPlan,(cufftDoubleComplex*) zin,(cufftDoubleComplex*) zout,CUFFT_FORWARD);
    else    cufftExecZ2Z(fftPlan,(cufftDoubleComplex*) zin,(cufftDoubleComplex*) zout,CUFFT_INVERSE);
#else
    if(fwd) cufftExecC2C(fftPlan,(cufftComplex*) zin,(cufftComplex*) zout,CUFFT_FORWARD);
    else    cufftExecC2C(fftPlan,(cufftComplex*) zin,(cufftComplex*) zout,CUFFT_INVERSE);
#endif
};

 
// ******************************************************

//split a complex array in two real arrays containing amplitude^2 and phase
__global__ void G_ampphase(int N,COMPLEX *A,REAL* amp2,REAL* phase)
{
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    REAL re,im;
    if(idx<N)
    {
        re=A[idx].re;im=A[idx].im;
        amp2[idx]=re*re+im*im;
        phase[idx]=atan2(im,re);
    }
};

//split a complex array in two real arrays containing real and imaginary parts
__global__ void G_splitreim(int N,COMPLEX *A,REAL* re,REAL* im)
{
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    if(idx<N)
    {
        im[idx]=A[idx].im;
        re[idx]=A[idx].re;
    }
};

// ******************************************************


//check for a CUDA error, use argument for identification
bool cerr(const char *s="n/a")
{
    cudaError_t err=cudaGetLastError();
    if(err==cudaSuccess)
        return false;
    printf("CUDA error [%s]: %s\n",s,cudaGetErrorString(err));
    return true;
};


//some function initializing a 2D complex array
__global__ void G_function(int Nx,int Ny, COMPLEX *f,REAL t) {
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    int i,j;
    if(idx<Nx*Ny) {
        i=idx%Nx;
        j=idx/Nx;
        
        f[idx].re=sin(0.1*i+t)*cos(0.1*t*j);
        f[idx].im=-sin(0.1*j+t)*cos(0.1*t*i);
    }
};


//function to calculate the x-derivate in Fourier space
__global__ void G_kernel(int Nx, COMPLEX *f,REAL dkx) {
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    int i,j,ki;
    REAL x,y,k;
    if(idx<Nx*Nx) {
        i=idx%Nx;
        j=idx/Nx;
        
        //calculate the x-derivative in Fourier space
        ki=k_INDEX(i,Nx); //the Fourier component
        //kj=k_INDEX(j,Ny);
        
        //multipy i*k_x to f(k_x,k_y)
        k=dkx*ki;
        y=k*f[idx].re;
        x=-k*f[idx].im;
        
        f[idx].re=x;
        f[idx].im=y;
    }
}

// ******************************************************

//outputs an NetPBM image based on a real array a, m is the minimum value of a and M the maximum, nx&ny the dimension of a and cgrad a color gradient
/* cgrad values
 0: rainbow
 1: rainbow 2
 2: rainbow 3
 3: dark rainbow
 4: temperature
 5: temperature 2
 6: thermo
 7: solar
 8: sunset
 9: neon
*/

void writeBM_real(string fn,REAL *a,REAL m,REAL M,int nx,int ny,int cgrad)
{
    int i,n;
    dcolor dcol;
    unsigned int col;
    REAL val,dci=1.0/(M-m);
    unsigned int *rgb;
    unsigned char *gray;
    PXMfile *Ifile;
    colorfunction *colors;
    
    n=nx*ny;
    
    colors= new colorfunction();
    colors->selectgradient(cgrad);
    
    rgb=new unsigned int[n];
    
    gray=(unsigned char *) rgb;
    for(i=0;i<n;i++)
    {
        val=(a[i]-m)*dci;
        if(cgrad<1) {
            col=(unsigned int) (256*val);if(col>255) col=255;
            gray[i]=col;}
        else {
            dcol=colors->getgradientcolor(val);col=colors->get32bitcolor(dcol);
            rgb[i]=col;}
    }
    
    Ifile=new PXMfile(fn,(cgrad<1?PXM_P5:PXM_P6));
    Ifile->setsize(nx,ny);
    if(cgrad<1) Ifile->writefile(gray,nx*ny);
    else Ifile->writefile(rgb,nx*ny);
    delete Ifile;
    
    delete[] rgb;
    delete colors;
};

// ******************************************************

int main(int argc,char *argv[])
{
    int N,i,n,dim;
    int threads,blocks;
    REAL t,dt,dkx,L,mval,Mval,x;
    size_t fmem,tmem;
    COMPLEX *GF,*f,*Gtmp;
    REAL *amp2,*phase,*Gamp2,*Gphase;
    cufftHandle fftPlan;
    
    //welcome info
    
    printf("template program using ");
#ifdef DOUBLE
    printf("double");
#else
    printf("single");
#endif
    printf(" precision arithmetics.\n");
    
    //default parameters
    //assume square
    dim=2;
    N=256;
    L=256.0;

    // check if arguments are present and read them
    
    if(argc > 1 ) N = atoi(argv[1]);
    
    //excute
    

    cudaSetDevice(0);
    
    cudaMemGetInfo(&fmem,&tmem);
    printf("GPU memory before allocation free: %lu, total: %lu\n",fmem,tmem);

    threads=MAXT;
    blocks=N*N/threads+(N*N%threads==0?0:1);
    
    
    cudaMalloc(&GF,N*N*sizeof(COMPLEX));
    cudaMalloc(&Gtmp,N*N*sizeof(COMPLEX));
    cudaMalloc(&Gamp2,N*N*sizeof(REAL));
    cudaMalloc(&Gphase,N*N*sizeof(REAL));
    f=new COMPLEX[N*N];
    amp2=new REAL[N*N];
    phase=new REAL[N*N];
    
    
    //for FFT
    dkx=TWOPI/L;
    
    //include normalization in dkx:
    dkx=dkx/(1.0*N*N);
    //we need a "plan"
#ifdef DOUBLE
         if(dim==1) cufftPlan1d(&fftPlan, N, CUFFT_Z2Z,1);
    else if(dim==2) cufftPlan2d(&fftPlan, N, N, CUFFT_Z2Z) ;
    else if(dim==3) cufftPlan3d(&fftPlan, N, N, N, CUFFT_Z2Z);
#else
         if(dim==1) cufftPlan1d(&fftPlan, N, CUFFT_C2C,1);
    else if(dim==2) cufftPlan2d(&fftPlan, N, N, CUFFT_C2C) ;
    else if(dim==3) cufftPlan3d(&fftPlan, N, N, N, CUFFT_C2C);
#endif
    cerr("FFT plan"); //check for error
    
    
    t=0.0;dt=0.1;
    for(n=0;n<100;n++) {
        G_function<<<blocks,threads>>>(N,N,GF,t);
        
        
        //output
        G_ampphase<<<blocks,threads>>>(N*N,GF,Gamp2,Gphase);
        cudaMemcpy(amp2,Gamp2,N*N*sizeof(REAL),cudaMemcpyDeviceToHost);
        cudaMemcpy(phase,Gphase,N*N*sizeof(REAL),cudaMemcpyDeviceToHost);
        writeBM_real("test_amp2_"+IntToStrF(n,4),amp2,0,2,N,N,5);
        writeBM_real("test_phase_"+IntToStrF(n,4),phase,-PI,PI,N,N,6);
        
        //FFT
        G_FFT(GF,Gtmp,fftPlan); //forward
        G_kernel<<<blocks,threads>>>(N,Gtmp,dkx);
        G_FFT(Gtmp,GF,fftPlan,false); //inverse
        
        //output
        G_splitreim<<<blocks,threads>>>(N*N,GF,Gamp2,Gphase);
        cudaMemcpy(amp2,Gamp2,N*N*sizeof(REAL),cudaMemcpyDeviceToHost);
        cudaMemcpy(phase,Gphase,N*N*sizeof(REAL),cudaMemcpyDeviceToHost);
        mval=Mval=amp2[0];
        for(i=1;i<N*N;i++) {x=amp2[i];if(x>Mval) Mval=x;else if(x<mval) mval=x;}
        printf("%f %f; ",mval,Mval);
        writeBM_real("test_dx_re_"+IntToStrF(n,4),amp2,mval,Mval,N,N,5);
        mval=Mval=phase[0];
        for(i=1;i<N*N;i++) {x=phase[i];if(x>Mval) Mval=x;else if(x<mval) mval=x;}
        printf("%f %f\n",mval,Mval);
        writeBM_real("test_dx_im_"+IntToStrF(n,4),phase,mval,Mval,N,N,5);
        
        
        t+=dt;
    }
    
    delete[] f;
    delete[] amp2;
    delete[] phase;
    cudaFree(GF);
    cudaFree(Gtmp);
    cudaFree(Gamp2);
    cudaFree(Gphase);

    
    return 0;
    }

// ******************************************************
