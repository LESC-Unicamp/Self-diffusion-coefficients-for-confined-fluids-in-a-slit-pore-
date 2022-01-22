/////////////// Self-diffusion tensor from the survival probability ////
/////////////// of the center of masses for confined fluids ////////////
/////////////// Required files: cmass.dat and density.xvg //////////////
/////////////// Arguments:                                     /////////
/////////////// 1) File cmass.dat                              /////////
/////////////// 2) hmin_par                                    /////////
/////////////// 3) hmax_par                                    /////////
/////////////// 4) hmax_perp                                   /////////
/////////////// 5) hmax_perp                                   /////////
/////////////// 6) File density.xvg                            /////////
/////////////// 7) Lower limit for regression (msd_over_sprob) /////////
/////////////// 8) Upper limit for regression (msd_over_sprob) /////////
/////////////// Spera - April, 2021 ////////////////////////////////////

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#define BUFFER_SIZE 512 
#define MAX 100000

float mmq(int t,float x[MAX],float y[MAX],float low, float high);
float trap(float dt, int tlim,float x[MAX]);

////////////////////////////////////////////////////////////////////////
/////////////// Main program ///////////////////////////////////////////

int main(int argc, char *argv[]){

/////////////// Global variables ///////////////////////////////////////

	int   i,j,k,l,m,n1,n2,t;
	float aux,aux0,aux1,aux2;
	int   atom;
	float rx[MAX],ry[MAX],rz[MAX];
	float vx[MAX],vy[MAX],vz[MAX];
	float *sum_x,*sum_y;
	float **matx,**maty,**matz;
	int   **tag;
	float *n;
	float sumx[MAX],sumy[MAX];
	float *pt,*n0;
	float *msd_x,*msd_y;

/////////////// System data ////////////////////////////////////////////

	int   totsteps = 5000000;
	int   interval = 100;
	int   steps    = totsteps/interval;
	float hmin_par = atof(argv[2]);
	float hmax_par = atof(argv[3]);

/////////////// Opening files //////////////////////////////////////////

	FILE *in,*in2,*in3,*in4,*in5,*in6,*in7,*in8,*in9;
	in  = fopen(argv[1],"r");
	in2 = fopen("results/msd.dat","w");	
	in3 = fopen("results/sprob_par.dat","w");	
	in4 = fopen("results/msd_over_sprob.dat","w");	
	in5 = fopen("results/sprob_perp.dat","w");	
	in6 = fopen(argv[6],"r");	
	in7 = fopen("results/lndensity.xvg","w");	
	in8 = fopen("results/diff_confined.dat","a");	
	in9 = fopen("results/alpha.dat","a");	

/////////////// Memory allocation //////////////////////////////////////
	
	n     = (float *)malloc(steps * sizeof (float));
	sum_x = (float *)malloc(steps * sizeof (float));
	sum_y = (float *)malloc(steps * sizeof (float));
	msd_x = (float *)malloc(steps * sizeof (float));
	msd_y = (float *)malloc(steps * sizeof (float));
	pt    = (float *)malloc(steps * sizeof (float));
	n0    = (float *)malloc(steps * sizeof (float));
	tag   = (int **)malloc(steps * sizeof (int *));

	matx  = (float **)malloc(steps * sizeof (float *));
	maty  = (float **)malloc(steps * sizeof (float *));
	matz  = (float **)malloc(steps * sizeof (float *));

/////////////// Reading data file cmass.dat and creating matrix ////////

	for(j=0;j<steps;j++){
		fscanf(in,"%5d",&n1);
		aux = n1;
		matx[j] = (float *)malloc(aux * sizeof (float));
		maty[j] = (float *)malloc(aux * sizeof (float));
		matz[j] = (float *)malloc(aux * sizeof (float));
		for(i=0;i<aux;i++){
			fscanf(in,"%d %f %f %f",&k,&rx[i],&ry[i],&rz[i]);
			matx[j][i]=rx[i];
			maty[j][i]=ry[i];
			matz[j][i]=rz[i];
		}
	}

////////////////////////////////////////////////////////////////////////
/////////////// Parallel self-diffusion coefficients ///////////////////
/////////////// Liu et al. 2004 ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/////////////// MSD from Equation 14 from Liu et al. 2004 //////////////

	int tlim = 350;
	int jmax = steps-tlim;
	int tmax; 

	for(j=0;j<jmax;j++){
		tmax=j+tlim;		
		for(t=0;t<tmax;t++){
			sum_x[t]=0.0;
			sum_y[t]=0.0;
			n[t]    =0.0;
		}
		tag[j] = (int *)malloc(aux * sizeof (int));
		for(i=0;i<aux;i++){
			tag[j][i]=i;
		}
		for(t=j;t<tmax;t++){
			m=t-j;
			n[t] = 0.0;
			for(i=0;i<aux;i++){ 					 
				if(matz[t][i]>=hmin_par && matz[t][i]<=hmax_par && tag[j][i]!='\0'){
					sum_x[t] = sum_x[t]+(matx[t][i]-matx[j][i])*(matx[t][i]-matx[j][i]);
					sum_y[t] = sum_y[t]+(maty[t][i]-maty[j][i])*(maty[t][i]-maty[j][i]);
					n[t]=n[t]+1;
				}
				else{
					tag[j][i]='\0';
				}
			}
			msd_x[m] = msd_x[m] + sum_x[t]/n[j];
			msd_y[m] = msd_y[m] + sum_y[t]/n[j];
			pt[m]    = pt[m]    + n[t]/n[j];
			n0[m]    = n0[m]    + n[j];
		}
	}


/////////////// Averaging over multiple time origins /////////////////// 

	float x[MAX],y[MAX];

	for(m=0;m<tlim;m++){
		
/////////////// Mean square displacement ///////////////////////////////

		msd_x[m]=msd_x[m]/jmax;
		msd_y[m]=msd_y[m]/jmax;
		fprintf(in2,"%d %f %f\n",m,msd_x[m],msd_y[m]);

/////////////// Survival probability ///////////////////////////////////
		
		pt[m]=pt[m]/jmax;
		n0[m]=n0[m]/jmax;
		fprintf(in3,"%d %f %f\n",m,pt[m],n0[m]);
	
/////////////// Ratio - Eq. 16 Liu et al. 2004  ////////////////////////

		x[m]     = msd_x[m]/pt[m];
		y[m]     = msd_y[m]/pt[m];
		fprintf(in4,"%d %f %f\n",m,x[m],y[m]);
	}

//////////////// Parallel self-diffusion - Eq. 16 Liu et al. 2004 //////

	float tempo[MAX];
	float DXX,DYY;
	float dt=0.002*interval; //0.002ps*100(intervalo)
	
/////////////// Region of the data considered for analysis (linear) ////

	float low  = (atof(argv[7]))*dt;
	float high = (atof(argv[8]))*dt;
	
	for(i=0;i<tlim;i++){
       		tempo[i]=i*dt;
	}

	DXX=mmq(tlim,tempo,x,low,high);
	DXX=DXX/2.0/1e06;
 	
	DYY=mmq(tlim,tempo,y,low,high);
	DYY=DYY/2.0/1e06;

////////////////////////////////////////////////////////////////////////	
/////////////// Perpendicular self-diffusion coefficient ///////////////
/////////////// Franco et al. 2016 /////////////////////////////////////
////////////////////////////////////////////////////////////////////////

	float hmin_perp = atof(argv[4]);
	float hmax_perp = atof(argv[5]);

/////////////// Survival probability ///////////////////////////////////

	tlim = 350;
	jmax = steps-tlim;

	for(i=0;i<steps;i++){
		pt[i]=0.0;
	}

	for(j=0;j<jmax;j++){
		tmax=j+tlim;		
		for(t=0;t<tmax;t++){
			sum_x[t]=0.0;
			sum_y[t]=0.0;
			n[t]    =0.0;
		}
		tag[j] = (int *)malloc(aux * sizeof (int));
		for(i=0;i<aux;i++){
			tag[j][i]=i;
		}
		for(t=j;t<tmax;t++){
			m=t-j;
			n[t] = 0.0;
			for(i=0;i<aux;i++){ 					 
				if(matz[t][i]>=hmin_perp && matz[t][i]<=hmax_perp && tag[j][i]!='\0'){
					n[t]=n[t]+1;
				}
				else{
					tag[j][i]='\0';
				}
			}
			pt[m]    = pt[m]    + n[t]/n[j];
		}
	}


/////////////// Averaging over multiple time origins /////////////////// 

	for(m=0;m<tlim;m++){
		
		pt[m]=pt[m]/jmax;
		fprintf(in5,"%d %f\n",m,pt[m]);
	}

/////////////// Residence time by trapezoidal rule /////////////////////
/////////////// Eq. 40 Franco et al. 2016 //////////////////////////////

	float tau;
	
	tau=trap(dt,tlim,pt);
	printf("%d %f\n",tlim,tau);

/////////////// Layer width ////////////////////////////////////////////

	float L;

	L = hmax_perp - hmin_perp;

/////////////// Derivative of the density with position (w) ////////////
/////////////// Eq. 18 Franco et al. 2016 //////////////////////////////
/////////////// It requires one density file per component /////////////
/////////////// i.e., density file must have only 2 columns ////////////	

	int   slab = 1000;
	float w;
	float z[MAX],rho1[MAX],rho2[MAX],lnrho[MAX];

	char  string[BUFFER_SIZE];

//	for(i=0;i<24;i++){
//		fgets(string,BUFFER_SIZE,in6);
//	}

	for(i=0;i<slab;i++){
		fscanf(in6,"%f %f %f\n",&z[i],&rho1[i],&rho2[i]);
		lnrho[i]=log(rho1[i]);
		fprintf(in7,"%f %f\n",z[i],lnrho[i]);
	}

	w = mmq(slab,z,lnrho,hmin_perp,hmax_perp);

/////////////// Calculation of Alpha^-1 ////////////////////////////////
/////////////// Eq. 27 Franco et al. 2016 //////////////////////////////

	int   lim = 100;
	float sum = 0.0;
	
	float alpha1;

	for(j=0;j<lim;j++){
		aux  = 2.0*j + 1.0;
		aux0 = pow(aux,4)*pow(M_PI,4) + 0.75*pow(w,2)*pow(L,2)*pow(aux,2)*pow(M_PI,2) - 0.25*pow(w,4)*pow(L,4);
		sum = sum + 1/aux0;
	}

	aux1   = exp(w*L)+1.0;
	aux2   = exp(w*L)-1.0;

	alpha1 = 4*w*L*(aux1/aux2)*sum;
	fprintf(in9,"%f %f\n",(hmax_perp+hmin_perp)*0.5,(1.0/alpha1));

/////////////// Eq. 26 Franco et al. 2016 //////////////////////////////

	float DZZ;

	DZZ = L*L*alpha1/tau*1e-06;
	
	fprintf(in8,"%.2f %.2e %.2e %.2e\n",(hmax_par+hmin_par)*0.5,DXX,DYY,DZZ);
	
/////////////// Free the memory space allocated ////////////////////////

free(n);
free(matx);
free(maty);
free(matz);
free(msd_x);
free(msd_y);
free(sum_x);
free(sum_y);
free(pt);
free(tag);

return 0;

}

/////////////// End of main program ////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/////////////// Trapezoidal rule for integration ///////////////////////

float trap(float dt, int tlim,float x[MAX]){

	int   i;
	float sum = 0.0;

	for(i=0;i<tlim;i++){
		sum = sum+x[i]+x[i+1];
	}

	sum = sum*dt/2.0;

return sum;

}

/////////////// Linear regression function /////////////////////////////

float mmq(int t,float x[MAX],float y[MAX],float low, float high){

	int   i;
	int   count=0;
	float a    = 0.0;
	float b    = 0.0;
	float r    = 0.0;
	float r2   = 0.0;
	float sumx = 0.0;
	float sumy = 0.0;
	float sumxx= 0.0;
	float sumxy= 0.0;
	float sumyy= 0.0;
	float numx = 0.0;
	float denx = 0.0;
	float numr = 0.0;
	float denr = 0.0;

	for(i=0;i<t;i++){
		if(x[i]>=low && x[i]<high){
			sumx      = sumx+x[i];
			sumy      = sumy+y[i];
			sumxx     = sumxx+x[i]*x[i];
			sumxy	  = sumxy+x[i]*y[i];
			sumyy	  = sumyy+y[i]*y[i];
			count     = count+1;
		}
	}

	numx = (count*sumxy)-(sumx*sumy);
	denx = (count*sumxx)-(sumx*sumx);

	a = numx/denx;
	b = (sumy - a*sumx)/count;

	numr = (count*sumxy-sumx*sumy);
	denr = sqrt(abs((count*sumxx-sumx*sumx)*(count*sumyy-sumy*sumy)));
	
	r  = numr/denr;
	r2 = r*r;
	printf("%f\n",r2);
	
return a;

}


