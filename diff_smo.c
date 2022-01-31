//////////////       Smoluchowski equation discretization      /////////
///////////////// Required files: cmass.dat and density.xvg ////////////
/////////////// Arguments:                                     /////////
/////////////// 1) File sprob_perp.dat                         /////////
/////////////// 2) hmin_perp                                   /////////
/////////////// 3) hmax_perp                                   /////////
/////////////// 4) File density.xvg                            /////////

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>


#define BUFFER_SIZE 512 
#define MAX 100000

float trap(int i, float dr, int slab,float **x);

////////////////////////////////////////////////////////////////////////
/////////////// Main program ///////////////////////////////////////////

int main(int argc, char *argv[]){


/////////////// Global variables ///////////////////////////////////////

	int   i,j,k,m,atom;
	int   conta,count1,passo;
	int   tlim,slab;
	float tempops,dt,dt2,L,dr;
	float z[MAX],rho1[MAX],rho2[MAX],rho[MAX];
	float con,con1,ta0,ta1,ta2,r;
	float sumd,a,b,c,d,tempo;
	float teta0,teta1,teta2;
	float mean,mean1,fun,fun1;
	float *pt2,*surP,*surP1;
        float **pt,**pt1;
	float *mn,*mn1;

/////////////// System data ////////////////////////////////////////////

        int   totsteps = 5000000;
        int   interval = 100;
        int   steps    = totsteps/interval;
	int   nsteps   = 28*steps;
	

/////////////// Opening files //////////////////////////////////////////

	FILE *in,*in2,*in3,*in5,*in10;
	in  = fopen(argv[1],"r");
	in2 = fopen(argv[4],"r");
	
	in3  = fopen("results/sprob_perp_smo.dat","w");
	in5  = fopen("results/sprob_perp_ps.dat","w");
	in10 = fopen("results/diff_smo.dat","a");

////////////////////////////////////////////////////////////////////////	
/////////////// Perpendicular self-diffusion coefficient ///////////////
/////////////// Franco et al. 2016 /////////////////////////////////////
////////////////////////////////////////////////////////////////////////

	float hmin_perp = atof(argv[2]);
        float hmax_perp = atof(argv[3]);

	L = hmax_perp - hmin_perp;

	pt2   = (float *)malloc(steps * sizeof (float));

	tlim = 350;

        dt2=0.002*interval;
	
	for(i=0;i<tlim;i++){
                fscanf (in,"%d %f",&m, &pt2[i]);
                tempops = i*dt2;
                fprintf(in5,"%f %f\n",tempops,pt2[i]);
        }


/////////////// Derivative of the density with position (w) ////////////
/////////////// Eq. 18 Franco et al. 2016 //////////////////////////////
/////////////// Density file must have 2 columns ///////////////////////	

	slab = 990;// from gmx density
	
	char  string[BUFFER_SIZE];


	for(i=0;i<24;i++){
		fgets(string,BUFFER_SIZE,in2);
	}

	for(j=0;j<slab;j++){
		fscanf(in2,"%f %f %f\n",&z[j],&rho1[j],&rho2[j]);
		rho[j] = rho1[j];
	}
	
	dr = (z[slab-1]-z[0])/slab*1e-09;// space in meters

////////////////////////////////////////////////////////////////////////        
/////////////// Perpendicular self-diffusion coefficient ///////////////
/////////////// Smoluchowski discretization ////////////////////////////
//////////////  Franco et al. 2016 - Eq.67 /////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////// Boundary conditions /////////////////////////////////
/////////////////Eq. 74 Franco et al. 2016//////////////////////////////


	pt    = (float **)malloc(nsteps * sizeof (float *));
        pt1   = (float **)malloc(nsteps * sizeof (float *));

        for(i=0;i<nsteps;i++){
        pt[i] = (float *)malloc(slab * sizeof (float));
        pt1[i] = (float *)malloc(slab * sizeof (float));
                for (j=0;j<slab;j++){
                         pt[i][j] = 0.0;
                         pt1[i][j] = 0.0;
                         }
        }


///////////////////    Initial condition     //////////////////////////
//////////////////  Eq. 75 Franco et al. 2016 /////////////////////////
/////////////// Trapeizoidal Integration of density ///////////////////

        for(j=0;j<slab-1;j++){
                if(z[j]>=hmin_perp && z[j]<=hmax_perp){
			    sumd = sumd+rho[j]+rho[j+1];
		    }
        }

        sumd = sumd*dr/2.0;

	for (j=0;j<slab;j++){
                 if(z[j]>=hmin_perp && z[j]<=hmax_perp){
                        pt[0][j]=rho[j]/sumd;
                        pt1[0][j]=rho[j]/sumd;
                 }
        }


	dt = 1.5e-05*interval*1e-12;//time in seconds
        //Stability condition : check if not reached


////////////////// Golden Minimization search///////////////////////////

	r = (-1+sqrt(5))/2;

	a = 1e-10;
	b = 1e-07;

	while (b-a >1e-13){

	c = r*a+(1-r)*b;
	d = (1-r)*a+r*b;

	con = c*dt/dr/dr;
	con1 = d*dt/dr/dr;
	
	ta0 = ta1 = ta2 = 0.0;
    
	for (i=0;i<nsteps-1;i++){
		for (j=1;j<slab-1;j++){
			if(z[j]>=hmin_perp && z[j]<=hmax_perp){
	        		teta0 = 1.0-0.25*(log(rho[j+1]/rho[j-1]));
				//Eq. 68
				teta1 = -(2.0+dr/2.0*log(rho[j+1]*rho[j-1]/rho[j]/rho[j]));
				//Eq.69
				teta2 = 1.0+0.25*log(rho[j+1]/rho[j-1]);
				//Eq.70

				pt[i+1][j]     = pt[i][j] + con*(pt[i][j+1]*teta0+pt[i][j]*teta1+pt[i][j-1]*teta2);
				pt1[i+1][j]    = pt1[i][j]+ con1*(pt1[i][j+1]*teta0+pt1[i][j]*teta1+pt1[i][j-1]*teta2);
                    
				ta0 = ta0+teta0;
				ta1 = ta1+teta1;
				ta2 = ta2+teta2;

				}
			}
	}

// 	printf("Dt0=%.2e Dt1=%.2e Dt2=%.2e\n", con*ta0, con*ta1, con*ta2);
//	Stability conditions


//////////////////// Survival probability//////////////////////////////

	surP   = (float *)malloc(nsteps * sizeof (float));
        surP1  = (float *)malloc(nsteps * sizeof (float));
	
	for(i=0;i<nsteps;i++){
		surP[i]  = trap(i,dr,slab,pt);
		surP1[i] = trap(i,dr,slab,pt1);
		tempo    = i*dt*1.0e+12;//time in ps
	}

///////////// Diference between methods ////////////////////////////////

	mn   = (float *)malloc(nsteps * sizeof (float));
	mn1  = (float *)malloc(nsteps * sizeof (float));

	conta = 0;

	mean  = mean1 = 0.0;

	passo = nsteps/tlim;
	
	for (i=0;i<nsteps;i++){
		m = i + passo;
	      	for (j=i;j<m-1;j++){
			mean  = mean+surP[j];
			mean1 = mean1+surP1[j];
		}
		mn[conta]  = mean/passo;
		mn1[conta] = mean1/passo;

		conta = conta+1;
		i     = i+passo-1;
		mean  = 0.0;
		mean1 = 0.0;
	}

	free(surP);
	free(surP1);

/////////////////// Minimization function//////////////////////////////

	for (i=0;i<tlim;i++){
		fun  = fun+(mn[i]-pt2[i])*(mn[i]-pt2[i]);
		fun1 = fun1+(mn1[i]-pt2[i])*(mn1[i]-pt2[i]);
	}

	if (fun > fun1){
  		a = c ;
		b = b;
	}
    	else{
       		a = a;
		b = d;
	}

	count1 = count1+1;

}

	fprintf(in3,"0 1\n");

	for (i=1;i<conta;i++){
		tempops = i*dt2;
		fprintf(in3,"%d %f\n",i,mn[i]);
	}

	
	fprintf(in10,"%f %f %.2e\n",hmin_perp,hmax_perp,a+b/2);

	printf("Number of iterations: %d \n", count1);
	
	printf("%f %f %.2e\n",hmin_perp,hmax_perp,a+b/2);

free(pt);
free(pt1);
free(pt2);
free(mn);
free(mn1);

return 0;

}

/////////////// End of main program ////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/////////////// Trapezoidal rule for integration ///////////////////////

float trap(int i,float dr, int slab, float **x){

	int j;

	float sum = 0.0;

        for(j=0;j<slab-1;j++){
                sum = sum+x[i][j]+x[i][j+1];

        }

        sum = sum*dr/2.0;

return sum;

}



