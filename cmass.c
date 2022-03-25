/////////////// Center of mass positions from GROMACS out.gro  /////////
///////gmx trjconv -f prod.trr -s prod.tpr -pbc nojump -o out.gro //////
/////////////// Arguments:                                     /////////
/////////////// 1) number of carbon atoms in the molecule      /////////
/////////////// 2) number of hydrogen atoms in the molecule    /////////
/////////////// 3) input file out.gro  			       /////////
/////////////// 4) output file cmass.dat		       /////////
/////////////// 5) nsteps in the simulation                    /////////
/////////////// 6) nstxout in the simulation                   ///////// 
/////////////// Spera and Braga - Jan, 2022 ////////////////////////////

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#define BUFFER_SIZE 256 
#define MAX 150000

int main(int argc, char *argv[]){

	int   i,j,k,aux,n1,n2;
	int   atom;
	float rx[MAX],ry[MAX],rz[MAX];//positions 
	float rx_cm[MAX],ry_cm[MAX],rz_cm[MAX];//center of mass position
	float vx[MAX],vy[MAX],vz[MAX];// velocities
	char  string[BUFFER_SIZE];
	char  *string1 = malloc(5);
	char  *string2 = malloc(5);
	float boxx,boxy,boxz;
	
	int   totsteps = atoi(argv[5]);
	int   interval = atoi(argv[6]); 
	int   steps    = totsteps/interval;

	int   c_number     = atoi(argv[1]);
	int   h_number     = atoi(argv[2]);
	int   atom_number  = c_number+h_number;

	int c_mass = 12.011;
	int h_mass = 1.008;

	FILE *in,*in2,*in3;
	in  = fopen(argv[3],"r");
	in2 = fopen(argv[4],"w");
//	in3 = fopen("teste.dat","w");

	for(k=0;k<=steps;k++){
	
		fgets(string,BUFFER_SIZE,in);

          	fscanf(in,"%5d",&n1);//number of particles

		aux = n1/atom_number;//number of molecules

		fprintf(in2,"%d\n",aux);

		for(j=0;j<aux;j++){
			
			for(i=0;i<atom_number;i++){
			
				fscanf(in,"%5d%5s%5s%5d%8f%8f%8f%8f%8f%8f\n",&n2,string1,string2,&atom,&rx[i],&ry[i],&rz[i],&vx[i],&vy[i],&vz[i]);
				// Read particles information for i lines
				if (strcmp(string2,"C") == 0){ //If string2=C 
					rx_cm[j] = rx_cm[j] + rx[i]*c_mass;//Sum of radius of all particles in a molecule
					ry_cm[j] = ry_cm[j] + ry[i]*c_mass;
					rz_cm[j] = rz_cm[j] + rz[i]*c_mass;
				//	fprintf(in3,"%5d %5s %8f\n",i,string2,rx[i]*12.011);
				}

				else{
					rx_cm[j] = rx_cm[j] + rx[i]*h_mass;//Sum of radius of all particles in a molecule
	                                ry_cm[j] = ry_cm[j] + ry[i]*h_mass;
        	                        rz_cm[j] = rz_cm[j] + rz[i]*h_mass;
				//	fprintf(in3,"%5d %5s %8f\n",i,string2,rx[i]*1.008);
				}

			
			}
			rx_cm[j]=rx_cm[j]/(c_number*c_mass+h_number*h_mass);//Divide total sum by number of molecules
			ry_cm[j]=ry_cm[j]/(c_number*c_mass+h_number*h_mass);
			rz_cm[j]=rz_cm[j]/(c_number*c_mass+h_number*h_mass);
		
			fprintf(in2,"%3d %8.3f %8.3f %8.3f\n",j,rx_cm[j],ry_cm[j],rz_cm[j]);// Write center of mass in new file

			rx_cm[j]=0.0;
                        ry_cm[j]=0.0;
                        rz_cm[j]=0.0;

		}
		
			fscanf(in,"%10f %10f %10f\n",&boxx,&boxy,&boxz);//Read size of the box 
	}

	free(string1);
	free(string2);
	fclose(in);
	fclose(in2);

return 0;

}

