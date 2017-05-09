#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

float f;
float C;
int N;
float tspan=1000;
float tstp = 500;
float vth = 0.3;
float vsl = 0.001;
int A = 3;
int B = 3;
float e = 0.01;
float c = -0.1;
float gsyn = 0.007;
float h = 0.001;

int* mat = NULL;
float* I=NULL;
float* E=NULL;
float** g = NULL;

void integrator_rk4(float dt,float t, float* p1,  float yout[N*(N+7)]);
void oscnetwork_opt(float t, float y[N*(N+7)],float dydt[N*(N+7)]);
int sign(float in);
float gate_rate(float Y,float alpha, float beta);

int main(int argc, char * argv[]){
	/* initializations*/

	N = atof(argv[1]);
	f = atof(argv[2]);
	C = atof(argv[3]);

	int i,j,iter;

//	initialize the timer
	srand(1);

//	printf("N*(N+7) :%d\n",N*(N+7));
	/* allocating space for all the global arrays of variable size*/
//	stimulation array
	if (I != 0) {
	    I = (float*) realloc(I, N * sizeof(float));
	} else {
	    I = (float*) malloc(N * sizeof(float));
	}
//	matrix giving the locations of voltage in the overall array	
	if (mat != 0) {
	    mat = (int*) realloc(mat, N * sizeof(int));
	} else {
	    mat = (int*) malloc(N * sizeof(int));
	}

// 	reversal potential array ( -ve  : inhibitory, +ve : excitatory)
	if (E != 0) {
	    E = (float*) realloc(E, N * sizeof(float));
	} else {
	    E = (float*) malloc(N * sizeof(float));
	}

//	weights for the adjacency matrix : g synaptic
	g = malloc(N * sizeof(float *));
	for(i = 0; i < N; i++){
		g[i] = malloc(N * sizeof(float));
		if(g[i] == NULL){
	        	fprintf(stderr, "out of memory\n");
        		break;
		}
	}
//	counter init		
	float count;
// 	fraction of neurons that are stimulated
	float stim= 1;
	

//	initialize a pointer		
	float *q;	

	float randomnum;
	
// 	total number of time iterations = tspan*step_size	
	int tot_time = (int) ceil(tspan/h);
	
//	Time array
	float T[tot_time];

// 	Time variable	
	float t;

//	Initialize the stimulation variable with magnitude
	float stimag = 1.0;	 
	for (i = 0;i<(int) floor(stim*N);i++){
		I[i]= stimag;	

	}
	
	for (i = (int) floor(stim*N);i<N;i++){
		I[i]= 0; 	

	}

	
//	randomize the stimulation array	
	for (i = 0; i < N; i++) 
        {
          j = i + rand() / (RAND_MAX / (N - i) + 1);
          t = I[j];
//	  printf(" j: I[%d] = %f  ",j,t);	
          I[j] = I[i];
          I[i] = t;
//	  printf("j: %d,I[%d]= %f ",j,i,I[i]);
        }
	printf("\n");

// 	initialize the reversal potential array
	for (i = 0;i<(int) floor(f*N);i++){
		E[i]= -5;	
	}	
	for (i =(int) floor(f*N);i<N;i++){
		E[i]= 5;	
	}

// 	randomize the reversal potential array
	for (i = 0; i < N; i++) 
        {
		j = i + rand() / (RAND_MAX / (N - i) + 1);
		t = E[j];
		E[j] = E[i];
		E[i] = t;
//		printf("E[i],%f ",E[i]);
        }

// 	initialize the adjacency matrix
	for (i=0;i<N;i++){	
		for(j=0;j<N;j++){
			randomnum = (float) rand()/(float) RAND_MAX;	
			if (randomnum<C){
				g[i][j] = gsyn;

			}
			else{ 	
				g[i][j] = 0;
			}
			if (i==j){
				g[i][j]  =0;
			} 	
		}
	}

// 	vector to hold values for each differential variable for all time iterations
	float Y[N*(N+7)];

//	initialize the vector of indices showing the locations of the voltage variables
//	The first location is 0 because the first variable in the differential vector will
//	always be v1, having the index 0.
	mat[0] = 0;	
//	printf("mat[0] = 0  ");
//	To get to the location of the next voltage variable, add 2 for vi and wi, 
//	add the sum of the adjacency matrix for the previous row	
//	printf("mat[0] = 0 ");
	for (i=1;i<N;i++){	
		mat[i] = mat[i-1]+N+7;	
	}
	printf("\n");
// 	initial conditions vector for time = 0
	for (i=0;i<N*(N+7);i++)
		Y[i] = 0;
	
//	set the time array
	T[0] = 0;
	
//	receive the result into this array
	float Y1[N*(N+7)];

// 	This loop calls the RK4 code
	for (i=0;i<tot_time-1;i++){
		q = Y;
//		call the RK4 integrator with current time value, and current 
//		values of voltage			
		integrator_rk4(h,T[i],Y,Y1);  
		
//		Return the time output of integrator into the next iteration of time
		T[i+1] = T[i]+h;	
		
//		copy the output of the integrator into the next iteration of voltage		
		q = memcpy(q, Y1, N*(N+7) * sizeof(float));
//		for (i=0;i<N*(N+7);i++)
//			Y[i] = t_y.y[i];

		printf("%f ",T[i]);
		count =0 ;
		for (iter = 0;iter<N;iter++){
			printf("%f ",Y1[mat[iter]]);
			if (Y1[mat[iter]] >-20)
				count+=1;			 
		}
		printf("%f ", count/100);	
		printf("\n");
	}

// 	free all the memory that was allocated for the arrays		
	free(E);
	free(I);

	for(i = 0; i < N; i++)
    		free(g[i]);
	free(g);
	free(mat);
	return 0;
}

float gate_rate(float Y,float alpha, float beta){
	float out;
	out = alpha*(1-Y) - beta*Y;
	return out;
}

int sign(float a){
	int out;
	if (a>0)
		out = 1;
	if (a == 0)
		out = 0;
	if (a < 0)
		out = -1;

	return out;
}

void integrator_rk4(float dt,float t, float y[N*(N+7)], float yout[N*(N+7)])
{	
//	initialize all the pointers
	float y1[N*(N+7)],y2[N*(N+7)],y3[N*(N+7)];
	float tout,dt_half;
	float k1[N*(N+7)],k2[N*(N+7)],k3[N*(N+7)],k4[N*(N+7)];
// 	initialize iterator
	int i;

	tout = t+dt;
	dt_half = 0.5*dt;
	float addition[N*(N+7)];

//	return the differential array into k1
	oscnetwork_opt(t,y,k1);

// 	multiply the array k1 by dt_half
	for(i=0;i<N*(N+7);i++)
		y1[i]=y[i]+(k1[i])*dt_half;	


// 	do the same thing 3 times
	oscnetwork_opt(t+dt_half,y1,k2);
	for(i=0;i<N*(N+7);i++)
		y2[i]=y[i]+(k2[i])*dt_half;	
	
	oscnetwork_opt(t+dt_half,y2,k3);
	for(i=0;i<N*(N+7);i++)
		y3[i]=y[i]+(k3[i])*dt;	

	oscnetwork_opt(tout,y3,k4);
//	Make the final additions with k1,k2,k3 and k4 according to the RK4 code
	for (i=0;i<N*(N+7);i++){
		addition[i] = ((k1[i]) + (k2[i])*2 + (k3[i])*2 + (k4[i])) *dt/6;
	}
//	add this to the original array
	for(i=0;i<N*(N+7);i++)
		yout[i]=y[i]+addition[i];		
}



// function to return the vector with coupled differential variables for each time iteration
void oscnetwork_opt(float t, float y[N*(N+7)],float dydt[N*(N+7)]){
	
	int i,j;

//	parameters for PR	
	float Vsyn,sum_syn_pr[N],theta,a_syn,b_syn;//,F_Pre;
	int count1_pr[N],count2[N];
	float g_c,p;
	g_c= 2.1;
	p  = 0.5;
//	float am_pr,bm_pr,an_pr,bn_pr,ah_pr,bh_pr,as_pr,bs_pr,aq_pr,bq_pr,ac_pr,bc_pr,m_inf_pr;
	float EL_pr,ENa_pr,EK_pr,ECa_pr;
	float Cm_pr;

	float gmax[6] = {0.1,30.0, 15.0,10.0,0.8,15.0};

//	float gNa_pr, gCa_pr, gDR_pr,gC_pr,gAHP_pr;
	EL_pr   = -60.0;
	ENa_pr  = 55.0;
	EK_pr   = -75.0;
	ECa_pr  = 80.0;  
	
//	parameters for WB
	float sum_syn_wb[N];
	int count1_wb[N];
	float gL_wb,gK_wb,gNa_wb;
	float EK_wb,EL_wb,ENa_wb,phi;
//	float am_wb,bm_wb,an_wb,bn_wb,ah_wb,bh_wb,m_inf_wb;

//	float IK_wb,IL_wb,INa_wb;
	float Cm_wb;

	gL_wb  = 0.1;
	gK_wb  = 9;
	gNa_wb = 35;

	EK_wb  = -90;
	EL_wb  = -65;
	ENa_wb =  55;

	phi = 5.0;

	theta=0;
	a_syn=12.0;
	b_syn=0.1;

	Vsyn =-75.0;
// 	i is reserved to iterate over the indices of voltage as mat[i].
//	i : 0...N-1
	#pragma omp parallel for private(j)
	for (i=0;i<N;i++){
//		count : 0... sum_column
		if (E[i]>0){
			sum_syn_pr[i]= 0;
			count1_pr[i] = 0;	
			for (j=0;j<N;j++){	
				if (j!=i){
					sum_syn_pr[i] += g[i][j]*y[i*(N+7)+8+count1_pr[i]]*(y[mat[j]]-Vsyn);
					count1_pr[i]+=1;
				}
			}
		// 	i is reserved to iterate over the indices of voltage as mat[i].
		//	i : 0...N-1
			Cm_pr  = 3.0;
			float am_pr = - (0.32*(y[i*(N+7)+0]+46.9))/(exp(-(y[i*(N+7)+0]+46.9)/4.0) - 1.0);
			float an_pr =  - (0.016*(y[i*(N+7)+0]+24.9))/(exp(-(y[i*(N+7)+0]+24.9)/5.0)-1.0);
			float ah_pr = 0.128 * exp((-43-y[i*(N+7)+0])/18.0);
			float as_pr = 1.6/(1.0 + exp(-(0.072*(y[i*(N+7)+1]-5.0))));
			float bm_pr  =  (0.28*(y[i*(N+7)+0]+19.9))/(exp((y[i*(N+7)+0]+19.9)/5.0) - 1.0);
			float bn_pr  = 0.25 * exp(-(y[i*(N+7)+0]+40.0)/40.0);
			float bh_pr  = 4.0/(1.0 + exp(-(y[i*(N+7)+0]+20)/5.0));
			float bs_pr  = (0.02*(y[i*(N+7)+1]+8.9))/(exp((y[i*(N+7)+1]+8.9)/5.0) - 1.0);
			float aq_pr,ac_pr,bc_pr,gC_pr; 
			if (0.00002*y[i*(N+7)+2]<0.01)
				aq_pr = 0.00002*y[i*(N+7)+2];
			else
				aq_pr =0.01;
			float bq_pr  = 0.001;
			if (y[i*(N+7)+1] <= -10.0){
				ac_pr = 0.0527*(exp( ((y[i*(N+7)+1]+50.0)/11.0) - ((y[i*(N+7)+1]+53.5)/27.0) ));        
				bc_pr  = 2.0 * exp(-(y[i*(N+7)+1]+53.5)/27.0) - ac_pr;
			}
			else{
				ac_pr = 2.0 * exp(-(y[i*(N+7)+1]+53.5)/27.0) ;
				bc_pr  = 0.0;
			}


			float m_inf_pr = am_pr/(am_pr+bm_pr);

			float gDR_pr = gmax[2]*y[i*(N+7)+3];
			float gNa_pr = gmax[1]*m_inf_pr*m_inf_pr*y[i*(N+7)+4];
			float gCa_pr = gmax[3]*y[i*(N+7)+5]*y[i*(N+7)+5];
			if (y[i*(N+7)+2]/250<1)
				gC_pr = gmax[5]*y[i*(N+7)+6]*y[i*(N+7)+2]/250.0;
			else
					gC_pr = gmax[5]*y[i*(N+7)+6];
			float gAHP_pr = gmax[4]*y[i*(N+7)+7]; 

			dydt[i*(N+7)+0] = (-gmax[0]*(y[i*(N+7)+0]-EL_pr)-gNa_pr*(y[i*(N+7)+0]-ENa_pr)-gDR_pr*(y[i*(N+7)+0]-EK_pr)+(g_c/p)*(y[i*(N+7)+1]-y[i*(N+7)+0])+((0.5+0.5*sign(t-tstp))*I[i])/p)/Cm_pr;
			dydt[i*(N+7)+1] = (-gmax[0]*(y[i*(N+7)+1]-EL_pr)-gCa_pr*(y[i*(N+7)+1]-ECa_pr)-gAHP_pr*(y[i*(N+7)+1]-EK_pr)-gC_pr*(y[i*(N+7)+1]-EK_pr)+(g_c/(1-p))*(y[i*(N+7)+0]-y[i*(N+7)+1])+ sum_syn_pr[i]/(1-p))/Cm_pr;
			dydt[i*(N+7)+2] =  -0.13*(gCa_pr*(y[i*(N+7)+1]-ECa_pr))-0.075*y[i*(N+7)+2];
			dydt[i*(N+7)+3] = gate_rate(y[i*(N+7)+3],an_pr,bn_pr);
			dydt[i*(N+7)+4] = gate_rate(y[i*(N+7)+4],ah_pr,bh_pr);
			dydt[i*(N+7)+5] = gate_rate(y[i*(N+7)+5],as_pr,bs_pr);
			dydt[i*(N+7)+6] = gate_rate(y[i*(N+7)+6],ac_pr,bc_pr);
			dydt[i*(N+7)+7] = gate_rate(y[i*(N+7)+7],aq_pr,bq_pr);

		}
		else if (E[i]<0){
			count1_wb[i]=0;
			sum_syn_wb[i]=0;
		//	initialize iterators
			for (j=0;j<N;j++){	
				if (j!=i){
					sum_syn_wb[i] += g[i][j]*y[i*(N+7)+8+count1_wb[i]]*(y[mat[j]]-Vsyn);
					count1_wb[i]+=1;
				}
			}	
			Cm_wb  = 1.5;
	
			float am_wb = -0.1*(y[i*(N+7)+0]+35)/(exp(-0.1*(y[i*(N+7)+0]+35))-1);
			float bm_wb  = 4*exp(-(y[i*(N+7)+0]+60)/18);
			float ah_wb = 0.07*exp(-(y[i*(N+7)+0]+58)/20);
			float bh_wb  = 1/(exp(-0.1*(y[i*(N+7)+0]+28))+1);
			float an_wb = -0.01*(y[i*(N+7)+0]+34)/(exp(-0.1*(y[i*(N+7)+0]+34))-1);
			float bn_wb  = 0.125*exp(-(y[i*(N+7)+0]+ 44)/80);

			float IK_wb = gK_wb*(y[i*(N+7)+1]*y[i*(N+7)+1]*y[i*(N+7)+1]*y[i*(N+7)+1])*(y[i*(N+7)+0]-EK_wb);
			float IL_wb = gL_wb*(y[i*(N+7)+0]-EL_wb);
			float m_inf_wb = am_wb/(am_wb+bm_wb);
			float INa_wb = gNa_wb*m_inf_wb*m_inf_wb*m_inf_wb*y[i*(N+7)+2]*(y[i*(N+7)+0]-ENa_wb);	
	
			dydt[i*(N+7)+0] = -(1/Cm_wb)*(INa_wb+IK_wb+IL_wb+sum_syn_wb[i]-((0.5+0.5*sign(t-tstp))*I[i]));	
			dydt[i*(N+7)+1] = phi*gate_rate(y[i*(N+7)+1],an_wb,bn_wb);
			dydt[i*(N+7)+2] = phi*gate_rate(y[i*(N+7)+2],ah_wb,bh_wb);
			dydt[i*(N+7)+3] = 0;
			dydt[i*(N+7)+4] = 0;
			dydt[i*(N+7)+5] = 0;
			dydt[i*(N+7)+6] = 0;
			dydt[i*(N+7)+7] = 0;
		}	
		count2[i] = 0;		

//		for each of the synaptic variables
		for (j = 0;j<N;j++){
			if (j!=i){	
				float F_Pre = 1.0/(1.0+exp(-(y[mat[j]]-theta)/2.0));
				dydt[i*(N+7)+8+count2[i]] = a_syn*F_Pre*(1.0-y[i*(N+7)+8+count2[i]])-b_syn*(y[i*(N+7)+8+count2[i]]); 
				count2[i] += 1;
			}			
		}
	
			
	}
}




	















	
	


 















