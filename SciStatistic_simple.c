#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <omp.h>
#include "tempDataAnalysisFunctions.h"
//#include <fftw3.h>
//#include "DataAnalysisFunctions.h"


//#define NUM_THREADS 4

int main();
void statistic(double *Nx, double *maxvx, double *minvx, double *meanvx, double *variancex,double *devpadx, double *meanqx,double *Ny, double *maxvy, double *minvy, double *meanvy, double *variancey, double *devpady, double *meanqy);

void readstats(double *Nx, double *maxvx, double *minvx, double *meanvx, double *variancex, double *devpadx, double *meanqx,double *Ny, double *maxvy, double *minvy, double *meanvy,double *variancey, double *devpady, double *meanqy);

void tratamento_vx_spat_file(double *N, double *meanv,double *devpad);
void tratamento_vy_spat_file(double *N, double *meanv,double *devpad);
void treatement_file(double *meanv);
void tratamento_vx_spat(double *data,double *N, double *meanv,double *devpad);
void tratamento_vy_spat(double *data,double *N, double *meanv,double *devpad);
void get_treated_data(double *data, char *dir);
void structure_23456(double *data,int size);
void evaluate_corr_xdir(FILE *fout_c, FILE *fout_spec,double *N, int size,double *data);
void evaluate_corr_ydir(FILE *fout_c, FILE *fout_spec,double *N, int size,double *data);

//---------------------------------------------------

int main ()
{

	//Statistics variables

	int i,flag=1,j;
	double l1,l2;
	double Nx,Ny;
	double aux;
	double maxvx,minvx,meanvx,variancex,devpadx,meanqx;
	double maxvy,minvy,meanvy,variancey,devpady,meanqy;
	double *datax = NULL,*datay = NULL;
	double t1,t2;
	char c;
	double hflat;

	//Correlation and spectrum variables

	FILE *fout_corr_a11 = NULL;
	FILE *fout_corr_a12 = NULL;
	FILE *fout_corr_a21 = NULL;
	FILE *fout_corr_a22 = NULL;
	FILE *fout_spec_a11 = NULL;
	FILE *fout_spec_a12 = NULL;
	FILE *fout_spec_a21 = NULL;
	FILE *fout_spec_a22 = NULL;


	fout_corr_a11 = fopen("../out_statistic/correlations/corrA11.dat","w");
	fout_corr_a12 = fopen("../out_statistic/correlations/corrA12.dat","w");
	fout_corr_a21 = fopen("../out_statistic/correlations/corrA21.dat","w");
	fout_corr_a22 = fopen("../out_statistic/correlations/corrA22.dat","w");

	fout_spec_a11 = fopen("../out_statistic/correlations/specA11.dat","w");
	fout_spec_a12 = fopen("../out_statistic/correlations/specA12.dat","w");
	fout_spec_a21 = fopen("../out_statistic/correlations/specA21.dat","w");
	fout_spec_a22 = fopen("../out_statistic/correlations/specA22.dat","w");

	//Histogram variables
	FILE *out_hist = NULL;
	int Nbins,displ;
	double *rearr_y = NULL, *rearr_x = NULL;
	char *fileout = NULL,str[4];
	fileout = malloc(43*sizeof(char));

	FILE *out_flat = NULL;
	FILE *out_kurt = NULL;
	out_flat = fopen("./out_statistic/hflat_bypdf_long.dat","a");
	out_kurt = fopen("./out_statistic/kurt_bypdf_long.dat","a");

	FILE *out_struc = NULL;
	out_struc = fopen("./out_statistic/StructureFunctions_transv.dat","a");

	//Structure variables	
	char fx[] = "../get/data/vx_long_isot_t25to200_.txt";
//	char fy[] = "./get/data/vy_isot_t60_N1024_completing.txt";
	char ftreatedx[] = "../out_statistic/treated_data/vy_transv_isot_t25to200.txt";
//	char ftreatedy[] = "../out_statistic/treated_data/vy_long_isot_t25to200.txt";
	
	FILE *fdatax = NULL;
//	FILE *fdatay = NULL;
	fdatax = fopen(ftreatedx,"r");
//	fdatay = fopen(ftreatedy,"r");

		t1=omp_get_wtime();
		printf("Starting statistics...\n");
		statistic(&Nx,&maxvx,&minvx,&meanvx,&variancex,&devpadx,&meanqx,&Ny,&maxvy,&minvy,&meanvy,&variancey,&devpady,&meanqy);
		t2=omp_get_wtime();
		printf("Statistic time = %lf s\n",t2-t1);



}

//O statistics serve para obter a m�dia, o num de entrada, o maximo e o minimo da velocidade, a variancia e o desvio padr�o
//--------------------------------------------
void statistic(double *Nx, double *maxvx, double *minvx, double *meanvx, double *variancex, 
double *devpadx, double *meanqx,double *Ny, double *maxvy, double *minvy, double *meanvy, 
double *variancey, double *devpady, double *meanqy){
	
	FILE *fleitura_x = NULL;
	FILE *fleitura_y = NULL;
	FILE *fstat = NULL;
	FILE *statfile = NULL;
//	double temp;
	float temp;

//	fleitura_x=fopen("./get/data/vx_isot_t60_N1024_completing.txt","r");
//	fleitura_y=fopen("./get/data/vy_isot_t60_N1024_completing.txt","r");
	fleitura_x = fopen("../get/data/vy_transv_isot_t25to200_.txt","r");
	fstat=fopen("../out_statistic/GeneralStats_transv.dat","w");
	statfile=fopen("../out_statistic/ReadyStats_transv.dat","w");

//	while(fscanf(fleitura_x,"%lf\n",&temp) != EOF){
	while(fscanf(fleitura_x,"%f\n",&temp) != EOF){
		*Nx=*Nx+1;
		if(temp>*maxvx) {*maxvx=temp;}
		if(*Nx==1) {*minvx=*maxvx;}		
		if(temp<*minvx) {*minvx=temp;}
		*meanvx+=temp;
		}

	*meanvx=(*meanvx)/(*Nx);	
	rewind(fleitura_x);

//	while(fscanf(fleitura_x,"%lf\n",&temp) != EOF){
	while(fscanf(fleitura_x,"%f\n",&temp) != EOF){
		*variancex+=(temp-*meanvx)*(temp-*meanvx)/(*Nx-1);
		*meanqx+=(temp-*meanvx)*(temp-*meanvx);
		}	

	*meanqx=(*meanqx)/(*Nx);
	*devpadx=sqrt(*variancex);

	fprintf(fstat,"Para a dire��o x,\nN=%2.0lf\nmax=%lf\nmin=%lf\nmean=%lf\nvariance=%lf\ndesvio padrao=%lf\nmeanq=%lf\n",*Nx,*maxvx,
*minvx,*meanvx,*variancex,*devpadx,*meanqx);

	fprintf(statfile,"%2.0lf %lf %lf %lf %lf %lf %lf ",*Nx,*maxvx,*minvx,*meanvx,*variancex,*devpadx,*meanqx);
/*
	while(fscanf(fleitura_y,"%lf\n",&temp) != EOF){
		*Ny=*Ny+1;
		if(temp>*maxvy) {*maxvy=temp;}
		if(*Ny==1) {*minvy=*maxvy;}		
		if(temp<*minvy) {*minvy=temp;}
		*meanvy+=temp;
		}

	*meanvy=(*meanvy)/(*Ny);	
	rewind(fleitura_y);

	while(fscanf(fleitura_y,"%lf\n",&temp) != EOF){
		*variancey+=(temp-*meanvy)*(temp-*meanvy)/(*Ny-1);
		*meanqy+=(temp-*meanvy)*(temp-*meanvy);
		}	

	*meanqy=(*meanqy)/(*Ny);
	*devpady=sqrt(*variancey);

	fprintf(fstat,"\n\nPara a dire��o y,\nN=%2.0lf\nmax=%lf\nmin=%lf\nmean=%lf\nvariance=%lf\ndesvio padrao=%lf\nmeanq=%lf\n",*Ny,*maxvy,*minvy,*meanvy,*variancey,*devpady,*meanqy);

	fprintf(statfile,"%2.0lf %lf %lf %lf %lf %lf %lf\n",*Ny,*maxvy,*minvy,*meanvy,*variancey,*devpady,*meanqy);
*/
	rewind(fleitura_x);
	fclose(fleitura_x);
//	fclose(fleitura_y);
	fclose(fstat);
	fclose(statfile);
	}

//-------------------------------------------

void readstats(double *Nx, double *maxvx, double *minvx, double *meanvx, double *variancex, double *devpadx, double *meanqx,double *Ny, double *maxvy, double *minvy, double *meanvy, double *variancey, double *devpady, double *meanqy){
	
	FILE *stats = NULL;
	stats=fopen("./out_statistic/ReadyStats.dat","r");
	//double *t = NULL;
	//t = malloc(14*sizeof(double));

	double t[14];

	if(stats = NULL) {printf("Deu pau na leitura\n");}


	while(fscanf(stats,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&t[0],&t[1],&t[2],&t[3],&t[4],&t[5],&t[6],&t[7],&t[8],&t[9],&t[10],&t[11],&t[12],&t[13]) != EOF){
		printf("Alo\n");
		printf("t[0]=%lf\n",t[0]);
		*Nx = t[0];
		*maxvx = t[1];
		*minvx = t[2];
		*meanvx = t[3];
		*variancex = t[4];
		*devpadx = t[5];
		*meanqx = t[6];
		*Ny = t[7];
		*maxvy = t[8];
		*minvy = t[9];
		*meanvy = t[10];
		*variancey = t[11];
		*devpady = t[12];
		*meanqy = t[13];

		}
//	fscanf(stats,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Nx, &maxvx, &minvx, &meanvx, &variancex, &devpadx, &meanqx, &Ny, &maxvy, &minvy, &meanvy, &variancey, &devpady, &meanqy);



	//free(t);
	fclose(stats);


}
//--------------------------------------- FUNCIONANDO

//O tratamento retira a m�dia espacial de cada amostra
//------------------------------------------
void tratamento_vx_spat(double *data,double *N, double *meanv,double *devpad){
	int i;
	double temp;	
	FILE *fleitura = NULL;
	FILE *record_data = NULL;
	fleitura=fopen("../get/data/vx_long_isot_t25to200_.txt","r");
	record_data = fopen("../out_statistic/treated_data/vx_long_isot_t25to200.txt","w");	

	i=0;
	while(fscanf(fleitura,"%lf",&temp)!=EOF){
		data[i]=temp-*meanv;
		fprintf(record_data,"%f\n",(float)data[i]);
		i++;		
		}
	fclose(fleitura);
	fclose(record_data);
}

void tratamento_vy_spat(double *data,double *N, double *meanv,double *devpad){
	int i;
	double temp;	
	FILE *fleitura = NULL;
	FILE *record_data = NULL;

	fleitura=fopen("./get/data/vy_isot_t60_N1024_completing.txt","r");
	record_data = fopen("./out_statistic/treated_data/vy_isot_t60_N1024_completing.txt","w");	
	
	i=0;
	while(fscanf(fleitura,"%lf",&temp)!=EOF){
		data[i]=temp-*meanv;
		fprintf(record_data,"%f\n",(float)data[i]);
		i++;
		}
	fclose(fleitura);
	fclose(record_data);
}
//----------------------------------------
void tratamento_vx_spat_file(double *N, double *meanv,double *devpad){

	int i;
//	double temp;	
	float temp;	
	FILE *fleitura = NULL;
	FILE *record_data = NULL;
//	fleitura=fopen("./get/data/vx_isot_t60_N1024_completing.txt","r");
	fleitura=fopen("../get/data/vx_long_isot_t25to200_.txt","r");
	record_data = fopen("../out_statistic/treated_data/vx_long_isot_t25to200.txt","w");	

	i=0;
//	while(fscanf(fleitura,"%lf",&temp)!=EOF){
	while(fscanf(fleitura,"%f",&temp)!=EOF){
		fprintf(record_data,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura);
	fclose(record_data);
}

void tratamento_vy_spat_file(double *N, double *meanv,double *devpad){
	int i;
	double temp;	
	FILE *fleitura = NULL;
	FILE *record_data = NULL;

	fleitura=fopen("./get/data/vy_isot_t60_N1024_completing.txt","r");
	record_data = fopen("./out_statistic/treated_data/vy_isot_t60_N1024_completing.txt","w");	
	
	i=0;
	while(fscanf(fleitura,"%lf",&temp)!=EOF){
		fprintf(record_data,"%f\n",(float)(temp-*meanv));
		i++;
		}
	fclose(fleitura);
	fclose(record_data);
}


//-----------------------------------------
/*void get_treated_data(double *data, char *dir){
	int i;
	double temp;	
	FILE *fleitura = NULL;
	fleitura=fopen(dir,"r");

	i=0;
	while(fscanf(fleitura,"%lf",&temp)!=EOF){
		data[i]=temp;
		i++;		
		}
	fclose(fleitura);
}
*/

//--------------------------------------- FUNCIONANDO
/*

void evaluate_corr_xdir(FILE *fout_c, FILE *fout_spec,double *N, int size,double *data){
	double *input = NULL;
	double *output = NULL;
	double *serie = NULL;
	double *serief = NULL;          //Fourier transform of serie
	int n,max,i,j;

	n=size;

	input = malloc( n*sizeof(double));
	output = malloc( n*sizeof(double));
	serie = malloc( n*sizeof(double));	
	serief = malloc( n*sizeof(double));	
	
	max=(int)(*N/size);
	
	for(i=0;i<n;i++) {serie[i]=0.;}

	for(j=0;j<max;j++){
		for(i=0;i<n;i++){
			input[i]=data[j*n+i];
			output[i]=0.;
			}
		autocorr_realdata(input,output,n);
		for(i=0;i<n;i++) {serie[i]+=output[i]/((double)n);}
	}

	//Take the spectrum
	fourier_spectrum(serie,serief,n);

	for(i=0;i<n;i++){
		fprintf(fout_c,"%d\t%lf\n",i,serie[i]/serie[0]);
		if(serief[i]<pow(10,13)){fprintf(fout_spec,"%d\t%lf\n",i,serief[i]);}
		}

	free(input);
	free(output);
	free(serie);
	//free(serief);
}
*/
//------------------------------------------------
/*
void evaluate_corr_ydir(FILE *fout_c, FILE *fout_spec, double *N, int size,double *data){
	double *input = NULL;
	double *output = NULL;
	double *serie = NULL;
	double *serief = NULL;          //Fourier transform of serie
	int n,max,i,j;

	n=size;

	input = malloc( n*sizeof(double));
	output = malloc( n*sizeof(double));
	serie = malloc( n*sizeof(double));
	serief = malloc( n*sizeof(double));	
	
	max=(int)(*N/n);
	
	for(i=0;i<n;i++) {serie[i]=0.;}

	for(i=0;i<max;i++) {
		for(j=0;j<n;j++){
			input[j]=data[j*n+i];
			output[j]=0.;
			}
		autocorr_realdata(input,output,n);
		for(j=0;j<n;j++) {serie[j]+=output[j]/((double)max);}
	}

	//Take the spectrum
	fourier_spectrum(serie,serief,n);

	for(i=0;i<n;i++){
		fprintf(fout_c,"%d\t%lf\n",i,serie[i]/serie[0]);
		if(serief[i]<pow(10,13)){fprintf(fout_spec,"%d\t%lf\n",i,serief[i]/(2*pow(M_PI,3)));}
		}

	free(input);
	free(output);
	free(serie);
	free(serief);
}
*/
//------------------------------------------------
