#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <omp.h>
#include <fftw3.h>
#include "DataAnalysisFunctions.h"

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
	out_struc = fopen("./out_statistic/StructureFunctions_long.dat","a");

	//Structure variables	
	char fx[] = "../get/data/vx_long_isot_t25to200_.txt";
//	char fy[] = "./get/data/vy_isot_t60_N1024_completing.txt";
	char ftreatedx[] = "../out_statistic/treated_data/vx_long_isot_t25to200.txt";
//	char ftreatedy[] = "../out_statistic/treated_data/vy_long_isot_t25to200.txt";
	
	FILE *fdatax = NULL;
//	FILE *fdatay = NULL;
	fdatax = fopen(ftreatedx,"r");
//	fdatay = fopen(ftreatedy,"r");

	//Evaluating the main
/*
	printf("Did you have already evaluated the statistics? (y/n)\n");
	scanf("%c",&c);
	if(c == 'y') {
		printf("c=y\n");
		readstats(&Nx,&maxvx,&minvx,&meanvx,&variancex,&devpadx,&meanqx,&Ny,&maxvy,&minvy,&meanvy,&variancey,&devpady,&meanqy);}
	else {
		t1=omp_get_wtime();
		printf("Starting statistics...\n");
		statistic(&Nx,&maxvx,&minvx,&meanvx,&variancex,&devpadx,&meanqx,&Ny,&maxvy,&minvy,&meanvy,&variancey,&devpady,&meanqy);
		t2=omp_get_wtime();
		printf("Statistic time = %lf s\n",t2-t1);
		}
*/

		t1=omp_get_wtime();
		printf("Starting statistics...\n");
		statistic(&Nx,&maxvx,&minvx,&meanvx,&variancex,&devpadx,&meanqx,&Ny,&maxvy,&minvy,&meanvy,&variancey,&devpady,&meanqy);
		t2=omp_get_wtime();
		printf("Statistic time = %lf s\n",t2-t1);


	//printf("Alo\n");

	t1=omp_get_wtime();
//		datax = malloc(Nx*sizeof(double));
//		datay = malloc(Ny*sizeof(double));
	printf("Starting treatment\n");
//		get_treated_data(datax,ftreatedx);
		tratamento_vx_spat_file(&Nx,&meanvx,&devpadx);       // <----------- DESBLOQUAR AQUI SE FOR ARQUIVO NOVO
//		tratamento_vx_spat(datax,&Nx,&meanvx,&devpadx);
//		tratamento_vy_spat(datay,&Ny,&meanvy,&devpady);
	t2=omp_get_wtime();
	printf("Treatment time = %lf s\n",t2-t1);
//	t1=omp_get_wtime();
		printf("Doing the correlations...\n");
		evaluate_corr_xdir(fout_corr_a11,fout_spec_a11,&Nx,1024,datax);
//		evaluate_corr_xdir(fout_corr_a12,fout_spec_a12,&Ny,1024,datay);
//		evaluate_corr_ydir(fout_corr_a21,fout_spec_a21,&Ny,1024,datax);
//		evaluate_corr_ydir(fout_corr_a22,fout_spec_a22,&Ny,1024,datay);

	t2=omp_get_wtime();
	printf("Structure function time = %lf s\n",t2-t1);
*/	

	//PDFs in X-direction
	//Nbins = 1500;
	Nbins = 3000;
/*
	j = 0;
	printf("Longitudinal PDFs are being evaluated...\n");
	for(displ=80;displ>4;displ=displ/2){
//	for(displ=800;displ<1023;displ=displ+50){
		for(i=0;i<43;i++) {fileout[i]=0;}
		sprintf(str,"%d",displ);
		strcat(fileout,"./out_statistic/histograms/full_cube_xy/long_x_l");
		strcat(fileout,str);
		strcat(fileout,".dat");

		printf("%s\n",fileout);
	
		out_hist = fopen (fileout,"w");
//		histogram(out_hist,displ, datax,&Nx,1024,Nbins,20,&hflat);
		histogram_from_file(out_hist,displ,fdatax,&Nx,1024,Nbins,15,&hflat,out_kurt);
		fprintf(out_flat,"%d\t%lf\n",j,hflat);
		j++;
		}
	fclose(out_flat);

	j = 0;
	printf("Longitudinal PDFs are being evaluated...\n");
//	for(displ=93;displ>2;displ=displ/1.5){
	for(displ=90;displ>2;displ=displ/1.5){
//	for(displ=800;displ<1023;displ=displ+50){

		for(i=0;i<43;i++) {fileout[i]=0;}
		sprintf(str,"%d",displ);
		strcat(fileout,"./out_statistic/histograms/full_cube_xy/long_x_l");
		strcat(fileout,str);
		strcat(fileout,".dat");

		printf("%s\n",fileout);

		out_hist = fopen (fileout,"w");
//		histogram(out_hist,displ, datax,&Nx,1024,Nbins,20,&hflat);
		histogram_from_file(out_hist,displ,fdatax,&Nx,1024,Nbins,20,&hflat,out_kurt);
		fprintf(out_flat,"%d\t%lf\n",j,hflat);
		j++;

		}

	fclose(out_flat);
*/

/*
j = 0;
	printf("Transversal PDFs are being evaluated...\n");
	for(displ=80;displ>4;displ=displ/2){
//	for(displ=800;displ<1023;displ=displ+50){
		for(i=0;i<43;i++) {fileout[i]=0;}
		sprintf(str,"%d",displ);
		strcat(fileout,"./out_statistic/histograms/full_cube_xy/transv_l");
		strcat(fileout,str);
		strcat(fileout,".dat");

		printf("%s\n",fileout);
	
		out_hist = fopen (fileout,"w");
//		histogram(out_hist,displ, datax,&Nx,1024,Nbins,20,&hflat);
		histogram_from_file(out_hist,displ,fdatax,&Nx,1024,Nbins,15,&hflat,out_kurt);

		fprintf(out_flat,"%d\t%lf\n",j,hflat);
		j++;
		}

	j = 0;
	printf("Transversal PDFs are being evaluated...\n");
//	for(displ=93;displ>2;displ=displ/1.5){
	for(displ=90;displ>2;displ=displ/1.5){
//	for(displ=800;displ<1023;displ=displ+50){

		for(i=0;i<43;i++) {fileout[i]=0;}
		sprintf(str,"%d",displ);
		strcat(fileout,"./out_statistic/histograms/full_cube_xy/transv_l");
		strcat(fileout,str);
		strcat(fileout,".dat");

		printf("%s\n",fileout);

		out_hist = fopen (fileout,"w");
//		histogram(out_hist,displ, datax,&Nx,1024,Nbins,20,&hflat);

		histogram_from_file(out_hist,displ,fdatax,&Nx,1024,Nbins,20,&hflat,out_kurt);
		fprintf(out_flat,"%d\t%lf\n",j,hflat);
		j++;

		}
*/
	fclose(out_flat);

//	for(displ=5;displ<200;displ=displ*1.3){
////		structure_function(out_struc,datax, displ, 1024, Nx, 1,8);
//		structure_function_file(out_struc,fdatax, displ, 1024, Nx, 2,3);
//		printf("displ = %d done\n",displ);
//		}

	//Closing the files
	fcloseall();
	printf("Processo concluido com sucesso\n");

}

//O statistics serve para obter a média, o num de entrada, o maximo e o minimo da velocidade, a variancia e o desvio padrão
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
	fleitura_x = fopen("./get/data/vx+vy_long_isot_t60_N1024_full.txt","r");
	fstat=fopen("./out_statistic/GeneralStats_long.dat","w");
	statfile=fopen("./out_statistic/ReadyStats_long.dat","w");

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

	fprintf(fstat,"Para a direção x,\nN=%2.0lf\nmax=%lf\nmin=%lf\nmean=%lf\nvariance=%lf\ndesvio padrao=%lf\nmeanq=%lf\n",*Nx,*maxvx,
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

	fprintf(fstat,"\n\nPara a direção y,\nN=%2.0lf\nmax=%lf\nmin=%lf\nmean=%lf\nvariance=%lf\ndesvio padrao=%lf\nmeanq=%lf\n",*Ny,*maxvy,*minvy,*meanvy,*variancey,*devpady,*meanqy);

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

//O tratamento retira a média espacial de cada amostra
//------------------------------------------
void tratamento_vx_spat(double *data,double *N, double *meanv,double *devpad){
	int i;
	double temp;	
	FILE *fleitura = NULL;
	FILE *record_data = NULL;
	fleitura=fopen("./get/data/vx_isot_t60_N1024_completing.txt","r");
	record_data = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_completing.txt","w");	

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
	fleitura = fopen("./get/data/vx+vy_long_isot_t60_N1024_full.txt","r");
	record_data = fopen("./out_statistic/treated_data/vx+vy_long_isot_t60_N1024_full.txt","w");	

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
void treatement_file(double *meanv){

	int i;
//	double temp;	
	float temp;	
	FILE *fleitura1 = NULL;
	FILE *fleitura2 = NULL;
	FILE *fleitura3 = NULL;
	FILE *fleitura4 = NULL;
	FILE *fleitura5 = NULL;
	FILE *fleitura6 = NULL;

	FILE *record_data1 = NULL;
	FILE *record_data2 = NULL;
	FILE *record_data3 = NULL;
	FILE *record_data4 = NULL;
	FILE *record_data5 = NULL;
	FILE *record_data6 = NULL;
	
	fleitura1=fopen("./get/data/vx_isot_t60_N1024_133planes.txt","r");
	fleitura2=fopen("./get/data/vx_isot_t60_N1024_500to599.txt","r");
	fleitura3=fopen("./get/data/vx_isot_t60_N1024_600to699.txt","r");
	fleitura4=fopen("./get/data/vx_isot_t60_N1024_700to790.txt","r");
	fleitura5=fopen("./get/data/vx_isot_t60_N1024_791to899.txt","r");
	fleitura6=fopen("./get/data/vx_isot_t60_N1024_900to1023.txt","r");

	record_data1 = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_133planes.txt","w");	
	record_data2 = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_500to599.txt","w");	
	record_data3 = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_600to699.txt","w");	
	record_data4 = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_700to790.txt","w");	
	record_data5 = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_791to899.txt","w");	
	record_data6 = fopen("./out_statistic/treated_data/vx_isot_t60_N1024_900to1023.txt","w");	

	while(fscanf(fleitura1,"%f",&temp)!=EOF){
		fprintf(record_data1,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura1);
	fclose(record_data1);

	while(fscanf(fleitura2,"%f",&temp)!=EOF){
		fprintf(record_data2,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura2);
	fclose(record_data2);


	while(fscanf(fleitura3,"%f",&temp)!=EOF){
		fprintf(record_data3,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura3);
	fclose(record_data3);


	while(fscanf(fleitura4,"%f",&temp)!=EOF){
		fprintf(record_data4,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura4);
	fclose(record_data4);


	while(fscanf(fleitura5,"%f",&temp)!=EOF){
		fprintf(record_data5,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura5);
	fclose(record_data5);


	while(fscanf(fleitura6,"%f",&temp)!=EOF){
		fprintf(record_data6,"%f\n",(float)(temp-*meanv));
		i++;		
		}
	fclose(fleitura6);
	fclose(record_data6);




}

//-----------------------------------------
void get_treated_data(double *data, char *dir){
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


//--------------------------------------- FUNCIONANDO


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

//------------------------------------------------

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

//------------------------------------------------
