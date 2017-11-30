#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<time.h>

void rand_mat(double *mat, int dim)
{
  int i,j; 
  for(i=0; i<dim; i++) 
    for(j=0; j<dim; j++) 
      *(mat+(j+i*dim)) = (100.0*rand())/RAND_MAX;
}

void mat_mult_thr(double *mat1, double *mat2, double *result, int dim, int nthr, int choice)
{

  
  int part_rows, th_id;
  part_rows = dim/nthr;
  omp_set_num_threads(nthr); 
  #pragma omp parallel shared(mat1,mat2,result,part_rows) private(th_id)
  {
    int i,j,k; 
    th_id = omp_get_thread_num(); 

    #pragma omp for schedule(guided,part_rows)
    for(i=0; i<dim; i++) 
    {
	if(choice)
      		printf("Thread #%d is doing row %d.\n",th_id,i); 
      for(j=0; j<dim; j++) 
      {
        *( result+(j+i*dim) ) = 0; 
        for(k=0; k<dim; k++)
          *( result+(j+i*dim) ) += *( mat1+(k+i*dim) )*( *( mat2+(j+k*dim) ));
      }
    }
  }
  

} 

void mat_mult_ser(double *mat1, double *mat2, double *result, int dim, int nthr)
{

 
  
    int i,j,k, count=0; 
    
    for(i=0; i<dim; i++)
    {
      
      for(j=0; j<dim; j++) 
      {
        *( result+(j+i*dim) ) = 0; 
      	
        
        for(k=0; k<dim; k++)
          *( result+(j+i*dim) ) += *( mat1+(k+i*dim) )*( *( mat2+(j+k*dim) ));
      }
    }
    
}


void print_matrix(double *mat, int dim)
{
	int i,j; 
 
  	for(i=0; i<dim; i++) 
    		for(j=0; j<dim; j++)
			printf("%lf ", *(mat+(j+i*dim)));
	printf("\n");
}


void check_results(double *result_ser, double *result_para, int dim)
{
	int flag = 1, i=0, j=0;
	for(i=0; i<dim; ++i)
		for(j=0; j<dim; ++j)
			if(*(result_ser+(j+i*dim)) != *(result_para+(j+i*dim)))
			{
				flag = 0;
				break;
			}
	if(flag)
		printf("\nResults match.\n");
	else
		printf("\nResults do not match.\n");
}

void dump_file(double *mat, int ch, int dim)
{
	int i, j;
	FILE *fp;
	if(ch==0)
		fp = fopen("wserial.txt", "w");
	else
		fp = fopen("wparallel.txt", "w");
	for(i=0; i<dim; ++i)
	{
		for(j=0; j<dim; ++j)
			fprintf(fp, "%lf  ", *(mat+j+i*dim));
		fprintf(fp, "\n");
	}

}

int main(int argc, char *argv[])
{
	int i, dim=1000, nthr=4, choice=0;
	double *mat1 = (double*)calloc(dim*dim, sizeof(double));
	double *mat2 = (double*)calloc(dim*dim, sizeof(double));
	
	double *result_ser = (double*)calloc(dim*dim, sizeof(double));
	double *result_para = (double*)calloc(dim*dim, sizeof(double));
	
	if(argc == 2)
		choice = atoi(argv[1]);
	else if(argc>2)
	{
    		dim=atoi(argv[1]);
		if(argc>2)
			nthr = atoi(argv[2]);
		else
			nthr=4;
		if(argc>3)
			choice = atoi(argv[3]);
		else
			choice = 0;
	
	}
	else
	{
		nthr = 4;
		dim = 1000;
		choice = 0;
	}
	
	printf("Program parameters:\nNumber of threads : %d\nMatrix dimensions : %d\n", nthr, dim);

	rand_mat(mat1, dim);
	rand_mat(mat2, dim);

	time_t start =0, end=0;
	time(&start);
	mat_mult_ser(mat1, mat2, result_ser, dim, nthr);
	time(&end);
	double time_taken_ser = difftime(end, start);
	dump_file(result_ser, 0, dim);
	
	time(&start);
	mat_mult_thr(mat1, mat2, result_para, dim, nthr, choice);
	time(&end);
	double time_taken_para = difftime(end, start);
	dump_file(result_para, 1, dim);

	//print_matrix(result, 10);
	check_results(result_ser, result_para, dim);
	printf("\nTime taken:\nSerial : %lf\nParallel : %lf\n", time_taken_ser, time_taken_para);

	//printf("%f\n", time_taken_para);
		
}
	


