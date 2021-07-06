#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"

#define N 16 //tamaño de la red

#define e 2.7132
#define IT 1e6 //numero de pasos monte carlo que va a dar
#define M 1 //valor que define el modo de inicio. 1 es ordenado y 2 aleatorio
#define C 3 //valor que le damos al spin -1 para que gnuplot plotee bien

gsl_rng*tau;


//FUNCION PRINCIPAL
int main()
{
  FILE *f1,*f2,*f3;
  FILE *fm,*fe,*fc,*ff;
  int i,j,k,l,a,s[N+1][N+1],auxi,auxj,contador;
  double x,y,z,p[N+1][N+1],min[N+1][N+1],E[N+1][N+1],g;

  //declaramos las variables de los parametros que hay que calcular
  double m_N,e_N,c_N,f[N+1],T,Ep,Ep2;
  double sumam,sumae,sumac,sumaf[N+1],sumafp[N+1],fcorrelacion,contador2,contadorprom;
  double mediae,mediae2,mediam,mediaf[N+1];//varibles donde se guardaran los valores medios necesarios para el calculo de errores
  int num,resto,auxf; //num=valor de i de la funcion de correlacion
  extern gsl_rng*tau;


    int semilla=141592;//cambiamos la semilla para cada valor de N para asegurarnos de que son configuraciones totalmente distintas
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);
    //ficheros donde se recogen los datos necesarios para hacer el calculo de errores
    f1=fopen("es_16.txt","w");
    f2=fopen("escuadrado_16.txt","w");
    f3=fopen("sumamagne_16.txt","w");
   
    //ficheros donde se recogen los resultados que usaremos en el tratamiento de datos
    fm=fopen("magne_16.txt","w");
    fe=fopen("energia_16.txt","w");
    fc=fopen("calor_16.txt","w");
    ff=fopen("correlacion_16.txt","w");
    
   
    T=1.5; //definimos la temperatura inicial
    for(x=0;x<10;x++)
      {
	//ELEGIMOS SI QUEREMOS QUE EMPICE ORDENADO O NO
	if(M==1)//ordenado
	  {
	    for(i=0;i<N+1;i++)
	      {
		for(j=0;j<N+1;j++)
		  {
		    s[i][j]=1;
        
		  }
	      }
	  }
	if(M==2)//aleatorio
	  {
	    for(i=0;i<N+1;i++)
	      {
		for(j=0;j<N+1;j++)
		  {
		    a=gsl_rng_uniform_int(tau,2);
		    s[i][j]=a;
		    if(s[i][j]==0)
		      s[i][j]=-1;
		  }
	      }
	  }
 
	contador=0;
	contador2=1;//contador que usaremos para ver cuándo se han hecho los 100PMC
	contadorprom=0.0;//contador que usaremos para calcular los promedios
	sumam=0.0;
	sumae=0.0;
	sumac=0.0;
	Ep=0.0;//igualamos a cero las variables donde acumulamos los valores necesarios de energia,magnetizacion y calor especifico cada vez que cambiamos de T
	m_N=0.0;
	Ep2=0.0;
	for(num=0;num<N+1;num++)
	  {
	    sumafp[num]=0.0;
	  }
   
	for(l=0;l<IT;l++)//BUCLE QUE NOS REPITE LOS PASOS DE MONTE CARLO
	  {
	    contador=contador+1;
	    num=gsl_rng_uniform_int(tau,N+1); //generamos el i de la funcion de correlacion de manera aleatorio desde 0 hasta N de forma que recorra toda la red
	    
	    for(k=0;k<=N*N;k++)//BUCLE QUE NOS HACE EL PASO MONTECARLO 
	      {
		i=gsl_rng_uniform_int(tau,N+1);//aqui seleccionamos un punto aleatorio de la red
		j=gsl_rng_uniform_int(tau,N+1);
		g=gsl_rng_uniform(tau);
		      		     
		if(i==0)//APLICAMOS LAS CONDICIONES DE CONTORNO PERIODICAS
		  {
		    s[i][j]=s[N][j];
		    auxi=s[N-1][j];
		  }
		if(i==N)
		  auxi=s[1][j];

		if(j==0)
		  {
		    s[i][j]=s[i][N];
		    auxj=s[i][N-1];
		  }
		if(j==N)
		  auxj=s[i][1];
	    
		//ENERGIAS PARA LAS ESQUINAS DE LA RED
		if(i==0 & j==0)
		  E[i][j]=2*s[i][j]*(s[i+1][j]+auxi+s[i][j+1]+auxj);
		if(i==N & j==N)
		  E[i][j]=2*s[i][j]*(s[i-1][j]+auxi+s[i][j-1]+auxj);
		if(i==0 & j==N)
		  {
		    E[i][j]=2*s[i][j]*(s[i+1][j]+auxi+s[i][j-1]+auxj);

		    //printf("%i\t%i\t%lf\t%i\t%i\t%i\t%i\t%i\n",i,j,E[i][j],s[i][j],s[i+1][j],auxi,s[i][j-1],auxj);
        
		  }
	     
		if(i==N & j==0)
		  E[i][j]=2*s[i][j]*(s[i-1][j]+auxi+s[i][j+1]+auxj);
		//ENERGIAS PARA LOS LADOS DE LA RED SIN LAS ESQUINAS
		if(i==0 & j!=0 & j!=N)
		  E[i][j]=2*s[i][j]*(s[i+1][j]+auxi+s[i][j+1]+s[i][j-1]);
		if(j==0 & i!=0 & i!=N)
		  E[i][j]=2*s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j+1]+auxj);
		if(j==N & i!=0 & i!=N)
		  E[i][j]=2*s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j-1]+auxj);
		if(i==N & j!=0 & j!=N)
		  E[i][j]=2*s[i][j]*(s[i-1][j]+auxi+s[i][j+1]+s[i][j-1]);
		//ENERGIAS PARA EL RESTO DE LA RED
		if(i!=0 & i!=N & j!=0 & j!=N)
		  E[i][j]=2*s[i][j]*(s[i-1][j]+s[i+1][j]+s[i][j+1]+s[i][j-1]);
	    	        
		p[i][j]=exp(-E[i][j]/T);
		if(p[i][j]<1)//con este if-else miramos cual es el minimo entre 1 y la exponencial
		  min[i][j]=p[i][j];
		else
		  min[i][j]=1.0;
		// printf("%lf\n",E[i][j]);
		if(g<min[i][j])//aqui miramos si hay que cambiar el sentido del spin
		  {
		    if(s[i][j]==1)
		      s[i][j]=-1;
		    else
		      if(s[i][j]==-1)
			s[i][j]=1;		
		  }	   
		
	      }//FINAL BUCLE PASOS MONTECARLO
	    
	    if(contador2==100) //hacemos los calculos despues de 100 PMC
	      {
		for(i=0;i<N+1;i++)
		  {
		     for(j=0;j<N+1;j++)
		      {
				//ENERGIAS PARA LAS ESQUINAS DE LA RED
			if(i==0 & j==0)
			  E[i][j]=s[i][j]*(s[i+1][j]+auxi+s[i][j+1]+auxj);
			if(i==N & j==N)
			  E[i][j]=s[i][j]*(s[i-1][j]+auxi+s[i][j-1]+auxj);
			if(i==0 & j==N)
			  {
			    E[i][j]=s[i][j]*(s[i+1][j]+auxi+s[i][j-1]+auxj);       
			  }	     
			if(i==N & j==0)
			  E[i][j]=2*s[i][j]*(s[i-1][j]+auxi+s[i][j+1]+auxj);
			//ENERGIAS PARA LOS LADOS DE LA RED SIN LAS ESQUINAS
			if(i==0 & j!=0 & j!=N)
			  E[i][j]=s[i][j]*(s[i+1][j]+auxi+s[i][j+1]+s[i][j-1]);
			if(j==0 & i!=0 & i!=N)
			  E[i][j]=s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j+1]+auxj);
			if(j==N & i!=0 & i!=N)
			  E[i][j]=s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j-1]+auxj);
			if(i==N & j!=0 & j!=N)
			  E[i][j]=s[i][j]*(s[i-1][j]+auxi+s[i][j+1]+s[i][j-1]);
			//ENERGIAS PARA EL RESTO DE LA RED
			if(i!=0 & i!=N & j!=0 & j!=N)
			  E[i][j]=s[i][j]*(s[i-1][j]+s[i+1][j]+s[i][j+1]+s[i][j-1]);

			sumae=sumae+E[i][j];//sumatoria de la ecuacion 12
		      }
		    
		  }
		sumae=-(1.0/2.0)*sumae;//valor de la energia E(S) ecuacion 12
		Ep=Ep+sumae;//valor de la energia la vamos acumulando
		fprintf(f1,"%lf\n",sumae);//guardamos el valor en un fichero para el calculo de errores
		for(i=0;i<N+1;i++)
		  {
		    for(j=0;j<N+1;j++)
		      {
			sumam=sumam+s[i][j];
		      }
		  }
		sumam=sumam/pow(N,2);//magnetizacion
		m_N=m_N+sumam;//valor de la magnetizacion la vamos acumulando
		fprintf(f3,"%lf\n",sumam);//guardamos el valor en un fichero para el calculo de errores
		sumac=pow(sumae,2);//valor de la energia al cuadrado
		Ep2=Ep2+sumac;//valor de la energia al cuadrado la vamos acumulando
		fprintf(f2,"%lf\n",sumac);//guardamos el valor en un fichero para el calculo de errores
		
		for(num=0;num<N+1;num++) //COMIENZO CALCULO DE SUMA FUNCION DE CORRELACION
		  {
		    sumaf[num]=0.0;
		    for(i=0;i<N+1;i++)
		      {
			for(j=0;j<N+1;j++)
			  {
			    if(i+num>N)
			      {
				resto=abs(i+num-N);
				auxf=s[resto][j];
				sumaf[num]=sumaf[num]+s[i][j]*auxf;
			      }else
			      {
				sumaf[num]=sumaf[num]+s[i][j]*s[i+num][j];
			       //acumulamos los valores de la sumatoria para cada i de f(i)
				
			      }
			  }		
		      } //FINAL CALCULO SUMA FUNCION DE CORRELACION
		  }
		for(num=0;num<N+1;num++)
		  {
		    sumafp[num]=sumafp[num]+sumaf[num];//acumulamos la suma del producto de spines de la funcion de correlacion
		  }
		contador2=0;
		contadorprom=contadorprom+1;
	      }
	    contador2=contador2+1;	     
	  }//FIN DEL BUCLE QUE ITERA LOS PMC
	
      
	//despues de los 10⁴ PMC calculamos los promedios
	Ep=Ep/contadorprom; //promedio del valor de la energia <E(S)>
	e_N=Ep/(2.0*N);//ENERGIA MEDIA
	m_N=m_N/contadorprom;//MAGNETIZACION PROMEDIO
	Ep2=Ep2/contadorprom;//promedio del valor de la energia al cuadrado <E(S)>²
	c_N=(Ep2-pow(Ep,2))/(pow(N,2)*T);//CALOR ESPECIFICO
	for(num=0;num<N+1;num++)
	  {
	    f[num]=sumafp[num]/(contadorprom*pow(N,2));//FUNCION DE CORRELACION	   
	  }
	
	fprintf(fm,"%lf\t%lf\n",m_N,T);
	fprintf(fe,"%lf\t%lf\n",e_N,T);
	fprintf(fc,"%lf\t%lf\n",c_N,T);
	for(num=0;num<N+1;num++)
	  {
	    fprintf(ff,"%lf\t%i\n",f[num],num);
	   
	  }
	fprintf(ff,"\n");

	fprintf(f1,"\n");
	fprintf(f2,"\n");
	fprintf(f3,"\n");
	
    
	
	T=T+0.2;//SALTAMOS DE TEMPERATURA
      }
    
    fclose(f1);
    fclose(f2);
    fclose(f3);
  
    fclose(fm);
    fclose(fe);
    fclose(fc);
    fclose(ff);
    
 
    return 0;
}
    













    
    


