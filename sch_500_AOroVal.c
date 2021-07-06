#include "complex.h"
#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"

#define N 500 //parametro que discretiza la posicion
#define n 90 //parametro que discretiza el tiempo
#define nc 30 //parametro que define el numero de ciclos
#define PI 3.141592
#define NUM_SIM 1000
#define CF 1000 //valor para que no se quede pillado el programa

gsl_rng*tau;


int main()
{
  FILE *f1,*f2;
  int i,j,k;
  fcomplex ao,c,Amenos,Amas,scomplejo,kocomplejo,auxgamma,auxgamma2,auxbeta1,auxbeta2;
  fcomplex b[N+1],phi[N+1],Vcomplejo[N+1],Ao[N+1],beta[N+1],chi[N+1],gamma[N+1],alfa[N+1];
  fcomplex auxnc,auxko,auxs1,auxs2,auxb,auxb2,auxchi,auxalfa,n1,n2;
  double mod,fase,suma;
  double ko,s,modphi[N+1],modphi2[N+1];
  double V[N+1],lambda,modbeta;
  extern gsl_rng*tau;

  int mt,num_sim;
  double x,y,z,transmision,norma,ND,n_d;

  int contadormt,cortafuegos;

  int semilla=141592;
  tau=gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(tau,semilla);

  f1=fopen("transmision_500.txt","w");
 
  for(z=0;z<5;z++)
    {
      ND=10.0;
      for(y=0;y<6;y++)
	{
	  if(y==0)
	    lambda=0.1;
	  if(y==1)
	    lambda=0.3;
	  if(y==2)
	    lambda=0.5;
	  if(y==3)
	    lambda=1.0;
	  if(y==4)
	    lambda=5.0;
	  if(y==5)
	    lambda=10.0;
	  mt=0.0;
	  for(num_sim=0;num_sim<NUM_SIM;num_sim++)
	    {
	      n1=Complex(1.0,0.0);
	      n2=Complex(-1.0,0.0);
	      ko=2*PI*nc/N;
	      kocomplejo=Complex(ko,0.0);
	      
	      s=1/(4*pow(ko,2));
	      scomplejo=Complex(s,0.0);
	      
	      ao=Complex(-2.0,2.0/s);
	      
	      Amas=Complex(1.0,0.0);
	      Amenos=Complex(1.0,0.0);
	      //CONDICIONES DE CONTORNO
	      phi[0]=Complex(0.0,0.0);
	      phi[N]=Complex(0.0,0.0);
	      modphi[0]=Cabs(phi[0]);
	      modphi2[0]=pow(modphi[0],2);
	      modphi[N]=Cabs(phi[N-1]);
	      modphi2[N]=pow(modphi[N-1],2);
	      
	      alfa[N-1]=Complex(0.0,0.0);
	      
	      auxb=Complex(0.0,4.0);
	      
	      beta[N-1]=Complex(0.0,0.0);
	      
	      
	      chi[0]=Complex(0.0,0.0);
	      chi[N]=Complex(0.0,0.0);
	      for(j=N-1;j>0;j--)
		{
		  if(j<=3*N/5 & j>=2*N/5)
		    {
		      V[j]=lambda*pow(ko,2);
		      Vcomplejo[j]=Complex(V[j],0.0);
		    }
		  else
		    {
		      V[j+1]=0.0;
		      Vcomplejo[j]=Complex(0.0,0.0);
		    }
		  
		  Ao[j]=Csub(ao,Vcomplejo[j]);
		  
		  auxgamma=Cmul(Amas,alfa[j]);
		  gamma[j]=Cadd(Ao[j],auxgamma);
		  
		  auxalfa=Cdiv(Amenos,gamma[j]);
		  alfa[j-1]=RCmul(-1,auxalfa);
		  
		  mod=exp(-8*pow(4*j-N,2)/pow(N,2));
		  fase=ko*j;		
		  phi[j].r=mod*cos(fase);
		  phi[j].i=mod*sin(fase);
		  
		}
	      norma=0.0;
	      for(j=0;j<=N;j++)//NORMLIZAMOS LA FUNCION INICIAL
		{
		  modphi[j]=Cabs(phi[j]);
		  modphi2[j]=modphi[j]*modphi[j];
		  norma=norma+modphi2[j];
		}
	      norma=sqrt(norma);
	      for(j=0;j<=N;j++)
		{
		  phi[j].r=phi[j].r/norma;
		  phi[j].i=phi[j].i/norma;
		  modphi[j]=Cabs(phi[j]);
		  modphi2[j]=pow(modphi[j],2);
		}
	      contadormt=0;
	      cortafuegos=0;
	      n_d=0.0;
	      while(contadormt==0)//ESCIRBIR LA CONDICION DEL WHILE
		{
		  for(j=N-1;j>0;j--)//CALCULO DE LAS BETAS CON RECURRENCIA DECRECIENTE
		    {
		      auxb2=Cdiv(phi[j],scomplejo);
		      b[j]=Cmul(auxb,auxb2);
		      
		      auxbeta1=Cmul(Amas,beta[j]);
		      auxbeta2=Csub(b[j],auxbeta1);
		      beta[j-1]=Cdiv(auxbeta2,gamma[j]);	  
		    }
		  for(j=0;j<=N;j++)//CALCULO DE LAS CHI Y DE LAS PHI CON RECURRENCIA CRECIENTE
		    {
		      auxchi=Cmul(alfa[j],chi[j]);
		      chi[j+1]=Cadd(auxchi,beta[j]);	  
		    }     
		  for(j=1;j<N;j++)
		    {
		      phi[j]=Csub(chi[j],phi[j]);//resultado de phi
		      modphi[j]=Cabs(phi[j]);
		      modphi2[j]=pow(modphi[j],2);//modulo de phi al cuadrado
		    }
		  suma=0.0;
		  for(j=0;j<=N;j++)
		    {
		      suma=modphi2[j]+suma;	  
		    }      
		  n_d=n_d+1.0;//EVOLCIONAMOS EL SISTEMA HASTA QUE HAYAMOS DADO n_D PASOS
		  if(n_d==ND)
		    {
		      //calculamos P_D	  
		      suma=0.0;
		      for(j=4*N/5;j<N+1;j++)
			{
			  suma=suma+modphi2[j];
			}
		      x=gsl_rng_uniform(tau);//generamos un numero aleatorio y comparamos con P_D
		      
		      if(x<suma)
			{
			  mt=mt+1;//aumentamos mt
			  contadormt=contadormt+1;
			  //Generamos la funcion de onda inicial
			}
		      else//si no recalculamos las funciones de onda
			{
			  for(j=4*N/5;j<N+1;j++)
			    {
			      phi[j].r=0;
			      phi[j].i=0;
			    }
			  k=0.0;
			  for(j=0;j<N;j++)
			    k=k+modphi2[i];
			  for(j=0;j<N;j++)
			    phi[j]=RCmul(1/sqrt(k),phi[j]);
			  //CALCULAMOS LA P_I
			  suma=0.0;
			  for(j=0;j<N/5+1;j++)
			    {
			      suma=suma+modphi2[j];
			    }
			  x=gsl_rng_uniform(tau);//generamos un numero aleatorio y comparamos con P_I
			  
			  if(x<suma)
			    {
			      contadormt=contadormt+1;//aumentamos mt
			      //RECALCULAMOS LA FUNCION DE ONDA INICIAL 
			    }
			  else
			    {
			      for(j=0;j<N/5+1;j++)
				{
				  phi[j].r=0;
				  phi[j].i=0;
				}
			      k=0.0;
			      for(j=0;j<N;j++)
				k=k+modphi2[i];
			      for(j=0;j<N;j++)
				phi[j]=RCmul(1/sqrt(k),phi[j]);
			    }
			}
		      n_d=0.0;    
		    }
		  cortafuegos=cortafuegos+1;
		  if(cortafuegos==CF)//SI HAY UN NUMERO CF DE VECES QUE NO DETECTA NADA RECALCULAMOS LA FUNCION DE ONDA INICIAL Y EMPEZAMOS DE NUEVO
		    {
		      for(j=N-1;j>0;j--)
			{ 
			  mod=exp(-8*pow(4*j-N,2)/pow(N,2));
			  fase=ko*j;		
			  phi[j].r=mod*cos(fase);
			  phi[j].i=mod*sin(fase);
			}
		       norma=0.0;
		       for(j=0;j<=N;j++)//NORMLIZAMOS LA FUNCION INICIAL
			 {
			   modphi[j]=Cabs(phi[j]);
			   modphi2[j]=modphi[j]*modphi[j];
			   norma=norma+modphi2[j];
			 }
		       norma=sqrt(norma);
		       for(j=0;j<=N;j++)
			 {
			   phi[j].r=phi[j].r/norma;
			   phi[j].i=phi[j].i/norma;
			   modphi[j]=Cabs(phi[j]);
			   modphi2[j]=pow(modphi[j],2);
			 }
		    }
		  
		}//FIN DEL WHILE
	    
	    }//FINAL DE LAS SIMULACIONES
	  //CALCULAMOS EL COEFICIENTE DE TRANSMISION
	  transmision=mt/NUM_SIM;
	  fprintf(f1,"%lf\t%lf\n",lambda,transmision);
	}
      fprintf(f1,"\n");//SEPARAMOS LOS RESULTADOS EN FUNCION DE n_D
      ND=ND+10;
    }
      fclose(f1);
      
      return 0;
}
