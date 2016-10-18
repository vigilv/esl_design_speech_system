
//|||||||||||||||||| IMPLEMENTATION FILE (hfilter.h) |||||||||||||||||

#include "hfilter.h"
#define PI 3.1428571


int N,stage;
float xr[1024],xi[1024],ablog[1024],ablogi[1024],res[256],cep[1024],fin[1024];

hfilter::hfilter(sc_module_name hfil_module):sc_module(hfil_module)
{
	SC_THREAD(bit_reverse_real);
	SC_THREAD(bit_reverse_img);
	SC_THREAD(num_stage);
	SC_THREAD(butterfly);
	SC_THREAD(abslog);
	SC_THREAD(bit_reverse_reali);
	SC_THREAD(butterflyi);
	SC_THREAD(divN);
	SC_THREAD(bit_reverse_real_fil);
	SC_THREAD(butterfly_fil);
	SC_THREAD(bit_reverse_img_fil);
	SC_THREAD(expo);
	SC_THREAD(bit_reverse_real_final);
	SC_THREAD(butterflyfinal);
	SC_THREAD(divNfinal);

}

//****************************** ANALYZER ****************************

//-------- BIT REVERSE REAL COEFFICIENTS FOR FFT OF ANALYZER ---------

void hfilter::bit_reverse_real(void)
{
	wait(SC_ZERO_TIME);
	int i,j,k;
	float tr;
	j=0;
	for(i=0;i<(N-1);i++)
	{
		if(i<j)
		{
			tr=xr[j];			
			xr[j]=xr[i];		
			xr[i]=tr;
		}
		k=N/2;
		while(k<=j)
		{
			j=j-k;
			k=k/2;
		}
		j=j+k;
	}
	bitrevre.notify();
}

//---------- FIND THE NUMBER OF STAGES FOR FFT OF ANALYZER -----------

void hfilter::num_stage(void)
{
	wait(SC_ZERO_TIME);
	wait(bitrevre);
	
	int irem;
	stage=0;
	irem=N;
	while (irem>1)
	{
		irem=irem/2;
		stage=stage+1;
	}
	sta.notify();
}

//------------- BUTTERFLY COMPUTATION FOR FFT OF ANALYZER ------------

void hfilter::butterfly(void)
{
	wait(SC_ZERO_TIME);
	wait(sta);

	int M,M2,l,j,i,ip,sign=1;
	float ur,ui,wr,wi,temp,tr,ti;
	M=1;
	for(l=1;l<=stage;l++)
	{
		M2=M;
		M=M*2;
		ur=1.0;
		ui=0;
		wr=cos(PI/M2);
		wi=sign*sin(PI/M2);
		for(j=0;j<M2;j++)
		{
			i=j;
			while(i<N)
			{
				ip=i+M2;
				tr=xr[ip]*ur-xi[ip]*ui;
				ti=xi[ip]*ur+xr[ip]*ui;
				xr[ip]=xr[i]-tr;
				xi[ip]=xi[i]-ti;
				xr[i]=xr[i]+tr;
				xi[i]=xi[i]+ti;
				i=i+M;
			}
			temp=ur*wr-ui*wi;
			ui=ui*wr+ur*wi;
			ur=temp;
		}
	}
butter.notify();	
}

//------ BIT REVERSE IMAGINARY COEFFICIENTS FOR FFT OF ANALYZER -------

void hfilter::bit_reverse_img(void)
{
	wait(SC_ZERO_TIME);
	wait(butter);
	int i;
	float ti;
	for(i=1;i<(N/2);i++)
	{		
			ti=xi[i];
			xi[i]=xi[N-i];
			xi[N-i]=ti;
	}
	revimag.notify();
}

//-------------- MAGNITUDE AND LOGARITHM IN ANALYZER ----------------

void hfilter::abslog(void)
{
	wait(SC_ZERO_TIME);
	wait(revimag);
	int i;
	float mag,temp;
	for(i=0;i<N;i++)
	{
		mag=sqrt(pow((xr[i]),2)+pow((xi[i]),2));
		if (mag==0)
		{
			mag=1;
		}
		temp=log(mag);
		ablog[i]=temp;
	}
	maglog.notify();
}

//-------- BIT REVERSE REAL COEFFICIENTS FOR IFFT OF ANALYZER --------

void hfilter::bit_reverse_reali(void)
{
	wait(SC_ZERO_TIME);
	wait(maglog);
	int i,j,k;
	float tr;
	j=0;
	for(i=0;i<(N-1);i++)
	{
		if(i<j)
		{
			tr=ablog[j];			
			ablog[j]=ablog[i];		
			ablog[i]=tr;						
		}
		k=N/2;
		while(k<=j)
		{
			j=j-k;
			k=k/2;
		}
		j=j+k;
	}
	for(i=0;i<N;i++)
	{
		ablogi[i]=0.0;
	}
	bitrevi.notify();
}

//----------- BUTTERFLY COMPUTATION FOR IFFT OF ANALYZER ------------

void hfilter::butterflyi(void)
{
	wait(SC_ZERO_TIME);
	wait(bitrevi);
	int M,M2,l,j,i,ip,sign=-1;
	float ur,ui,wr,wi,temp,tr,ti;
	M=1;
	for(l=1;l<=stage;l++)
	{
		M2=M;
		M=M*2;
		ur=1.0;
		ui=0;
		wr=cos(PI/M2);
		wi=sign*sin(PI/M2);
		for(j=0;j<M2;j++)
		{
			i=j;
			while(i<N)
			{
				ip=i+M2;
				tr=ablog[ip]*ur-ablogi[ip]*ui;
				ti=ablogi[ip]*ur+ablog[ip]*ui;
				ablog[ip]=ablog[i]-tr;
				ablogi[ip]=ablogi[i]-ti;
				ablog[i]=ablog[i]+tr;
				ablogi[i]=ablogi[i]+ti;
				i=i+M;
			}
			temp=ur*wr-ui*wi;
			ui=ui*wr+ur*wi;
			ur=temp;
		}
	}
butteri.notify();	
}

//------------ DIVIDE EACH TERM BY N FOR IFFT OF ANALYZER ------------

void hfilter::divN(void)
{
	wait(SC_ZERO_TIME);
	wait(butteri);
	int i;
	float temp;
	for(i=0;i<N;i++)
	{
		ablog[i]=ablog[i]/N;
		temp=ablog[i];
		cep[i]=temp;
	}
	div.notify();
}

//**************************** SYNTHESIZER ***************************

//------- BIT REVERSE REAL COEFFICIENTS FOR FFT OF SYNTHESIZER --------

void hfilter::bit_reverse_real_fil(void)
{
	wait(SC_ZERO_TIME);
	wait(div);
	int i,j,k;
	float tr;
	j=0;
	for(i=0;i<(N-1);i++)
	{
		if(i<j)
		{
			tr=ablog[j];			
			ablog[j]=ablog[i];		
			ablog[i]=tr;
		}
		k=N/2;
		while(k<=j)
		{
			j=j-k;
			k=k/2;
		}
		j=j+k;
	}
	for(i=0;i<N;i++)
	{
		xr[i]=ablog[i];
	}
	for(i=0;i<N;i++)
	{
		xi[i]=0.0;
	}
	bitrevfil.notify();
}

//---------- BUTTERFLY COMPUTATION FOR FFT OF SYNTHESIZER ------------

void hfilter::butterfly_fil(void)
{
	wait(SC_ZERO_TIME);
	wait(bitrevfil);

	int M,M2,l,j,i,ip,sign=1;
	float ur,ui,wr,wi,temp,tr,ti;
	M=1;
	for(l=1;l<=stage;l++)
	{
		M2=M;
		M=M*2;
		ur=1.0;
		ui=0;
		wr=cos(PI/M2);
		wi=sign*sin(PI/M2);
		for(j=0;j<M2;j++)
		{
			i=j;
			while(i<N)
			{
				ip=i+M2;
				tr=xr[ip]*ur-xi[ip]*ui;
				ti=xi[ip]*ur+xr[ip]*ui;
				xr[ip]=xr[i]-tr;
				xi[ip]=xi[i]-ti;
				xr[i]=xr[i]+tr;
				xi[i]=xi[i]+ti;
				i=i+M;
			}
			temp=ur*wr-ui*wi;
			ui=ui*wr+ur*wi;
			ur=temp;
		}
	}
butterfil.notify();	
}

//----- BIT REVERSE IMAGINARY COEFFICIENTS FOR FFT OF SYNTHESIZER ----

void hfilter::bit_reverse_img_fil(void)
{
	wait(SC_ZERO_TIME);
	wait(butterfil);
	int i;
	float ti;
	for(i=1;i<(N/2);i++)
	{		
			ti=xi[i];
			xi[i]=xi[N-i];
			xi[N-i]=ti;
	}
	revimagfil.notify();
}

//------------------- EXPONENTIATION IN SYNTHESIZER ------------------

void hfilter::expo(void)
{
	wait(SC_ZERO_TIME);
	wait(revimagfil);
	int i;
	float temp,mag;
	for(i=0;i<N;i++)
	{
		mag=sqrt(pow((xr[i]),2)+pow((xi[i]),2));
		temp=exp(mag);
		xr[i]=temp;
	}
	expon.notify();
}

//------ BIT REVERSE REAL COEFFICIENTS FOR IFFT OF SYNTHESIZER -------

void hfilter::bit_reverse_real_final(void)
{
	wait(SC_ZERO_TIME);
	wait(expon);
	int i,j,k;
	float tr;
	j=0;
	for(i=0;i<(N-1);i++)
	{
		if(i<j)
		{
			tr=xr[j];			
			xr[j]=xr[i];		
			xr[i]=tr;
		}
		k=N/2;
		while(k<=j)
		{
			j=j-k;
			k=k/2;
		}
		j=j+k;
	}

	for(i=0;i<N;i++)
	{
		xi[i]=0.0;
	}
	bitrevfinal.notify();
}

//---------- BUTTERFLY COMPUTATION FOR IFFT OF SYNTHESIZER -----------

void hfilter::butterflyfinal(void)
{
	wait(SC_ZERO_TIME);
	wait(bitrevfinal);
	int M,M2,l,j,i,ip,sign=-1;
	float ur,ui,wr,wi,temp,tr,ti;
	M=1;
	for(l=1;l<=stage;l++)
	{
		M2=M;
		M=M*2;
		ur=1.0;
		ui=0;
		wr=cos(PI/M2);
		wi=sign*sin(PI/M2);
		for(j=0;j<M2;j++)
		{
			i=j;
			while(i<N)
			{
				ip=i+M2;
				tr=xr[ip]*ur-xi[ip]*ui;
				ti=xi[ip]*ur+xr[ip]*ui;
				xr[ip]=xr[i]-tr;
				xi[ip]=xi[i]-ti;
				xr[i]=xr[i]+tr;
				xi[i]=xi[i]+ti;
				i=i+M;
			}
			temp=ur*wr-ui*wi;
			ui=ui*wr+ur*wi;
			ur=temp;
		}
	}
butterfinal.notify();	
}

//---------- DIVIDE EACH TERM BY N FOR IFFT OF SYNTHESIZER -----------

void hfilter::divNfinal(void)
{
	wait(SC_ZERO_TIME);
	wait(butterfinal);
	int i;
	for(i=0;i<N;i++)
	{
		fin[i]=xr[i]/N;
	
	}

}

//****************************** MAIN PROGRAM ************************ 

int sc_main(int argc, char* argv[])
{

	int i,Fs;
	float temp;
	FILE *fp;
	FILE *fg;

	struct in
	{
		float real;
		int real1;
	};
	struct in invalues;

	struct ceps
	{
		float cep;
	};
	struct ceps cepstrum;

	struct res
	{
		float re;
		int   re1;
	};
	struct res final;

// Get values from the speech file

	fp= fopen("speech_values.txt","r");	
	fscanf(fp," %d",&invalues.real1);
	N=invalues.real1;

	for(i=1;i<N;i++)
	{		
		if(i==1)
		{
			fscanf(fp," %d",&invalues.real1);
			Fs=invalues.real1;
		}
		else
		{
		fscanf(fp," %f",&invalues.real);
		xr[i]=invalues.real;
		}
	}
	fclose(fp);

	for(i=0;i<N;i++)
	{
		temp=xr[i+2];
		xr[i]=temp;
	}

//Start the SYSTEMC simulation engine after instantiating the module

	hfilter hfilter_inst("hfilter_inst");
	sc_start();

// Write cepstrum values to file

	fp= fopen("cepstrum.txt","w");

	for(i=0;i<N;i++)
	{
		cepstrum.cep=cep[i];
		fprintf(fp,"  %f",cepstrum.cep);
	}
	fclose(fp);

// Write final values to file
	
	fg= fopen("finvalues.txt","w");

	final.re1=N;
	fprintf(fg,"  %d",final.re1);
	final.re1=Fs;
	fprintf(fg,"  %d",final.re1);

	for(i=0;i<N;i++)
	{
			
		final.re=fin[i];
		fprintf(fg,"  %f",final.re);
	}
	fclose(fg);
	
	return 0;
}

