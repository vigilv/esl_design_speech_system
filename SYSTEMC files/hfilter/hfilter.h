
//||||||||||||||||||||| HEADER FILE (hfilter.h) ||||||||||||||||||||||

#include<systemc.h>
#include<math.h>
#include<stdio.h>

SC_MODULE(hfilter)
{
	sc_event bitrevre,sta,butter,revimag,maglog,butteri,bitrevi;
	sc_event div,bitrevfil,butterfil,revimagfil;
	sc_event expon,bitrevfinal,butterfinal;

	SC_CTOR(hfilter);

	void bit_reverse_real(void);
	void num_stage(void);
	void butterfly(void);
	void bit_reverse_img(void);
	void abslog(void);
	void butterflyi(void);
	void bit_reverse_reali(void);
	void divN(void);
	void bit_reverse_real_fil(void);
	void butterfly_fil(void);
	void bit_reverse_img_fil(void);
	void expo(void);
	void bit_reverse_real_final(void);
	void butterflyfinal(void);
	void divNfinal(void);	
};

