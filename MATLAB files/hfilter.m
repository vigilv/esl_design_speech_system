
% Get the speech values stored in the file

fid=fopen('speech_values.txt','r');
N1=fscanf(fid,'%d',1);
Fs1=fscanf(fid,'%d',1);
val=fscanf(fid,'%f');

%Speech analysis

trans=fft(val,N1);
mag=abs(trans);
sep=log(mag);
cepstrum=ifft(sep,N1);

%Store cepstrum values in file

fid=fopen('cepstrum.txt','w');
fprintf(fid,' %f ',cepstrum);
fclose(fid);

%Speech synthesis

inver=fft(cepstrum);
ab=abs(inver);
com=exp(ab);
fin=ifft(com);

%Store final values in file
fid=fopen('finvalues.txt','w');
fprintf(fid,' %d',N1);
fprintf(fid,' %d',Fs1);
fprintf(fid,' %f ',fin);
fclose(fid);
