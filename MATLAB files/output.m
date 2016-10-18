
%Get cepstrum values from file

fid=fopen('cepstrum.txt','r');
cepo=fscanf(fid,'%f');
fclose(fid);

%Get final output values from file

fid=fopen('finvalues.txt','r');
No=fscanf(fid,'%d',1);
Fso=fscanf(fid,'%d',1);
fino=fscanf(fid,'%f');
fclose(fid);

Tso=(1/Fso);
no=[(1*Tso):Tso:(No*Tso)];

% Plot cepstrum

figure,plot(no,cepo);
title('Cepstrum');
xlabel('Quefrency in seconds');
ylabel('Amplitude');

%Plot final output waveform

figure,plot(no,fino);
title('Output speech waveform');
xlabel('Time in seconds');
ylabel('Amplitude');

wavplay(fino);
