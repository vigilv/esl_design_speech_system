
%  synthetic speech waveform 

N=1024; % N for FFT and IFFT
Fs=8000; % Sampling rate = 8000 Hz

% vocal tract response: v[n]                
				% poles
p1 = .1*exp(j*.3*pi);
p2 = -.8*exp(-j*.8*pi);
p3 = -.2*exp(j*.2*pi);
p4 = .32*exp(-j*.6*pi);

                % zeros
z1 = .83*exp(j*.5*pi);
z2 = -.73*exp(-j*.1*pi);

K = 1;

P = [p1, p2, p3, p4]';
Z = [z1, z2]';

[B,A] = zp2tf(Z,P,K);

impulse = zeros(1, 100);
impulse(1) = 1.0;
v = filter(B,A,impulse);
t = [v(3:length(v)), 0, 0];
v = t;

% glottal pulse: g[n] 
% two poles outside unit circle

eq = (.9).^(0:99);
g = conv(eq, eq);
g = g(length(g):-1:1);
g = g(180:199);

% impulse train: p[n]

p = zeros(1, 1000);
beta = 0.70;
p(1) = 1.0;
p(201) = beta^1;
p(401) = beta^2;
p(601) = beta^3;
p(801) = beta^4;

% glottal pulse train: x[n]

x = conv(p, g);

% speech: y[n] = x[n]*v[n]

y = filter(B,A,x);

wavplay(y);

% Quantize the speech values 

smax=max(y);
slen=length(y);

for i=1:slen
    y(i)=y(i)/(smax+1);
end

%Sample the speech with a sampling rate of 8000 hz

wavwrite(y,Fs,16,'synthetic');
speech=wavread('synthetic');

pause(.75);

original=(smax+1)*speech;
wavplay(original);

%set N for FFT and IFFT

if slen<N
    original(slen:N)=0; 
end
if slen>N
     original=original(1:N);
end

% Display speech waveform

Ts=(1/Fs);
n=[(1*Ts):Ts:(N*Ts)];

plot(n,original);
title('Synthetic speech waveform');
xlabel('Time in seconds');
ylabel('Amplitude');

%Save the values in a file

fid=fopen('speech_values.txt','w');
fprintf(fid,' %d ',N);
fprintf(fid,' %d ',Fs);
fprintf(fid,' %f ',original);
fclose(fid);

