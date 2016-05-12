
clear
close
clc

%Run for all mic and speaker positions
for ii=1:6 
    for jj=1:3 


fs=2000;%sample rate


%mic and speaker positions
M1=[1.475,2.030,1.13];
M2=[0.725,2.03,1.13];
M3=[2.225,2.03,1.65];
M4=[0.725,3.03,1.65];
M5=[1.475,3.03,1.65];
M6=[2.225,3.03,1.13];


S1=[0.7,0.45,1.33];
S2=[1.730,0.45,1.33];
S3=[1.2,1.3,1.03];

%mic and speaker positions in vectors
MALL=zeros(6,3); 
MALL(1,:)=M1;
MALL(2,:)=M2;
MALL(3,:)=M3;
MALL(4,:)=M4;
MALL(5,:)=M5;
MALL(6,:)=M6;

SALL=zeros(6,3); 
SALL(1,:)=S1;
SALL(2,:)=S2;
SALL(3,:)=S3;

%room dim
Lx=2.950;
Ly=6.105;
Lz=2.260;


Dims = [Lx,Ly,Lz];        % Room dimensions [lx, ly, lz] in meters
speaker = SALL(jj,:);        % Source position1 [x,y,z] in meters 
mic = MALL(ii,:);

%absorption coefficients

Alphaav31=(Lx*Ly*Lz*24*log10(10))/(343*RT60AT31*2*((Lx*Ly)+(Lx*Lz)+(Ly*Lz)));
Alphaav63=(Lx*Ly*Lz*24*log10(10))/(343*RT60AT63*2*((Lx*Ly)+(Lx*Lz)+(Ly*Lz)));
Alphaav125=(Lx*Ly*Lz*24*log10(10))/(343*RT60AT125*2*((Lx*Ly)+(Lx*Lz)+(Ly*Lz)));
Alphaav250=(Lx*Ly*Lz*24*log10(10))/(343*RT60AT250*2*((Lx*Ly)+(Lx*Lz)+(Ly*Lz)));

A31 = Alphaav31; 
A63 = Alphaav63;
A125 = Alphaav125; 
A250 = Alphaav250;
B31 = Alphaav31; 
B63 = Alphaav63;
B125 = Alphaav125; 
B250 = Alphaav250;
C31 = Alphaav31; 
C63 = Alphaav63;
C125 = Alphaav125; 
C250 = Alphaav250;
D31 = Alphaav31; 
D63 = Alphaav63;
D125 = Alphaav125; 
D250 = Alphaav250;
E31 = Alphaav31; 
E63 = Alphaav63;
E125 = Alphaav125; 
E250 = Alphaav250;
F31 = Alphaav31; 
F63 = Alphaav63;
F125 = Alphaav125; 
F250 = Alphaav250;

%linear interpolation of values between these points

A1=zeros(1,250);
B1=zeros(1,250);
C1=zeros(1,250);
D1=zeros(1,250);
E1=zeros(1,250);
F1=zeros(1,250);

gradA1=(A63-A31)/32;
gradA2=(A125-A63)/62;
gradA3=(A250-A125)/125;

gradB1=(A63-A31)/32;
gradB2=(A125-A63)/62;
gradB3=(A250-A125)/125;

gradC1=(A63-A31)/32;
gradC2=(A125-A63)/62;
gradC3=(A250-A125)/125;

gradD1=(A63-A31)/32;
gradD2=(A125-A63)/62;
gradD3=(A250-A125)/125;

gradE1=(A63-A31)/32;
gradE2=(A125-A63)/62;
gradE3=(A250-A125)/125;

gradF1=(A63-A31)/32;
gradF2=(A125-A63)/62;
gradF3=(A250-A125)/125;


for i=1:63
A1(i)=i*gradA1;
end
for i=64:125
A1(i)=i*gradA2;
end
for i=126:250
A1(i)=i*gradA3;
end

for i=1:63
B1(i)=i*gradB1;
end
for i=64:125
B1(i)=i*gradB2;
end
for i=126:250
B1(i)=i*gradB3;
end

for i=1:63
C1(i)=i*gradC1;
end
for i=64:125
C1(i)=i*gradC2;
end
for i=126:250
C1(i)=i*gradC3;
end

for i=1:63
D1(i)=i*gradD1;
end
for i=64:125
D1(i)=i*gradD2;
end
for i=126:250
D1(i)=i*gradD3;
end

for i=1:63
E1(i)=i*gradE1;
end
for i=64:125
E1(i)=i*gradE2;
end
for i=126:250
E1(i)=i*gradE3;
end

for i=1:63
E1(i)=i*gradF1;
end
for i=64:125
E1(i)=i*gradF2;
end
for i=126:250
E1(i)=i*gradF3;
end


%run for freq of 1 to 250
Ztot=zeros(2,10000);

for i=1:250
w=i;
Z=FDTDmethod(w, Dims, speaker, mic, A1(w), B1(w), C1(w), D1(w), E1(w), F1(w)); 
Ztot=Ztot+Z;
end



% take fft
Zfft = fft(Ztot(2,:));
      
% keep only meaningful frequencies
NFFT = length(Z(2,:));
if mod(NFFT,2)==0
Nout = (NFFT/2)+1;
else
Nout = (NFFT+1)/2;
end
Zfft = Zfft(1:Nout);
freq1 = ((0:Nout-1)./NFFT).*fs;
% put into dB
Zfft = 20*log10(abs(Zfft)./NFFT);
%       
% smooth
Noct = 12;
R = smooth_spectrum(Zfft,freq1,Noct);
       
% plot
figure(1)
axis([1 250 -40 20])
plot(freq1,Zfft,freq1,R)
grid on

%save data 
%filenamestring=sprintf('FDTDMX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);
filenamestring=sprintf('FDTDVARCORNERMX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);
save(filenamestring,'freq1','Zfft','R');

    end
end
