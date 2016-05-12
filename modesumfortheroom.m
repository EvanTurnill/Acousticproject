%This is the model of the room measured in the report

%Run for all 18 combinations of positions
for ii=1:6
    for jj=1:3


Lx=2.950;
Ly=6.105;
Lz=2.260;

M1=[1.475,2.030,1.13];
M2=[0.725,2.03,1.13];
M3=[2.225,2.03,1.65];
M4=[0.725,3.03,1.65];
M5=[1.475,3.03,1.65];
M6=[2.225,3.03,1.13];


S1=[0.7,0.45,1.33];
S2=[1.730,0.45,1.33];
S3=[1.2,1.3,1.03];

%Absoprtion coefficients calculated from the RT60 measured in the room
A31 = 0.04; 
A63 = 0.06;
A125 = 0.08; 
A250 = 0.1;
B31 = 0.04; 
B63 = 0.06;
B125 = 0.08; 
B250 = 0.1;
C31 = 0.04; 
C63 = 0.06;
C125 = 0.08; 
C250 = 0.1;
D31 = 0.04; 
D63 = 0.06;
D125 = 0.08; 
D250 = 0.1;
E31 = 0.1; 
E63 = 0.16;
E125 = 0.2; 
E250 = 0.25;
F31 = 0.1; 
F63 = 0.16;
F125 = 0.2; 
F250 = 0.25;

A1=zeros(1,250);
B1=zeros(1,250);
C1=zeros(1,250);
D1=zeros(1,250);
E1=zeros(1,250);
F1=zeros(1,250);


%Linear interpolation of the absorption coefficients 
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


%Set up vectors for all mic and speaker positions 
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

Rx=MALL(ii,1);
Ry=MALL(ii,2);
Rz=MALL(ii,3);
Sx=SALL(jj,1);
Sy=SALL(jj,2);
Sz=SALL(jj,3);



%calculate for frequencies between 1 and 250 Hz
X=1:250;

for F=1:length(X)
X(2,F)= 70+10*abs(modesummethod(Lx,Ly,Lz,X(1,F), A1(w), B1(w), C1(w), D1(w), E1(w), F1(w), Rx,Ry,Rz,Sx,Sy,Sz));
end

figure(1)
axis([1 250 0 100])
plot(X(1,:),X(2,:))
grid on

%save data created
filenamestring=sprintf('MODALVARWALLS_MX%dSX%d_1to250_X_row1freq_row2amp.mat',ii,jj);
save(filenamestring,'X');

    end
end