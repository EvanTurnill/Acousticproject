function [u] = modesummethod(Lx,Ly,Lz,F, A1, B1, C1, D1, E1, F1, Rx,Ry,Rz,Sx,Sy,Sz)

% Lx=length of room in x direction
% Ly=length of room in y direction
% Lz=length of room in z direction
% F=frequency
% A1=absorption coefficient on wall at x=0 perp to x-axis
% B1=absorption coefficient on wall at x=Lx perp to x-axis
% C1=absorption coefficient on wall at y=0 perp to y-axis
% D1=absorption coefficient on wall at y=Ly perp to y-axis
% E1=absorption coefficient on wall at z=0 perp to z-axis
% F1=absorption coefficient on wall at z=Lz perp to z-axis
% Rx=microphone position x-coordinate
% Ry=microphone position y-coordinate
% Rz=microphone position z-coordinate
% Sx=speaker position x-coordinate
% Sy=speaker position y-coordinate
% Sz=speaker position z-coordinate

%density of air at 20 degrees
Ro=1.2041;  %Kg/m^3 (from wiki at 20 degrees)- https://en.wikipedia.org/wiki/Acoustic_impedance)

%speed of sound at 20 degrees
C=343.21; %m/s (from wiki at 20 degrees - https://en.wikipedia.org/wiki/Acoustic_impedance)

%room volume
V=Lx*Ly*Lz; 

%Calculate angular frequency
W=2*pi*F;

M=20; %number of modes

%En are scaling factors dependent on the order of the mode

E=ones(1,M)*2;

E(1)=1;

surfx=Ly*Lz*2; %surface area of boundaries perp to x axis
surfy=Lx*Lz*2; %surface area of boundaries perp to y axis
surfz=Lx*Ly*2; %surface area of boundaries perp to z axis

ax=((A1+B1)/2)*surfx;%total surface absorption of the room boundaries perp to x 
ay=((C1+D1)/2)*surfy;%total surface absorption of the room boundaries perp to y
az=((E1+F1)/2)*surfz;%total surface absorption of the room boundaries perp to z

r=sqrt(((Rx-Sx)^2)+((Ry-Sy)^2)+((Rz-Sz)^2)); %radial distance from sound source to receiver

sum=0;
for nx=1:M
for ny=1:M
for nz=1:M   
Kn = (C/(8*V))*((0.5*E(nx)*ax)+(0.5*E(ny)*ay)+(0.5*E(nz)*az));
Wn=pi*C*sqrt((((nx-1)/Lx)^2)+(((ny-1)/Ly)^2)+(((nz-1)/Lz)^2));
Phis=cos((nx-1)*pi*Sx/Lx)*cos((ny-1)*pi*Sy/Ly)*cos((nz-1)*pi*Sz/Lz);
Phir=cos((nx-1)*pi*Rx/Lx)*cos((ny-1)*pi*Ry/Ly)*cos((nz-1)*pi*Rz/Lz);
sum = sum + ((E(nx)*E(ny)*E(nz)*Phis*Phir)/(W^2 - Wn^2 - 2*1i*Kn*Wn));
end
end
end
PR=sum*Ro*(C^2)/V; %mode sum 


u=(PR);


