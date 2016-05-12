function [Z] = FDTDmethod(w, Dims, speaker, mic, A1, B1, C1, D1, E1, F1) 
%Dims = room dimensions
% w=angular frequency
% A1=absorption coefficient on wall at x=0 perp to x-axis
% B1=absorption coefficient on wall at x=Lx perp to x-axis
% C1=absorption coefficient on wall at y=0 perp to y-axis
% D1=absorption coefficient on wall at y=Ly perp to y-axis
% E1=absorption coefficient on wall at z=0 perp to z-axis
% F1=absorption coefficient on wall at z=Lz perp to z-axis
% mic=microphone positions [x,y,z] 
% speaker= speaker positions[x,y,z]


fs=2000;                   % Sample rate (Hz)
lam=sqrt(1/3);             % Courant number
c=343.5;                   % Speed of sound at 20 degrees
simLength = 5;             % Simulation length 
T=1/fs;                    % time step
X=c*T/lam;                 % spatial step

 
DimsN = round(Dims/X)+2;        % Room dimensions [Nx,Ny,Nz] in nodes
speakerN = round(speaker/X);    % Source position in nodes
micN = round(mic/X);            % Microphone position in nodes
maxN = round(simLength/T);      % Sim length in time steps  

t=(0:maxN-1)*T;                       %time steps as a vector      
source=70*cos(t*w);                   %Source at single freq w  

% AllOCATE MATRICES

p = zeros([DimsN 3]);       % Pressure matrix (Nx x Ny x Nz x 3)
                            % Fourth dimension is discrete time
                            % (:,:,:,1) is (n+1)
                            % (:,:,:,2) is (n)
                            % (:,:,:,3) is (n-1)
                           
                         
% Define room length in nodes 
Mx1=DimsN(1);
My1=DimsN(2);
Mz1=DimsN(3);

%Last interior node
Mx=DimsN(1)-1;
My=DimsN(2)-1;
Mz=DimsN(3)-1;


% BOUNDARY CONDITIONS

%specific impedance of each boundary
ImA1=(1+sqrt(1-A1))/(1-sqrt(1-A1));
ImB1=(1+sqrt(1-B1))/(1-sqrt(1-B1));
ImC1=(1+sqrt(1-C1))/(1-sqrt(1-C1));
ImD1=(1+sqrt(1-D1))/(1-sqrt(1-D1));
ImE1=(1+sqrt(1-E1))/(1-sqrt(1-E1));
ImF1=(1+sqrt(1-F1))/(1-sqrt(1-F1));


Z=zeros(2,maxN);

for nn=2:maxN
    
   
%BC

    
%Face1
for j = 2:My
    for k = 2:Mz  
p(1,j,k,1) = (((lam^2)*(2*(p(1+1,j,k,2))+p(1,j+1,k,2)+...
        p(1,j-1,k,2)+p(1,j,k+1,2)+p(1,j,k-1,2)-6*p(1,j,k,2)))+...
        2*p(1,j,k,2)-((1+(1/ImA1*lam))*p(1,j,k,3)))/(1+(1/(ImA1*lam)));  

   end 
end

%Face2 
for j = 2:My
    for k = 2:Mz  
    p(Mx1,j,k,1) = (((lam^2)*(2*(p(Mx,j,k,2))+ ...
        p(Mx1,j+1,k,2)+p(Mx1,j-1,k,2)...
        +p(Mx1,j,k+1,2)+p(Mx1,j,k-1,2)- ...
        6*p(Mx1,j,k,2)))+2*p(Mx1,j,k,2)+ ... 
        ((-1+(1/ImB1*lam))*p(Mx1,j,k,3)))/(1+(1/(ImB1*lam)));

    end 
end

%Face3
for i = 2:Mx
    for k = 2:Mz  
    p(i,1,k,1) = (((lam^2)*((p(i+1,1,k,2)+p(i-1,1,k,2)+(2*p(i,1+1,k,2))...
        +p(i,1,k+1,2)+p(i,1,k-1,2)-6*p(i,1,k,2)))+2*p(i,1,k,2)-...
        ((1+(1/ImC1*lam))*p(i,1,k,3))))/(1+(1/(ImC1*lam)));
    end 
end

%Face4
for i = 2:Mx
    for k = 2:Mz  
    p(i,My1,k,1) = (((lam^2)*((p(i+1,My1,k,2)+ ...
        p(i-1,My1,k,2)+(2*p(i,My,k,2))...
        +p(i,My1,k+1,2)+p(i,My1,k-1,2)- ...
        6*p(i,My1,k,2)))+2*p(i,My1,k,2)+ ...
        ((-1+(1/ImD1*lam))*p(i,My1,k,3))))/(1+(1/(ImD1*lam)));

    end 
end

%Face5
for i = 2:Mx
    for j = 2:My  
    p(i,j,1,1) = (((lam^2)*((p(i+1,j,1,2)+p(i-1,j,1,2)+p(i,j+1,1,2)...
        +p(i,j-1,1,2)+(2*(p(i,j,1+1,2)))-6*p(i,j,1,2)))+2*p(i,j,1,2)-...
        ((1+(1/ImE1*lam))*p(i,j,1,3))))/(1+(1/(ImE1*lam)));
    end 
end

%Face6
for i = 2:Mx
    for j = 2:My  
    p(i,j,1,Mz1) = (((lam^2)*((p(i+1,j,Mz1,2)+ ...
        p(i-1,j,Mz1,2)+p(i,j+1,Mz1,2)...
        +p(i,j-1,Mz1,2)+(2*(p(i,j,Mz,2))-...
        6*p(i,j,Mz1,2))))+2*p(i,j,Mz1,2)+...
        ((-1+(1/ImF1*lam))*p(i,j,Mz1,3))))/(1+(1/(ImF1*lam)));
    end 
end

% %corners
% 
%(0,0,0)
p(1,1,1,1) = (((lam^2)*((2*(p(1+1,1,1,2)))+(2*(p(1,1+1,1,2)))+...
        (2*p(1,1,1+1,2))-(6*p(1,1,1,2))))+...
        (2*p(1,1,1,2))-((1+(1/ImA1*lam)+(1/ImC1*lam)+(1/ImE1*lam))*p(1,1,1,3)))/(1+(1/(ImA1*lam)+(1/ImC1*lam)+(1/ImE1*lam)));  

%(Mx1,0,0)
p(Mx1,1,1,1) = (((lam^2)*((2*(p(Mx,1,1,2)))+(2*(p(Mx1,1+1,1,2)))+...
        (2*p(Mx1,1,1+1,2))-(6*p(Mx1,1,1,2))))+...
        (2*p(Mx1,1,1,2))-((1-(1/ImB1*lam)+(1/ImC1*lam)+(1/ImE1*lam))*p(Mx1,1,1,3)))/(1+(1/(ImB1*lam)+(1/ImC1*lam)+(1/ImE1*lam)));  

    
%(0,My1,0)
p(1,My1,1,1) = (((lam^2)*((2*(p(1+1,My1,1,2)))+(2*(p(1,My,1,2)))+...
        (2*p(1,My1,1+1,2))-(6*p(1,My1,1,2))))+...
        (2*p(1,My1,1,2))-((1+(1/ImA1*lam)-(1/ImD1*lam)+(1/ImE1*lam))*p(1,My1,1,3)))/(1+(1/(ImA1*lam)+(1/ImD1*lam)+(1/ImE1*lam)));  
     
%(0,0,Mz1)
p(1,1,Mz1,1) = (((lam^2)*((2*(p(1+1,1,Mz1,2)))+(2*(p(1,1+1,Mz1,2)))+...
        (2*p(1,1,Mz,2))-(6*p(1,1,Mz1,2))))+...
        (2*p(1,1,Mz1,2))-((1+(1/ImA1*lam)+(1/ImC1*lam)-(1/ImF1*lam))*p(1,1,Mz1,3)))/(1+(1/(ImA1*lam)+(1/ImC1*lam)+(1/ImF1*lam)));  

%(Mx1,My1,0)
p(Mx1,My1,1,1) = (((lam^2)*((2*(p(Mx,My1,1,2)))+(2*(p(Mx1,My,1,2)))+...
        (2*p(Mx1,My1,1+1,2))-(6*p(Mx1,My1,1,2))))+...
        (2*p(Mx1,My1,1,2))-((1-(1/ImB1*lam)-(1/ImD1*lam)+(1/ImE1*lam))*p(Mx1,My1,1,3)))/(1+(1/(ImB1*lam)+(1/ImD1*lam)+(1/ImE1*lam))); 
    
%(Mx1,0,Mz1)
p(Mx1,1,Mz1,1) = (((lam^2)*((2*(p(Mx,1,Mz1,2)))+(2*(p(Mx1,1+1,Mz1,2)))+...
        (2*p(Mx1,1,Mz,2))-(6*p(Mx1,1,Mz1,2))))+...
        (2*p(Mx1,1,Mz1,2))-((1-(1/ImB1*lam)+(1/ImC1*lam)-(1/ImF1*lam))*p(Mx1,1,Mz1,3)))/(1+(1/(ImB1*lam)+(1/ImC1*lam)+(1/ImF1*lam)));  

%(0,My1,Mz1)
p(1,My1,Mz1,1) = (((lam^2)*((2*(p(1+1,My1,Mz1,2)))+(2*(p(1,My,Mz1,2)))+...
        (2*p(1,My1,Mz,2))-(6*p(1,My1,Mz1,2))))+...
        (2*p(1,My1,Mz1,2))-((1+(1/ImA1*lam)-(1/ImD1*lam)-(1/ImF1*lam))*p(1,My1,Mz1,3)))/(1+(1/(ImA1*lam)+(1/ImD1*lam)+(1/ImF1*lam)));  
    
%(Mx1,My1,Mz1)
p(Mx1,My1,Mz1,1) = (((lam^2)*((2*(p(Mx,My1,Mz1,2)))+(2*(p(Mx1,My,Mz1,2)))+...
        (2*p(Mx1,My1,Mz,2))-(6*p(Mx1,My1,Mz1,2))))+...
        (2*p(Mx1,My1,Mz1,2))-((1-(1/ImB1*lam)-(1/ImD1*lam)-(1/ImF1*lam))*p(Mx1,My1,Mz1,3)))/(1+(1/(ImB1*lam)+(1/ImD1*lam)+(1/ImF1*lam)));  


%Edges

%(0,0,0) to (Mx1,0,0)

for i=2:Mx
p(i,1,1,1) = (((lam^2)*...
        (p(i+1,1,1,2)+p(i-1,1,1,2)+(2*(p(i,1+1,1,2)))+(2*p(i,1,1+1,2))-(6*p(i,1,1,2)))+...
        (2*p(i,1,1,2))-((1+(1/ImC1*lam)+(1/ImE1*lam))*p(i,1,1,3)))/...
        (1+(1/ImC1*lam)+(1/ImE1*lam)));   
    
end

%(0,My1,0) to (Mx1,My1,0)
for i=2:Mx
p(i,My1,1,1) = (((lam^2)*...
        (p(i+1,My1,1,2)+p(i-1,My1,1,2)+(2*(p(i,My,1,2)))+(2*p(i,My1,1+1,2))-(6*p(i,My1,1,2)))+...
        (2*p(i,My1,1,2))-((1-(1/ImD1*lam)+(1/ImE1*lam))*p(i,My1,1,3)))/...
        (1+(1/ImD1*lam)+(1/ImE1*lam)));   
    
end
%(0,0,Mz1) to (Mx1,0,Mz1)
for i=2:Mx
p(i,1,1,Mz1) = (((lam^2)*...
        (p(i+1,1,Mz1,2)+p(i-1,1,Mz1,2)+(2*(p(i,1+1,Mz1,2)))+(2*p(i,1,Mz,2))-(6*p(i,1,Mz1,2)))+...
        (2*p(i,1,Mz1,2))-((1+(1/ImC1*lam)-(1/ImF1*lam))*p(i,1,1,3)))/...
        (1+(1/ImC1*lam)+(1/ImF1*lam)));   
    
end


%(0,My1,Mz1) to (Mx1,My1,Mz1)

for i=2:Mx
p(i,My1,Mz1,1) = (((lam^2)*...
        (p(i+1,My1,Mz1,2)+p(i-1,My1,Mz1,2)+(2*(p(i,My,Mz1,2)))+(2*p(i,My1,Mz,2))-(6*p(i,My1,Mz1,2)))+...
        (2*p(i,My1,Mz1,2))-((1+(1/ImD1*lam)+(1/ImF1*lam))*p(i,1,1,3)))/...
        (1-(1/ImD1*lam)-(1/ImF1*lam)));   
end



%(0,0,0)to (0,My1,0)
for j=2:My
p(1,j,1,1) = (((lam^2)*...
        (2*p(1+1,j,1,2)+(p(1,j+1,1,2))+(p(1,j-1,1,2))+(2*p(1,j,1+1,2))-(6*p(1,j,1,2)))+...
        (2*p(1,j,1,2))-((1+(1/ImA1*lam)+(1/ImE1*lam))*p(1,j,1,3)))/...
        (1+(1/ImA1*lam)+(1/ImE1*lam)));   
end

%(Mx1,0,0) to (Mx1,My1,0)
for j=2:My
p(Mx1,j,1,1) = (((lam^2)*...
        (2*p(Mx,j,1,2)+(p(Mx1,j+1,1,2))+(p(Mx1,j-1,1,2))+(2*p(Mx1,j,1+1,2))-(6*p(Mx1,j,1,2)))+...
        (2*p(Mx1,j,1,2))-((1-(1/ImB1*lam)+(1/ImE1*lam))*p(Mx1,j,1,3)))/...
        (1+(1/ImB1*lam)+(1/ImE1*lam)));   
end


%(0,0,Mz1) to (0,My1,Mz1)
for j=2:My
p(1,j,Mz1,1) = (((lam^2)*...
        (2*p(1+1,j,Mz1,2)+(p(1,j+1,Mz1,2))+(p(1,j-1,Mz1,2))+(2*p(1,j,Mz,2))-(6*p(1,j,Mz1,2)))+...
        (2*p(1,j,Mz1,2))-((1+(1/ImA1*lam)-(1/ImF1*lam))*p(1,j,Mz1,3)))/...
        (1+(1/ImA1*lam)+(1/ImE1*lam)));   
end

%(Mx1,0,Mz1) to (Mx1,My1,Mz1)
for j=2:My
p(Mx1,j,Mz1,1) = (((lam^2)*...
        (2*p(Mx,j,Mz1,2)+(p(Mx1,j+1,Mz1,2))+(p(Mx1,j-1,Mz1,2))+(2*p(Mx1,j,Mz,2))-(6*p(Mx1,j,Mz1,2)))+...
        (2*p(Mx1,j,Mz1,2))-((1-(1/ImB1*lam)-(1/ImF1*lam))*p(1,j,Mz1,3)))/...
        (1+(1/ImA1*lam)+(1/ImE1*lam)));   
end


%(Mx1,My1,Mz1) to (Mx1,0,Mz1)
for j=2:My
p(Mx1,j,Mz1,1) = (((lam^2)*...
        (2*p(Mx,j,Mz1,2)+(p(Mx1,j+1,Mz1,2))+(p(Mx1,j-1,Mz1,2))+(2*p(Mx1,j,Mz,2))-(6*p(Mx1,j,Mz1,2)))+...
        (2*p(Mx1,j,Mz1,2))-((1-(1/ImB1*lam)-(1/ImF1*lam))*p(1,j,1,3)))/...
        (1+(1/ImB1*lam)+(1/ImF1*lam)));   
    
end

%(0,0,0)to (0,0,Mz1)
for k=2:Mz
p(1,1,k,1) = (((lam^2)*...
        (2*p(1+1,1,k,2)+(2*(p(1,1+1,k,2)))+p(1,1,k+1,2)+p(1,1,k-1,2)-(6*p(1,1,k,2)))+...
        (2*p(1,1,k,2))-((1+(1/ImA1*lam)+(1/ImC1*lam))*p(1,1,k,3)))/...
        (1+(1/ImA1*lam)+(1/ImC1*lam)));   
end

%(Mx1,0,0) to (Mx1,0,Mz1)
for k=2:Mz
p(Mx1,1,k,1) = (((lam^2)*...
        (2*p(Mx,1,k,2)+(2*(p(Mx1,1+1,k,2)))+p(Mx1,1,k+1,2)+p(Mx1,1,k-1,2)-(6*p(Mx1,1,k,2)))+...
        (2*p(Mx1,1,k,2))-((1-(1/ImB1*lam)+(1/ImC1*lam))*p(1,1,k,3)))/...
        (1+(1/ImB1*lam)+(1/ImC1*lam)));   
end

%(0,My1,0) to (0,0,Mz1)
for k=2:Mz
p(1,My1,k,1) = (((lam^2)*...
        (2*p(1+1,My1,k,2)+(2*(p(1,My,k,2)))+p(1,My1,k+1,2)+p(1,My1,k-1,2)-(6*p(1,My1,k,2)))+...
        (2*p(1,My1,k,2))-((1+(1/ImA1*lam)-(1/ImD1*lam))*p(1,My1,k,3)))/...
        (1+(1/ImA1*lam)+(1/ImC1*lam)));   
end

%(Mx1,My1,Mz1) to (Mx1,My1,0)
for k=2:Mz
p(Mx1,My1,k,1) = (((lam^2)*...
        (2*p(Mx,My1,k,2)+(2*(p(Mx1,My,k,2)))+p(Mx1,My1,k+1,2)+p(Mx1,My1,k-1,2)-(6*p(Mx1,My1,k,2)))+...
        (2*p(Mx1,My1,k,2))-((1-(1/ImB1*lam)-(1/ImD1*lam))*p(Mx1,My1,k,3)))/...
        (1+(1/ImA1*lam)+(1/ImC1*lam)));  
end
    
    
    % Update the grid
    for i=2:Mx
        for j=2:My
            for k=2:Mz
        p(i,j,k,1) = ((lam^2)*((p(i+1,j,k,2)+p(i-1,j,k,2)+p(i,j+1,k,2)+...
        p(i,j-1,k,2)+p(i,j,k+1,2)+p(i,j,k-1,2)-6*p(i,j,k,2)))+...
        2*p(i,j,k,2)-p(i,j,k,3));  
            end
        end
    end
    

    % Add the source
    p(speakerN(1),speakerN(2),speakerN(3),2) = p(speakerN(1),speakerN(2),speakerN(3),2) + source(nn);
    %p(srcPosD2(1),srcPosD2(2),srcPosD2(3),2) = p(srcPosD2(1),srcPosD2(2),srcPosD2(3),2) + srcFn2(nn);
    
    % Take a time step
    p(:,:,:,3) = p(:,:,:,2);
    p(:,:,:,2) = p(:,:,:,1);
    

    Z(1,nn)=nn;
    Z(2,nn)=(p(micN(1),micN(2),micN(3),1));
end


