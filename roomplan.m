%room plan

Lx=2.950;
Ly=6.105;
Lz=2.260;

Mcom=zeros(6,3);
Mcom(1,:)=M1;
Mcom(2,:)=M2;
Mcom(3,:)=M3;
Mcom(4,:)=M4;
Mcom(5,:)=M5;
Mcom(6,:)=M6;

M1=[1.475,2.030,1.13];
M2=[0.725,2.03,1.13];
M3=[2.225,2.03,1.65];
M4=[0.725,3.03,1.65];
M5=[1.475,3.03,1.65];
M6=[2.225,3.03,1.13];


S1=[0.7,0.45,1.33];
S2=[1.730,0.45,1.33];
S3=[1.2,1.3,1.03];

Scom=zeros(3,3);
Scom(1,:)=S1;
Scom(2,:)=S2;
Scom(3,:)=S3;

xmin=0;
xmax=Lx;
ymin=0;
ymax=Ly;
zmin=0;
zmax=Lz;


clf;
figure(1);
format compact 
h(1) = axes('Position',[0.2 0.2 0.6 0.6]);
vert = [xmax ymax zmin; 
        xmin ymax zmin; 
        xmin ymax zmax; 
        xmax ymax zmax; 
        xmin ymin zmax;
        xmax ymin zmax; 
        xmax ymin zmin;
        xmin ymin zmin];


fac = [1 2 3 4; 
       4 3 5 6; 
       6 7 8 5; 
       1 2 8 7; 
       6 7 1 4; 
       2 3 5 8];

   
xmin=0;
xmax=Lx;
ymin=1;
ymax=1.3;
zmin=Lz-0.2;
zmax=Lz;

% 
vert2 = [xmax ymax zmin; 
        xmin ymax zmin; 
        xmin ymax zmax; 
        xmax ymax zmax; 
        xmin ymin zmax;
        xmax ymin zmax; 
        xmax ymin zmin;
        xmin ymin zmin];

fac2 = fac;

xmin=Lx-2;
xmax=Lx;
ymin=Ly-1;
ymax=Ly;
zmin=0;
zmax=.8;

% 
vert3 = [xmax ymax zmin; 
        xmin ymax zmin; 
        xmin ymax zmax; 
        xmax ymax zmax; 
        xmin ymin zmax;
        xmax ymin zmax; 
        xmax ymin zmin;
        xmin ymin zmin];

fac3 = fac;

xmin=0;
xmax=0.2;
ymin=(0.5*Ly);
ymax=(0.5*Ly)+1;
zmin=0;
zmax=1;

% 
vert4 = [xmax ymax zmin; 
        xmin ymax zmin; 
        xmin ymax zmax; 
        xmax ymax zmax; 
        xmin ymin zmax;
        xmax ymin zmax; 
        xmax ymin zmin;
        xmin ymin zmin];

fac4 = fac;

xmin=0;
xmax=1.5;
ymin=0;
ymax=0.6;
zmin=0;
zmax=0.9;


vert5 = [xmax ymax zmin; 
        xmin ymax zmin; 
        xmin ymax zmax; 
        xmax ymax zmax; 
        xmin ymin zmax;
        xmax ymin zmax; 
        xmax ymin zmin;
        xmin ymin zmin];

fac5 = fac;

patch('Faces',fac,'Vertices',vert,'FaceColor','none');  % patch function
axis([-4, 4, -4, 4, -4, 4]);
axis equal;
grid on;
grid minor
%set(gca,'xtick',[0:0.6:Lx])
set(gca,'XTickLabel','')
set(gca,'ytick',[0:0.6:Ly])
set(gca,'ztick',[0:0.6:Lz])

xlabel('x','FontSize',16) % x-axis label
ylabel('y','FontSize',16) % y-axis label
zlabel('z','FontSize',16) % y-axis label

hold on;

patch('Faces', fac2, 'Vertices', vert2, 'FaceColor', [0,0,0]);

patch('Faces', fac3, 'Vertices', vert3, 'FaceColor', [0,0,0]);

patch('Faces', fac4, 'Vertices', vert4, 'FaceColor', [0,0,0]);

patch('Faces', fac5, 'Vertices', vert5, 'FaceColor', [0,0,0]);

micstring=sprintf('Mic%d',i);
sourcestring=sprintf('Source%d',j);


txt1=('\otimes mic1');
txt2=('\otimes mic2');
txt3=('\otimes mic3');
txt4=('\otimes mic4');
txt5=('\otimes mic5');
txt6=('\otimes mic6');

txt7=('\oplus source1');
txt8=('\oplus source2');
txt9=('\oplus source3');

x1=zeros(1,6);
y1=zeros(1,6);
z1=zeros(1,6);


x2=zeros(1,3);
y2=zeros(1,3);
z2=zeros(1,3);

for i=1:6
x1(1,i) = Mcom(i,1);
y1(1,i) = Mcom(i,2);
z1(1,i) = Mcom(i,3);
end



tx1=text(x1(1),y1(1),z1(1), txt1);
tx1.Color = [1 0 0]; 
tx1.FontSize = 9;       
tx1.FontWeight = 'bold'; 
tx2=text(x1(2),y1(2),z1(2), txt2)
tx2.Color = [1 1 0]; 
tx2.FontSize = 9;       
tx2.FontWeight = 'bold'; 
tx3=text(x1(3),y1(3),z1(3), txt3)
tx3.Color = [0 1 1]; 
tx3.FontSize = 9;       
tx3.FontWeight = 'bold'; 
tx4=text(x1(4),y1(4),z1(4), txt4)
tx4.Color = [0 1 0]; 
tx4.FontSize = 9;       
tx4.FontWeight = 'bold'; 
tx5=text(x1(5),y1(5),z1(5), txt5)
tx5.Color = [0 0 1]; 
tx5.FontSize = 9;       
tx5.FontWeight = 'bold'; 
tx6=text(x1(6),y1(6),z1(6), txt6)
tx6.Color = [1 0 1]; 
tx6.FontSize = 9;       
tx6.FontWeight = 'bold'; 

for j=1:3
x2(1,j) = Scom(j,1);
y2(1,j) = Scom(j,2);
z2(1,j) = Scom(j,3);
end


tx7=text(x2(1),y2(1),z2(1), txt7);
tx7.Color = [0 0 0]; 
tx7.FontSize = 9;       
tx7.FontWeight = 'bold'; 
tx8=text(x2(2),y2(2),z2(2), txt8);
tx8.Color = [1 0.6 0]; 
tx8.FontSize = 9;       
tx8.FontWeight = 'bold'; 
tx9=text(x2(3),y2(3),z2(3), txt9);
tx9.Color = [0.6 1 0]; 
tx9.FontSize = 9;       
tx9.FontWeight = 'bold'; 



material metal;
alpha('color');
alphamap('rampdown');
view(3);