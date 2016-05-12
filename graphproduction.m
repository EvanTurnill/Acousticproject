%delete figrues between runs

figs2keep = [ ];


all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));

%set default formatting
set(0,'defaultAxesFontSize', 14)

plot(rand(4))
h_legend=legend('One','Two','Three','Four');
set(h_legend,'FontSize',14);

%graph production for all 18 positions


for jj=1:3
    for ii=1:6

%strings for 18 combinations
modalstring=sprintf('MODAL_MX%dSX%d_1to250_X_row1freq_row2amp.mat',ii,jj);
modalvarstring=sprintf('MODALVARDAMP_MX%dSX%d_1to250_X_row1freq_row2amp.mat',ii,jj);
modalvarnodirstring=sprintf('MODALVARDAMPNODIR_MX%dSX%d_1to250_X_row1freq_row2amp.mat',ii,jj);
modalvarwallsstring=sprintf('MODALVARWALLS_MX%dSX%d_1to250_X_row1freq_row2amp.mat',ii,jj);

FDTDstring=sprintf('FDTDMX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);
FDTDvarstring=sprintf('FDTDVARDAMPMX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);
FDTDhighstring=sprintf('FDTDVARHIGHMX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);
FDTDhigh40string=sprintf('FDTDVARHIGH40MX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);
FDTDvarcornerstring=sprintf('FDTDVARCORNERMX%dSX%d_1to250_freq1ZfftR.mat',ii,jj);

Actualstring=sprintf('S%dM%dfrq.csv',jj,ii);

%title for graphs
Titlestring=sprintf('Microphone position %d speaker position %d',ii,jj);

%upload actual data
Y = csvread(Actualstring);

%create vector of all integer freq values for actual data
Yinteger=zeros(2,200);

for i=1:200
val = i; %value to find
tmp = abs(Y(:,1)-val);
[idx idx] = min(tmp); %index of closest value
Yinteger(1,i)=i;
Yinteger(2,i) = Y(idx,2);
end

%label figures
hFig = figure(((jj-1)*6)+ii);
set(hFig, 'Position', [100 100 1200 600])

figure(((jj-1)*6)+ii)

%set axes
xlim manual
ylim manual
hold on
xlim ([12.5,200])
ylim ([0,50])

%plot actual
plot(Y(:,1),Y(:,2),'r')

%load(modalstring);
%plot(X(1,:),(X(2,:)-55))
% Noct = 12;
% Xsmooth = smooth_spectrum((X(2,:)-55),X(1,:),Noct);
% plot(X(1,:),Xsmooth,'y')

%plot modal data
load(modalvarstring);
Noct = 12;
Xsmooth = smooth_spectrum((X(2,:)-55),X(1,:),Noct);
plot(X(1,:),Xsmooth,'g')
% plot(X(1,:)/1.05,(X(2,:)-55))
%plot(X(1,:)/1.0,(X(2,:)-55))
MODAL200=X(:,(1:200));

%find peaks of modal data and compare to modal
ACTUALNORM=norm(Yinteger(2,40:140));
MODALNORM=norm(MODAL200(2,40:140));
ACTUALAV=mean(Yinteger(2,40:140));
MODALAV=mean(MODAL200(2,40:140));

[pks,locs] = findpeaks(Xsmooth(40:140),X(1,40:140));
locsround=round(locs);

max=length(locs);
Ypeaksm=zeros(2,max);
for n=1:max;
Ypeaksm(1,n)=Yinteger(1,locsround(n));
Ypeaksm(2,n)=Yinteger(2, locsround(n));
end
ANVEC=zeros(1,max);
MOVEC=zeros(1,max);
ANVECAV=zeros(1,max);
MOVECAV=zeros(1,max);
for m=1:max
ANVEC(1,m)=ACTUALNORM;
MOVEC(1,m)=MODALNORM;
ANVECAV(1,m)=ACTUALAV;
MOVECAV(1,m)=MODALAV;
end

%quality value for comparison evaluated at peaks of modal(turned out to not
%reflect matching of peaks)
Q1=norm(((Ypeaksm(2,:)-ANVECAV)/ACTUALAV)-((pks-MOVECAV)/MODALAV));

%This is the comparison of all points (also unhelpful)
%Q1=norm((Ypeaksm(2,:)-ANVEC)-(pks-MOVEC));


% load(modalvarwallsstring);
% Noct = 12;
% Xsmooth = smooth_spectrum((X(2,:)-55),X(1,:),Noct);
% plot(X(1,:),Xsmooth,'y')

%load(modalvarnodirstring);
% Noct = 12;
% Xsmooth = smooth_spectrum((X(2,:)-55),X(1,:),Noct);
% plot(X(1,:),Xsmooth)

%plot(X(1,:)/1.05,(X(2,:)-55))

% load(FDTDstring);
% plot(freq1,(R+17.5),'y')
%plot(freq1,(Zfft+17.5))



%load and plot FDTD
load(FDTDvarstring);
plot(freq1,(R+17.5),'b')

%find peaks of FDTD
FDTD200=zeros(2,200);
FDTD15=zeros(2,200);
for k=1:200
FDTD200(1,k)=k;
FDTD200(2,k)=R(5*k); 
FDTD15(1,k)=k*1.1;
FDTD15(2,k)=R(5*k); 

end

ACTUALNORM=norm(Yinteger(2,40:140));
FDTDNORM=norm(FDTD200(2,40:140));
ACTUALAV=mean(Yinteger(2,40:140));
FDTDAV=mean(FDTD200(2,40:140));

[pks,locs] = findpeaks(FDTD200(2,40:140),FDTD200(1,40:140));
locsround=round(locs);

max=length(locs);
Ypeaksf=zeros(2,max);
for n=1:max;
Ypeaksf(1,n)=Y(locsround(n),1);
Ypeaksf(2,n)=Y(locsround(n),2);
end
ANVEC=zeros(1,max);
FDVEC=zeros(1,max);
ANVECAV=zeros(1,max);
FDVECAV=zeros(1,max);
for m=1:max
ANVEC(1,m)=ACTUALNORM;
FDVEC(1,m)=FDTDNORM;
ANVECAV(1,m)=ACTUALAV;
FDVECAV(1,m)=FDTDAV;
end

%Quality number comparing differences at peaks of FDTD and actual (not
%helpful)
Q2=norm(((Ypeaksf(2,:)-ANVECAV)./ACTUALAV)-((pks-FDVECAV)./FDTDAV));

%Quality number comparing actual and FDTD at all points
%Q2=norm((Ypeaksf(2,:)-ANVEC)-(pks-FDVEC));

%comparing the FDTD model with freq multiplied by 15 (decided in the end
%not any better
FDTDNORM=norm(FDTD200(2,40:140));
FDTD15NORM=norm(FDTD15(2,40:140));
ACTUALAV=mean(Yinteger(2,40:140));
FDTD15AV=mean(FDTD15(2,40:140));

[pks,locs] = findpeaks(FDTD15(2,40:140),FDTD15(1,40:140));
locsround=round(locs);

max=length(locs);
Ypeaksf15=zeros(2,max);
for n=1:max;
Ypeaksf15(1,n)=Y((locsround(n)),1);
Ypeaksf15(2,n)=Y((locsround(n)),2);
end
ANVEC=zeros(1,max);
FD15VEC=zeros(1,max);
ANVECAV=zeros(1,max);
FD15VECAV=zeros(1,max);
for m=1:max
ANVEC(1,m)=ACTUALNORM;
FD15VEC(1,m)=FDTD15NORM;
ANVECAV(1,m)=ACTUALAV;
FD15VECAV(1,m)=FDTD15AV;
end

Q3=norm(((Ypeaksf15(2,:)-ANVECAV)./ACTUALAV)-((pks-FD15VECAV)./FDTD15AV));

%Q3=norm((Ypeaksf15(2,:)-ANVEC)-(pks-FD15VEC));

%load(FDTDvarcornerstring);


%plot(1*freq1,(R+17.5),'b')


% load(FDTDhighstring);
% plot(freq1,(R+32.5))

% Noct=12;
% Zfftsmooth = smooth_spectrum((Zfft+17.5),freq1,Noct);
% % plot(X(1,:)/1.07,Xsmooth)
% plot(freq1,(Zfftsmooth+17.5))

% load(FDTDhigh40string);
% plot(freq1,(R+32.5))


% The following code changes the minor grid 
% spacing by adjusting the tick spacing:


% Q1=norm((Yinteger(2,20:200)-ANVEC(20:200))-(MODAL200(2,20:200)-MOVEC(20:200)));
% Q2=norm((Yinteger(2,20:200)-ANVEC(20:200))-(FDTD200(2,20:200)-FDVEC(20:200)));
% Q3=norm((Yinteger(2,20:200)-ANVEC(20:200))-(FDTD15(2,20:200)-FD15VEC(20:200)));
% Q4=norm((MODAL200(2,20:200)-MOVEC(20:200))-(FDTD200(2,20:200)-FDVEC(20:200)));
% Q5=norm((MODAL200(2,20:200)-MOVEC(20:200))-(FDTD15(2,20:200)-FD15VEC(20:200)));

Qstring=sprintf('Q1=%d Q2=%d Q3=%d',Q1,Q2,Q3);

ALLQ=zeros(18,3);
ALLQ((((jj-1)*6)+i),1)=Q1;
ALLQ((((jj-1)*6)+i),2)=Q2;
ALLQ((((jj-1)*6)+i),3)=Q3;

grid on
grid minor
set(gca,'xtick',[0:10:200])
set(gca,'ytick',[0:10:100])
xlabel('frequency, Hz')
ylabel('dB relative to arbitary reference')
title(Titlestring)
%title(Qstring)
legend('Actual measurement','Modal summation method','FDTD','Location','southwest')
%legend('Actual measurement','FDTD','Location','southwest')


hold off 

    end
end
%Show all quality numbers
ALLQ