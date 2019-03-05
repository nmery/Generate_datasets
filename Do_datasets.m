%% Inputs
load('mayo_data_for_denis.mat');
cap=30;     % capping value
lc=3;       % lenght composites
dx=5;       % block model size x
dy=5;       % block model size y
dz=5;       % block model size z

%% Statistics raw data
stats_ddh_raw=[mean(ddh(:,4)), var(ddh(:,4)),std(ddh(:,4))]
ststs_rc_raw=[mean(rc(:,4)) var(rc(:,4)),std(rc(:,4))]

%% Composites
[ddh3]=compos_rosebel(ddh,lc);

%% Capping
for i=1:length(rc)
    if rc(i,4)>=cap
        rc(i,4)=cap;
    else
    end
end

for i=1:length(ddh3)
    if ddh3(i,4)>=cap
        ddh3(i,4)=cap;
    else
    end
end

save('rc_cap.dat','rc','-ascii')

% Statistics cap RC
ststs_rc_cap=[mean(rc(:,4)) var(rc(:,4)),std(rc(:,4))]

%% Define RC domain
xmin=floor(min(rc(:,1)));
xmax=round(max(rc(:,1)));
ymin=floor(min(rc(:,2)));
ymax=round(max(rc(:,2)));
zmin=floor(min(rc(:,3)));
zmax=round(max(rc(:,3)));

BM_rc=grille3(xmin-dx/2,xmax,dx,ymin-dy/2,ymax,dy,zmin-dz/2,zmax,dz);

[x0s_rc,~,~,~]=invd_new_3d(rc,BM_rc,1,6.3,20);
x0s_rc(any(isnan(x0s_rc),2),:)=[];

%% Define DDH dataset (according to RC domain)
xx=invd_new_3d(x0s_rc,ddh3(:,1:3),1,20,1);
ddh3=ddh3(~isnan(xx(:,end)),:);
save('ddh.dat','ddh3','-ascii')

% Statistics DDH dataset
stats_ddh3_cap=[mean(ddh3(:,4)), var(ddh3(:,4)),std(ddh3(:,4))]

%% Turning bands simulation (RC reality)

% Add block discretization (3x3x3)
x0s_d=[];
for i=1:length(x0s_rc)
    x0s_d1=[x0s_rc(i,1)-(5/3),x0s_rc(i,2)-(5/3),x0s_rc(i,3)-(5/3)];
    x0s_d2=[x0s_rc(i,1),x0s_rc(i,2)-(5/3),x0s_rc(i,3)-(5/3)];
    x0s_d3=[x0s_rc(i,1)+(5/3),x0s_rc(i,2)-(5/3),x0s_rc(i,3)-(5/3)];
    x0s_d4=[x0s_rc(i,1)-(5/3),x0s_rc(i,2),x0s_rc(i,3)-(5/3)];
    x0s_d5=[x0s_rc(i,1),x0s_rc(i,2),x0s_rc(i,3)-(5/3)];
    x0s_d6=[x0s_rc(i,1)+(5/3),x0s_rc(i,2),x0s_rc(i,3)-(5/3)];
    x0s_d7=[x0s_rc(i,1)-(5/3),x0s_rc(i,2)+(5/3),x0s_rc(i,3)-(5/3)];
    x0s_d8=[x0s_rc(i,1),x0s_rc(i,2)+(5/3),x0s_rc(i,3)-(5/3)];
    x0s_d9=[x0s_rc(i,1)+(5/3),x0s_rc(i,2)+(5/3),x0s_rc(i,3)-(5/3)];
    x0s_d10=[x0s_rc(i,1)-(5/3),x0s_rc(i,2)-(5/3),x0s_rc(i,3)];
    x0s_d11=[x0s_rc(i,1),x0s_rc(i,2)-(5/3),x0s_rc(i,3)];
    x0s_d12=[x0s_rc(i,1)+(5/3),x0s_rc(i,2)-(5/3),x0s_rc(i,3)];
    x0s_d13=[x0s_rc(i,1)-(5/3),x0s_rc(i,2),x0s_rc(i,3)];
    x0s_d14=[x0s_rc(i,1),x0s_rc(i,2),x0s_rc(i,3)];
    x0s_d15=[x0s_rc(i,1)+(5/3),x0s_rc(i,2),x0s_rc(i,3)];
    x0s_d16=[x0s_rc(i,1)-(5/3),x0s_rc(i,2)+(5/3),x0s_rc(i,3)];
    x0s_d17=[x0s_rc(i,1),x0s_rc(i,2)+(5/3),x0s_rc(i,3)];
    x0s_d18=[x0s_rc(i,1)+(5/3),x0s_rc(i,2)+(5/3),x0s_rc(i,3)];
    x0s_d19=[x0s_rc(i,1)-(5/3),x0s_rc(i,2)-(5/3),x0s_rc(i,3)+(5/3)];
    x0s_d20=[x0s_rc(i,1),x0s_rc(i,2)-(5/3),x0s_rc(i,3)+(5/3)];
    x0s_d21=[x0s_rc(i,1)+(5/3),x0s_rc(i,2)-(5/3),x0s_rc(i,3)+(5/3)];
    x0s_d22=[x0s_rc(i,1)-(5/3),x0s_rc(i,2),x0s_rc(i,3)+(5/3)];
    x0s_d23=[x0s_rc(i,1),x0s_rc(i,2),x0s_rc(i,3)+(5/3)];
    x0s_d24=[x0s_rc(i,1)+(5/3),x0s_rc(i,2),x0s_rc(i,3)+(5/3)];
    x0s_d25=[x0s_rc(i,1)-(5/3),x0s_rc(i,2)+(5/3),x0s_rc(i,3)+(5/3)];
    x0s_d26=[x0s_rc(i,1),x0s_rc(i,2)+(5/3),x0s_rc(i,3)+(5/3)];
    x0s_d27=[x0s_rc(i,1)+(5/3),x0s_rc(i,2)+(5/3),x0s_rc(i,3)+(5/3)];
    x0s_d=[x0s_d;x0s_d1;x0s_d2;x0s_d3;x0s_d4;x0s_d5;x0s_d6;x0s_d7;x0s_d8;x0s_d9;x0s_d10;x0s_d11;x0s_d12;x0s_d13;x0s_d14;
        x0s_d15;x0s_d16;x0s_d17;x0s_d18;x0s_d19;x0s_d20;x0s_d21;x0s_d22;x0s_d23;x0s_d24;x0s_d25;x0s_d26;x0s_d27];
end

% Obtain anamorphosis of rc
y=anamor(rc(:,4));

% Apply TB
x0=x0s_d;
x=[rc(:,1:3),y(:,2)];
model=[1 1; 4 25; 4 60];
c=[0.49;0.46;0.05];
datasim=tourband_cond_large(x,x0,model,c,123456,1);

% Back-transforming TB
yz=[rc(:,1:3),y];
[z]=anamorinv(yz,datasim(:,4));
za=z(:,2);

% Average of discretization
zf=[];
for ii=1:27:length(za)
    zf=[zf;mean(za(ii:ii+26,:))];
end
real=[x0s_rc(:,1:3),zf];
save('real_tb.dat','real','-ascii')

% Statistics of real (TB)
stats_real_tb=[mean(real(:,4)),var(real(:,4)),std(real(:,4))]

%% Define panel model UC

% Create grid panel model
xmin=floor(min(real(:,1)));
xmax=round(max(real(:,1)));
dxp=50;
ymin=floor(min(real(:,2)));
ymax=round(max(real(:,2)));
dyp=50;
zmin=floor(min(real(:,3)));
zmax=round(max(real(:,3)));
dzp=30;
BM_uc=grille3(xmin-dx/2+dxp/2,xmax+dx/2,dxp,ymin-dy/2+dyp/2,ymax+dy/2,dyp,zmin-dz/2+dzp/2,zmax+dz/2,dzp);

% Remove panels far to RC domain
[x0s_uc,~,~,~]=invd_new_3d(x0s_rc,BM_uc,1,40,1);
x0s_uc(any(isnan(x0s_uc),2),:)=[];
bm_uc=x0s_uc(:,1:3);
save('BM_uc.dat','bm_uc','-ascii')

%% Declustered dataset

% Define declustered fine grid
xmin=floor(min(real(:,1)));
xmax=round(max(real(:,1)));
dxd=2;
ymin=floor(min(real(:,2)));
ymax=round(max(real(:,2)));
dyd=2;
zmin=floor(min(real(:,3)));
zmax=round(max(real(:,3)));
dzd=2;
BM_fine=grille3(xmin-dx/2+dxd/2,xmax+dx/2,dxd,ymin-dy/2+dyd/2,ymax+dy/2,dyd,zmin-dz/2+dzd/2,zmax+dz/2,dzd);

% Apply NN to obtain declustered data
[x0s_dec,~,~]=invd_new_3d(ddh3,BM_fine,1,600,1);

% Remove blocks far to RC domain
[x0s_dec1,~,~,~]=invd_new_3d(x0s_rc,x0s_dec(:,1:3),1,7,1);
x0s_dec2=x0s_dec(~isnan(x0s_dec1(:,end)),:);

save('ddh_dec.dat','x0s_dec2','-ascii')

% Stats declustered ddh
stats_ddh_dec=[mean(x0s_dec2(:,4)), var(x0s_dec2(:,4)),std(x0s_dec2(:,4))]

%% Define DDH given the RC samples (ranking approach)
[n,p]=size(ddh3);
vq=[[0:0.005:99.9]/100]';
[y,ii]=sort(ddh3(:,4));
r=[1:n]';
rank=r/(n+1);
z=quantile(rc(:,end),vq);
w=interp1(vq,z,rank,'linear',z(end,1));
ddh3(ii,4)=w;

% Statistics ddh based on rc
stats_ddh_rc=[mean(ddh3(:,4)), var(ddh3(:,4)),std(ddh3(:,4))]

save('ddh_rc.dat','ddh3','-ascii')
