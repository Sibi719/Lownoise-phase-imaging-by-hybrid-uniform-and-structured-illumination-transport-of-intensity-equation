clc;
clear all;
close all;
N=1024;
wo=N/2;
dpix=4.4*10^-6;
lambda=0.0633*10^-9;
k=2*pi/lambda;
d=0.4;
x=(-N/2:N/2-1)*dpix;
y=(-N/2:N/2-1)*dpix;
[X,Y]=meshgrid(x,y);
umax=1/(2*dpix); 
df=1/(N*dpix);
u=(-N/2:N/2-1)*df;
v=(-N/2:N/2-1)*df;
[U,V]=meshgrid(u,v);
freq_sq=U.^2+V.^2;
nn=1;
vari=10^-4;
noise=randn(N).*sqrt(vari);

fm=14.2*10^3;
m=0.5;
a=0.5;

Io=ones(N);
Imx=a + 0.5*m*sin(2*pi*X.*fm);
Iqx=a + 0.5*m*cos(2*pi*X.*fm);
Imy=a + 0.5*m*sin(2*pi*Y.*fm);
Iqy=a + 0.5*m*cos(2*pi*Y.*fm);

W=zeros(N);
W= sqrt(U.^2 + V.^2)<2*fm;

figure
imagesc(W);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['W']);


P = imread('image2.jpg');
P=rgb2gray(P);
P = double(P);
P=rescale(P,-2.6,2.6);

figure
imagesc(x*10^3,y*10^3,P);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['Input Phase profile']);

Object=exp(1j.*P);

Io=ones(N,N);
Id=FresProp(dpix,d,lambda,N,Object);

figure
imagesc(x*10^3,y*10^3,Id);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['Id']);

[UP, deld]=uniform(dpix,d,lambda,N,Object,k,freq_sq,noise,Id,Io,nn);

figure
imagesc(x*10^3,y*10^3,UP);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['Recovered phase using uniform illumination']);


Imxd=FresProp(dpix,d,lambda,N,Object.*sqrt(Imx));
Iqxd=FresProp(dpix,d,lambda,N,Object.*sqrt(Iqx));
Imyd=FresProp(dpix,d,lambda,N,Object.*sqrt(Imy));
Iqyd=FresProp(dpix,d,lambda,N,Object.*sqrt(Iqy));

[phi_x]=phix(k,d,m,fm,Id,Imx,Iqx,Imxd,Iqxd,X,Y,W);

figure
imagesc(x*10^3,y*10^3,real(phi_x));
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['phi_x']);


[phi_y]=phiy(k,d,m,fm,Id,Imy,Iqy,Imyd,Iqyd,X,Y,W);

figure
imagesc(x*10^3,y*10^3,real(phi_y));
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['phi_y']);

phi_x= real(phi_x);
phi_y=real(phi_y);
[Lx,~]=gradient(phi_x,dpix,dpix);
[~,Ly]=gradient(phi_y,dpix,dpix);

L=Lx+Ly;
T1=fftshift(fft2(L));
T2= -(4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
SP= ifft2(ifftshift(T2.*T1));


figure
imagesc(x*10^3,y*10^3,SP);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['SP']);

%%Hybrid

T= (U.^2 + V.^2)./((fm*0.5*m)^2 + (U.^2 + V.^2));
figure
mesh(T);
colorbar
colormap(jet);
axis on
title(['filter']);

F_HP=fftshift(fft2(UP)).*T + fftshift(fft2(SP)).*(1-T);
HP=ifft2(ifftshift(F_HP));
figure
imagesc(x*10^3,y*10^3,HP);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['HP']);


function[phi_x]=phix(k,d,m,fm,Id,Imx,Iqx,Imxd,Iqxd,X,Y,W)
Smx=(-k./(d*m*pi*fm)).*(Imxd - (Id.*Imx));
Sqx=((1j.*k)./(d*m*pi*fm)).*(Iqxd - (Id.*Iqx));
Smx=Smx.*exp(-1j*2*pi*fm.*X);
Sqx=Sqx.*exp(-1j*2*pi*fm.*X);
F_Smx=fftshift(fft2(Smx));
F_Sqx=fftshift(fft2(Sqx));
F_phi_x=(F_Smx+F_Sqx).*W;
phi_x=ifft2(ifftshift(F_phi_x));
end
function[phi_y]=phiy(k,d,m,fm,Id,Imy,Iqy,Imyd,Iqyd,X,Y,W)
Smy=(-k./(d*m*pi*fm)).*(Imyd - (Id.*Imy));
Sqy=((1j.*k)./(d*m*pi*fm)).*(Iqyd - (Id.*Iqy));
Smy=Smy.*exp(-1j*2*pi*fm.*Y);
Sqx=Sqy.*exp(-1j*2*pi*fm.*Y);
F_Smy=fftshift(fft2(Smy));
F_Sqy=fftshift(fft2(Sqy));
F_phi_y=(F_Smy+F_Sqy).*W;
phi_y=ifft2(ifftshift(F_phi_y));
end
function [uni_rec_ph, deld]=uniform(dpix,d,lambda,N,Object,k,freq_sq,noise,Id,Io,nn)
deld=(Id-Io)+nn.*noise;
T1=fftshift((fft2((k/d).*(deld))));
T2= (4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
uni_rec_ph= ifft2(ifftshift(T2.*T1));
end
function [R]=FresProp(dpix,z,lambda,Hsize,Hcrop)

%Spatial frequencies
Xsize = Hsize*dpix; %Hsize is the number of pixel, dpix is the length of one pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
%Nyquist cut-off for Sampling Hologram
umax = 1/(2*dpix); %define the k space 
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
clear  u V  du;
 
%Evanescent cut-off 
uev = 1/lambda; %???????
 
%Nyquist cut-off for Fresnel Propagation Kernel
unp = uev*(Xsize/(2*abs(z)));
clear Xsize;
 
%Circular window
A = U.^2+(U').^2;
clear U;
if uev>=unp
    ucut = unp;
end
if unp>uev
    ucut = uev;
end
W= sqrt(A);
W = (W<=ucut); 
% disp(['Cutoff =',num2str(ucut),' Evansecent Cutoff =',num2str(uev),...
%' Nyquist Cutoff =', num2str(unp),'u max =',num2str(umax)])
clear ucut uev unp
 
%Fresnel kernel: paraxial approximation
H = exp((-1i*pi*lambda* z).*(A));
clear A;
 
%Truncate kernel
H = W.*H;
clear W;
 
%Hologram Spectrum
Htemp = fft2(Hcrop);
HH = fftshift(Htemp);
clear Htemp;
 
%Propagate field
RR = HH.*H;
clear H HH;
RR =ifftshift(RR);
Ri = ifft2(RR); 
R=abs(Ri.^2);

clear RR;

end
