clc;
clear all;
close all;
Hsize=2^10;
lambda=0.0633*10^-9;
dpix=4.4*10^-6;
k=(2*pi)/lambda;
d=0.4;
n=0.5;
fm=14.2*10^3;
x= (-Hsize/2:Hsize/2-1)*dpix;
y= (-Hsize/2:Hsize/2-1)*dpix;
umax = 1/(2*dpix); 
df=1/(Hsize*dpix);
% fx=-umax:df:umax-df;
% fy=-umax:df:umax-df;
u=(-Hsize/2:Hsize/2-1)*df;
v=(-Hsize/2:Hsize/2-1)*df;
[X,Y]=meshgrid(x,y);
[U,V]=meshgrid(u,v);
spat_sq = X.^2+Y.^2;
freq_sq= U.^2+V.^2;
Ioz=ones(Hsize);

Imxz=ones(Hsize) + (n.*sin(2.*pi.*fm.*X));
Iqxz=ones(Hsize) + (n.*cos(2.*pi.*fm.*X));
Imyz=ones(Hsize) + (n.*sin(2.*pi.*fm.*Y));
Iqyz=ones(Hsize) + (n.*cos(2.*pi.*fm.*Y));
vari=10^-4;
noise=randn(Hsize).*sqrt(vari);
cut=1;
Hcrop=imread("image2.jpg");
Hcrop=rgb2gray(Hcrop);
Hcrop = padarray(Hcrop,[5 5],0,'both');
Hcrop=imresize(Hcrop,[Hsize Hsize]);
Hcrop=im2double(Hcrop);
Hcrop=Hcrop.*5.6;

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,Hcrop); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Input Phase Object ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"1_Input_phase_object.png");

%%Illumination
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,Imxz); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Structured illumination with x variation ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"1_SIX.png");


figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,Imyz); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Structured illumination with y variation ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"1_SIY.png");
%%Noise
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,noise); 
view(2);   
colorbar
colormap(jet);
shading interp;

%%%Uniform illumination

[uni_rec_ph, deld]=uniform(dpix,d,lambda,Hsize,Hcrop,k,Ioz,freq_sq,noise);

% figure('units','normalized','outerposition',[0 0 1 1]);
% imagesc(x,y,deld); 
% colormap(gray);
% colorbar
% xlabel('x(m)') ;
% ylabel('y(m)');
% axis on
% title(['Intensity differential with unifrom illumination ']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"2_Intensity_differential_with_unifrom_illumination.png");

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,(uni_rec_ph)); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Recovered Phase with uniform illumination ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"2_Recovered_Phase_with_uniform_illumination.png");

%%%Phi_x
[rphi_x,phi_x, deld,delmxd,delqxd]=phix(dpix,d,lambda,Hsize,Hcrop,k,Ioz,Imxz,Iqxz,freq_sq,noise,fm,u,n,cut);
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,real(rphi_x)); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Recovered Phase with structured illumination along x direction ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"3_Recovered_Phase_with_structured_illumination_along_x_direction.png");

%%%Phi_y
[rphi_y,phi_y, deld,delmyd,delqyd]=phiy(dpix,d,lambda,Hsize,Hcrop,k,Ioz,Imyz,Iqyz,freq_sq,noise,fm,u,n,cut);
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,real(rphi_y)); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Recovered Phase with structured illumination along y direction ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"4_Recovered_Phase_with_structured_illumination_along_y_direction.png");

%%%Laplacian
[d_phix_x, d_phix_y] = gradient((phi_x),dpix);
[d_phiy_x, d_phiy_y] = gradient((phi_y),dpix);

LFT_in=fftshift(fft2(d_phix_x+d_phiy_y));
T2= (4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
com_rec_ph= ifft2(ifftshift(T2.*LFT_in));
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,real(com_rec_ph)); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Recovered Phase using laplacian equation with structured illumination along x and y direction']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"5_Recovered_Phase_using_laplacian_equation_with_structured_illumination_along_x_and_y_direction.png");

%%Error analysis

[U,V]=meshgrid(u,v);
T2= (4.*pi.*pi.*(U.^2+V.^2));
T2=T2./(T2.^2 + eps^2);
E_phi_1= (k/d)^2.*(T2.^2).*vari;
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log10(E_phi_1));   
colorbar
colormap(jet);
shading interp;

% 
% [U,V]=meshgrid(u,v);
% T2= ((U.^2));
% T2=T2./(T2.^2 + eps^2);
% E_phi_2= (k/(d*4*pi*pi*fm*n))^2.*(T2).*vari;figure('units','normalized','outerposition',[0 0 1 1]);
% imagesc(u*(d/k)^0.5,u*(d/k)^0.5,log10(E_phi_2));   
% colorbar
% colormap(jet);
% shading interp;
% 

[U,V]=meshgrid(u,v);
T2= ((U.^2+V.^2));
T2=T2./(T2.^2 + eps^2);
E_phi_3= (k/(d*4*pi*pi*fm*n))^2.*(T2).*vari;
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u*(d/k)^0.5,u*(d/k)^0.5,log10(E_phi_3));   
colorbar
colormap(jet);
shading interp;

E_phi_h=(E_phi_1.*E_phi_3)./( E_phi_1 +E_phi_3)
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u*(d/k)^0.5,u*(d/k)^0.5,log10(E_phi_h));   
colorbar
colormap(jet);
shading interp;
%%%Filter
T_ft=( 4.*(U.^2+V.^2))./(4.*(U.^2+V.^2) + (fm*n)^2);
T_ft(isinf(T_ft))=0;
T_ft(isnan(T_ft))=0;


figure
mesh(u,u,abs(T_ft-1));   
colorbar
colormap(jet);
shading interp;

figure('units','normalized','outerposition',[0 0 1 1]);
mesh(u,v,abs(T_ft-1)); 
colormap(jet);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Low pass filter - Frequency response']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6_lp.png");


figure('units','normalized','outerposition',[0 0 1 1]);
mesh(u,v,abs(T_ft)); 
colormap(jet);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['High pass filter - Frequency response']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6_hp.png");


%%Hybrid
[fggg,ix] = min(abs(u-(2*fm)));
index_2fm=ix-(513);
Mcut=zeros(Hsize);
x= (-Hsize/2:Hsize/2-1);
y= (-Hsize/2:Hsize/2-1);
[X,Y]=meshgrid(x,y);
AA=(X.^2+Y.^2)<=index_2fm.^2;
Mcut(AA)=1;

L_ft=fftshift(fft2(real(com_rec_ph))).*imcomplement(T_ft);
Uni_ft=fftshift(fft2(real(uni_rec_ph))).*(1-imcomplement(T_ft));


figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,v,log(abs(L_ft)).*imcomplement(T_ft)); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['\phi_u*t(x,y)']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6_hu.png");


figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,v,log(abs(Uni_ft)).*(T_ft)); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['\phi_s*(1-t(x,y))']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6_ls.png");

final_ft=L_ft+Uni_ft;

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,v,log(abs(final_ft))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['\phi_h']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6_h.png");

hybrid_image=ifft2(ifftshift(final_ft));
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,real(hybrid_image)); 
colormap(gray);
colorbar
xlabel('x(m)') ;
ylabel('y(m)');
axis on
title(['Recovered Phase with Hybridd structured and uniform illumination ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6_Recovered_Phase_with_Hybridd_structured_and_uniform_illumination.png");


function [uni_rec_ph, deld]=uniform(dpix,d,lambda,Hsize,Hcrop,k,Ioz,freq_sq,noise)

Iod=FresProp(dpix,d,lambda,Hsize,exp(1j.*Hcrop));
deld=(Iod-Ioz)+noise;
T1=fftshift((fft2((k/d).*(deld))));
T2= (4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
uni_rec_ph= ifft2(ifftshift(T2.*T1));
end


function [rphi_x,phi_x, deld,delmxd,delqxd]=phix(dpix,d,lambda,Hsize,Hcrop,k,Ioz,Imxz,Iqxz,freq_sq,noise,fm,u,n,cut)


Imxd=FresProp(dpix,d,lambda,Hsize,sqrt(Imxz).*exp(1j.*Hcrop));
Iqxd=FresProp(dpix,d,lambda,Hsize,sqrt(Iqxz).*exp(1j.*Hcrop));
Iod= FresProp(dpix,d,lambda,Hsize,sqrt(Ioz).*exp(1j.*Hcrop));

deld=(Iod-Ioz)+noise;
delmxd=(Imxd-Imxz)+noise;
delqxd=(Iqxd-Iqxz)+noise;

Smxd=(k/d).*((deld.*Imxz)-delmxd);
fft_smx=fftshift(fft2(Smxd));
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log(abs(fft_smx))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Fourier transform of modulation signal ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"14_Fourier_transform_of_modulation_signal.png");


Sqxd=(k/d).*((deld.*Iqxz)-delqxd);

[fggg,ix] = min(abs(u-fm));
index_fm=ix-(513);
[fggg,ix] = min(abs(u-(2*fm)));
index_2fm=ix-(513);
Mcut=zeros(Hsize);
x= (-Hsize/2:Hsize/2-1);
y= (-Hsize/2:Hsize/2-1);
[X,Y]=meshgrid(x,y);
AA=(X.^2+Y.^2)<=index_2fm.^2;
Mcut(AA)=1;
FT_in1=fftshift(fft2(1j.*Smxd));
FT_in1=circshift(FT_in1,index_fm,2);
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log(abs(FT_in1))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Fourier transform of modulation signal ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"14_Fourier_transform_of_modulation_signal_shifted.png");

FT_in2=fftshift(fft2(1j.*Sqxd));
FT_in2=circshift(FT_in2,index_fm,2);
FT_in=FT_in1-FT_in2;
[U,V]=meshgrid(u,u);
HH=(2.*n.*pi.*fm);
H=HH./(HH.^2 + eps^2);
if(cut==1)
    FT_in=(H).*FT_in.*Mcut ;
else
    FT_in=(H).*FT_in;
end
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log(abs(FT_in))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Fourier transform of modulation signal ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"14_Fourier_transform_of_modulation_signal_shifted_cut.png");


phi_x=(ifft2(ifftshift(FT_in)));
HH=(2.*n.*pi.*fm.*1j.*2.*pi.*U);
H=HH./(HH.^2 + eps^2);
if(cut==1)
    FT_in=(H).*FT_in.*Mcut ;
else
    FT_in=(H).*FT_in;
end
rphi_x=imcomplement(ifft2(ifftshift(FT_in)));
end



function [rphi_y,phi_y, deld,delmyd,delqyd]=phiy(dpix,d,lambda,Hsize,Hcrop,k,Ioz,Imyz,Iqyz,freq_sq,noise,fm,u,n,cut)

Imyd=FresProp(dpix,d,lambda,Hsize,sqrt(Imyz).*exp(1j.*Hcrop));
Iqyd=FresProp(dpix,d,lambda,Hsize,sqrt(Iqyz).*exp(1j.*Hcrop));
Iod= FresProp(dpix,d,lambda,Hsize,sqrt(Ioz).*exp(1j.*Hcrop));

deld=(Iod-Ioz)+noise;
delmyd=(Imyd-Imyz)+noise;
delqyd=(Iqyd-Iqyz)+noise;

Smyd=(k/d).*((deld.*Imyz)-delmyd);
fft_smy=fftshift(fft2(Smyd));
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log(abs(fft_smy))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Fourier transform of modulation signal ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"14_Fourier_transform_of_modulation_signal_y.png");

Sqyd=(k/d).*((deld.*Iqyz)-delqyd);

[fggg,ix] = min(abs(u-fm));
index_fm=ix-(513);
[fggg,ix] = min(abs(u-(2*fm)));
index_2fm=ix-(513);
Mcut=zeros(Hsize);
x= (-Hsize/2:Hsize/2-1);
y= (-Hsize/2:Hsize/2-1);
[X,Y]=meshgrid(x,y);
AA=(X.^2+Y.^2)<=index_2fm.^2;
Mcut(AA)=1;

FT_in1=fftshift(fft2(1j.*Smyd));
FT_in1=circshift(FT_in1,index_fm,1);
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log(abs(FT_in1))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Fourier transform of modulation signal ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"14_Fourier_transform_of_modulation_signal_shifted_y.png");

FT_in2=fftshift(fft2(1j.*Sqyd));
FT_in2=circshift(FT_in2,index_fm,1);

FT_in=FT_in1-FT_in2;

[U,V]=meshgrid(u,u);

HH=(2.*n.*pi.*fm);
H=HH./(HH.^2 + eps^2);
if(cut==1)
    FT_in=(H).*FT_in.*Mcut ;
else
    FT_in=(H).*FT_in;
end
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(u,u,log(abs(FT_in))); 
colormap(gray);
colorbar
xlabel('u') ;
ylabel('v');
axis on
title(['Fourier transform of modulation signal ']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"14_Fourier_transform_of_modulation_signal_shifted_cut_y.png");

phi_y=(ifft2(ifftshift(FT_in)));

HH=(2.*n.*pi.*fm.*1j.*2.*pi.*V);
H=HH./(HH.^2 + eps^2);
if(cut==1)
    FT_in=(H).*FT_in.*Mcut ;
else
    FT_in=(H).*FT_in;
end

rphi_y=imcomplement(ifft2(ifftshift(FT_in)));

end




function [R]=FresProp(dpix,d,lambda,Hsize,Hcrop)
 z=d;

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