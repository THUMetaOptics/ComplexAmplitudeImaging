%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                                    
%    This code performs a simulation of transverse shearing interferometry using      
%    a metalens.                                                                                                              
                                                                                                                                   
%    Author: Liu Li     Date: 2024.9.3

%    If you find this code useful, please cite the following article:
%    L. Li et al. Single-shot deterministic complex amplitude imaging with a 
%    single-layer metalens. Sci. Adv. 10, eadl0501 (2024).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

%% Parameters of metalens %%
f = 1.5e-2;   % Focal length
D = 2000e-6;    % Diameter
d1 = 2*f;     % Object distance
d2 = d1*f/(d1-f);   % Image distance
lambda = 800e-9;    % wavelength
k=2*pi/lambda;   % wavevector
L0  = 2e-3;  % simulation area
delta_s = 15e-6; % shearing distance
delta_y = ((d2/d1)+1)*delta_s; % real shearing distance
delta_x = 0.2e-3; % image distance between x and y- shearing image

%% Object  %%
N0 = 300;
N1 = round(N0*2/6);

x0 = linspace(-L0/2, L0/2,N0);
y0 = linspace(-L0/2, L0/2,N0);
[X0, Y0] = meshgrid(x0, y0);
scaler = 1.2e-4;
f_scaler = 0.1;
f_gauss = 1/sqrt(2*pi)*exp(- ((2*X0./L0).^2+(2*Y0./L0).^2)./2/f_scaler);
f_gauss = f_gauss./max(max(f_gauss));
Amplitude = f_gauss;

Phase = 0.2*1/2*(3*(1-X0/scaler).^2.*exp(-(X0/scaler).^2-(Y0/scaler+1).^2)-10*(X0/scaler/5-(X0/scaler).^3-(Y0/scaler).^5).*exp(-(X0/scaler).^2-(Y0/scaler).^2)-1/3*exp(-(X0/scaler+1).^2-(Y0/scaler).^2));
Phase = Phase-min(min(Phase));
Uo = Amplitude.*exp(-1i*Phase);

N2 = round(N1*1.3);
figure(1)
subplot(121);imshow(Phase(round(N0/2-N2/2):round(N0/2-N2/2)+N2-1,round(N0/2-N2/2):round(N0/2-N2/2)+N2-1),[]);title('Phase')
subplot(122);imshow(Amplitude(round(N0/2-N2/2):round(N0/2-N2/2)+N2-1,round(N0/2-N2/2):round(N0/2-N2/2)+N2-1),[]);title('Amplitude')

%%  Captured Intensity %%
[I_0,I_90,I_45,I_135] = shearing_interferometry(lambda,d1,L0,d2,Uo,D,delta_s,f,N0,delta_x);
N_cutx = round(N0*6/7);
N_cuty = round(N0*1.3/3);
I_0_temp = I_0(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2-N_cutx/2):round(N0/2-N_cutx/2)+N_cutx-1);
I_45_temp = I_45(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2-N_cutx/2):round(N0/2-N_cutx/2)+N_cutx-1);
I_90_temp = I_90(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2-N_cutx/2):round(N0/2-N_cutx/2)+N_cutx-1);
I_135_temp = I_135(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2-N_cutx/2):round(N0/2-N_cutx/2)+N_cutx-1);
I_max_sum = [I_0_temp, I_45_temp, I_90_temp, I_135_temp];
I_max = max(I_max_sum(:));

figure(2)
subplot(221);imshow(I_0_temp./I_max,[],'border','tight','initialmagnification','fit'); title('I0')
subplot(222);imshow(I_45_temp./I_max,[],'border','tight','initialmagnification','fit'); title('I45')
subplot(223);imshow(I_90_temp./I_max,[],'border','tight','initialmagnification','fit'); title('I90')
subplot(224);imshow(I_135_temp./I_max,[],'border','tight','initialmagnification','fit'); title('I135')

S1 = (I_0 - I_90);S2 = (I_45 - I_135);
phase_gradient = atan2(S2, S1);

%% Reference phase gradient %%
Uo = f_gauss;
[I_0ref,I_90ref,I_45ref,I_135ref] = shearing_interferometry(lambda,d1,L0,d2,Uo,D,delta_s,f,N0,delta_x);
S1 = (I_0ref - I_90ref);S2 = (I_45ref - I_135ref);
phase_gradient_ref = atan2(S2, S1);

%% Calculated phase gradient %%
figure(3)
phase_net = unwrap_phase(phase_gradient-phase_gradient_ref)./(delta_y*1e6);
phase_net(abs(phase_net)>0.02) = nan;
imshow(phase_net(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2-N_cutx/2):round(N0/2-N_cutx/2)+N_cutx-1),[]);
caxis([-max(phase_net(:)) max(phase_net(:))])
title('Phase Gradient')

%% Phase reconstruction %%
N_cut1 = round(N0 * 2*delta_x/L0);
phase_L = phase_net(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2-N_cut1*2):round(N0/2));
phase_R = phase_net(round(N0/2-N_cuty/2):round(N0/2-N_cuty/2)+N_cuty-1,round(N0/2)+1:round(N0/2+N_cut1*2)+1);
dL = L0/N0;
x_length = linspace(0,dL*size(phase_L,2),size(phase_L,2));
y_length = linspace(0,dL*size(phase_L,1),size(phase_L,1));
[X_length,Y_length] = meshgrid(x_length,y_length);
Reconstructed_phase = hfli2(-phase_L,-phase_R,X_length,Y_length)*1e6;
Reconstructed_phase = rot90(Reconstructed_phase-min(Reconstructed_phase(:)),2); 

figure(4)
subplot(121);imshow(Reconstructed_phase,[]);title('Reconstructed Phase')
subplot(122);imshow(Phase,[]);title('Ground Truth')


