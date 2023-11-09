%% Cleaning
clear all; clc;

%% Spatial coordinates
L = 30*pi; N = 400; h = 2*L/N; x  = (-L+(0:N-1)*h)';
[X,Y] = meshgrid(x,x);

%% Parameters 
alpha_1 = 0.1692;
alpha_2 = 0.1215;
b_1 = 0.76;
b_2 = 0.675;
s_1 = 0.5812;
s_2 = 0.4172;
q = 2*cos(pi/8);

theta = 0.09;
beta = 4.1;

p0(1) = theta;
p0(2) = beta;
p0(3) = L;
p0(4) = alpha_1;
p0(5) = alpha_2;
p0(6) = b_1;
p0(7) = b_2;
p0(8) = s_1;
p0(9) = s_2;
p0(10) = q;

%% Connectivity function
p_fun = @(r,b,s,q) exp(-s*r).*((cos(q*r))+b*sin(q*r));
kernel = @(r, alpha_1, alpha_2, b_1, b_2, s_1, s_2, q) alpha_1 * p_fun(r,b_1,s_1,1) + alpha_2 * p_fun(r,b_2,s_2,q);

%% Kernel and its fft
W = kernel(sqrt(X.^2 + Y.^2), alpha_1, alpha_2, b_1, b_2, s_1, s_2, q);
wHat = fft2(W); 

%% Initial condition
load('u0_N8.mat');

%% Plot synaptic kernel
% figure;
% subplot(1,2,1); plot(x,W(N/2,:),'.-');
% subplot(1,2,2); imagesc(x,x,W); axis square;
% drawnow;

%% function handles
problemHandle = @(t,u) Amarimodel2D(u,p0,wHat);

hFig = figure; set(hFig, 'Position', [230 250 570 510]);

%% solve
tspan = [0 200];
[T,U] = ode45(problemHandle,tspan,u0);

p = p0; u = reshape(U(end,1:N^2),N,N);

%% Plot final state
surf(x,x,u), view(2), 
axis square,  shading interp; colorbar;
set(gca,  'XLim', [-20*pi 20*pi]); set(gca,  'YLim', [-20*pi 20*pi]);

%% save solution
% save('sol2d.mat','sol','p');
