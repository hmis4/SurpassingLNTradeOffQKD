close all

%% Mg0:LiNbO3

% Sellmeier
% From: Gayer, O. et al. Appl. Phys. B 91, 343–348 (2008).
a1_e = 5.756; a1_o = 5.653;
a2_e = 0.0983; a2_o =  0.1185;
a3_e =  0.2020; a3_o =  0.2091;
a4_e =  189.32; a4_o =  89.61;
a5_e =  12.52; a5_o =  10.85;
a6_e =  1.32e-2; a6_o =  1.97e-2;
b1_e =  2.860e-6; b1_o =  7.941e-7;
b2_e =  4.700e-8; b2_o =  3.134e-8;
b3_e =  6.113e-8; b3_o =  -4.641e-9;
b4_e =  1.516e-4; b4_o =  -2.188e-6;

f = @(T) (T-24.5)*(T+570.82);

% refractive indices 
ne = @(x,T) sqrt(a1_e + b1_e*f(T) + (a2_e + b2_e*f(T))./(x.^2 - (a3_e - b3_e*f(T)).^2) + (a4_e + b4_e*f(T))./(x.^2 - a5_e^2) - a6_e*x.^2);
no = @(x,T) sqrt(a1_o + b1_o*f(T) + (a2_o + b2_o*f(T))./(x.^2 - (a3_o - b3_o*f(T)).^2) + (a4_o + b4_o*f(T))./(x.^2 - a5_o^2) - a6_o*x.^2);


l = linspace(1.53,1.56,2500);
[X,Y] = meshgrid(l,l);

sinc = @(x) sin(x)./x;

c = 299792458;

P = 1./(1./X + 1./Y); % pump wavelength from energy conservation
Pw = (c./X + c./Y)*1e6; % pump frequency

pump = 'e'; % e
signal = 'e'; % o 
idler = 'e'; % e

lp = 0.7722*1e-6; % wavelength
fp = c/lp;
dtau = 1/6*3.5e-12/1.76; % pump duration;
dtau2 = 1*3.5e-12/1.76; % 
df = 1/(4*dtau); % frequency width
df2 = 1/(4*dtau2); 

alpha_p = sech((fp-Pw)/df); % pump envelope
alpha_p2 = sech((fp-Pw)/df2);

% crystal temp
T = 82;


%% refractive indices %%
switch signal
    case 'e'
        ns = ne(X,T);
    case 'o'
        ns = no(X,T);
end

switch idler
    case 'e'
        ni = ne(Y,T);
    case 'o'
        ni = no(Y,T);
end

switch pump
    case 'e'
        np = ne(P,T);
    case 'o'
        np = no(P,T);
end

%%

texp = 13.3e-6;
texp_pp = 9.7e-9;

L = 3e3 * (1+(T-25)*texp); % crystal length
L2 = 40e3 * (1+(T-25)*texp); % reference crystal length
pp = 1*19.1 * (1+(T-25)*texp_pp); % poling period

dk = 2*pi*(np./P - ns./X - ni./Y - 1/pp);

% phase matching func
pm = sin(L/2 * dk)./(L/2 * dk);
pm2 = sin(L2/2 * dk)./(L2/2 * dk);

tks = linspace(min(l),max(l),5);

fs = 22;
tfs = 14;
%% plotting
% ax1 = figure;
% surf(X,Y,abs(pm).^1);shading flat;colormap jet;axis equal;
% % ylim([min(l) max(l)]);xlim([min(l) max(l)]);
% view([0 90]);
% ax = gca;ax.FontSize = tfs;
% xlabel('$\lambda_s \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);ylabel('$\lambda_i  \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);
% xticks(tks);yticks(tks);
% ax1.Position = [441 379 450 430];set(gca,'position',[0.18 0.087 0.78 0.88])

%%

% ax2 = figure;
% surf(X,Y,abs(alpha_p).^1);shading flat;colormap jet;axis equal;
% ylim([min(l) max(l)]);xlim([min(l) max(l)]);view([0 90]);xticks(linspace(1.54,1.56,5));yticks(linspace(1.54,1.56,5));
% ax = gca;ax.FontSize = tfs;
% xlabel('$\lambda_s \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);ylabel('$\lambda_i  \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);
% xticks(tks);yticks(tks);
% ax2.Position = [441 379 450 430];set(gca,'position',[0.16 0.085 0.80 0.9])
% 
% 
% %%
% ax3 = figure;
% surf(X,Y,(abs(pm)+abs(alpha_p)).^1);shading flat;colormap jet;axis equal;
% ylim([min(l) max(l)]);xlim([min(l) max(l)]);view([0 90]);xticks(linspace(1.54,1.56,5));yticks(linspace(1.54,1.56,5));
% ax = gca;ax.FontSize = tfs;
% xlabel('$\lambda_s \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);ylabel('$\lambda_i  \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);
% xticks(tks);yticks(tks);
% ax3.Position = [441 379 450 430];set(gca,'position',[0.16 0.085 0.80 0.9])
% %%
% ax4 = figure;
% surf(X,Y,abs(pm.*alpha_p).^1);shading flat;colormap jet;axis equal;
% ylim([min(l) max(l)]);xlim([min(l) max(l)]);view([0 90]);xticks(linspace(1.54,1.56,5));yticks(linspace(1.54,1.56,5));
% ax = gca;ax.FontSize = tfs;
% xlabel('$\lambda_s \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);ylabel('$\lambda_i  \hspace{1mm}(\mu m)$','interpreter','latex','fontsize',fs);
% xticks(tks);yticks(tks);
% ax4.Position = [441 379 450 430];set(gca,'position',[0.16 0.085 0.80 0.9])
% 
% %%
% figure;
% JSI = (pm.*alpha_p).^2;
% plot(l,sum(JSI,1));hold on;plot(l,sum(JSI,2));
% % colormap jet;

%%
sgm = 1e-3;

alpha_pX = exp(-(1/2 * (X-2*lp*1e6+5e-3).^2/(sgm)^2 ));
alpha_pY = exp(-(1/2 * (Y-2*lp*1e6-5e-3).^2/(sgm)^2 ));
alpha_pX2 = exp(-(1/2 * (X-2*lp*1e6-5e-3).^2/(sgm)^2 ));
alpha_pY2 = exp(-(1/2 * (Y-2*lp*1e6+5e-3).^2/(sgm)^2 ));
ax4 = figure;
surf(X,Y,abs(pm.*alpha_p.*alpha_pX.*alpha_pY).^1+(pm.*alpha_p.*alpha_pX2.*alpha_pY2).^1);shading flat;colormap jet;axis equal;
ylim([min(l) max(l)]);xlim([min(l) max(l)]);view([0 90]);xticks(linspace(min(l),max(l),5));yticks(linspace(min(l),max(l),5));

JSI_flt = abs(pm.*alpha_p.*alpha_pX.*alpha_pY).^2+(pm.*alpha_p.*alpha_pX2.*alpha_pY2).^2;
%%
figure
surf(X,Y,abs(alpha_p.*pm).^1+abs(alpha_pY)+abs(alpha_pX));shading flat;colormap jet;axis equal;
% ylim([min(l) max(l)]);xlim([min(l) max(l)]);
view([0 90]);xticks(linspace(min(l),max(l),5));yticks(linspace(min(l),max(l),5));

%% top-hat bandpass filters
w = 1e-3;
idx1 = l < 2*lp*1e6-w/2;
idx2 = l > 2*lp*1e6+w/2;
[alpha_pX,alpha_pY] = deal(ones(numel(l),numel(l)));
alpha_pX(idx1,:) = 0;alpha_pX(idx2,:) = 0;
alpha_pY(:,idx1) = 0;alpha_pY(:,idx2) = 0;
figure
surf(X,Y,abs(pm2.*alpha_p.*alpha_pX.*alpha_pY).^1);shading flat;colormap jet;axis equal;
ylim([min(l) max(l)]);xlim([min(l) max(l)]);view([0 90]);xticks(linspace(min(l),max(l),5));yticks(linspace(min(l),max(l),5));
JSI = abs(pm2.*alpha_p.*alpha_pX.*alpha_pY).^2;

%
disp(['Relative generation rate for 1nm bandwidth (accounting for crystal length) ' num2str((trapz(JSI_flt(:))/40*3)/trapz(JSI(:)))])