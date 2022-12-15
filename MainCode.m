clear all; 
close all;
fontsize = 15;
clc

% constant parameters
QYdye=0.99; %quantum yield of the dye
P_r2=0.98; %Ideality of the lambertian back reflector/walls
h=6.6*10^-34; %Planck's constant
c=3*10^8; %speed of light
fibre=1; %The view angle correction factor for the fibre which is later used to calculate the total incoming irradiance by integrating the spectro-angular reflected intensity due to a Lambertian

%% Calculating the incoming irradiance by integrating the spectro-angular reflected intensity due to a Lambertian
a=load('Lambertian.mat');a=a.a;
b(:,1) = linspace(510,800,151); %decide the step size
b(:,2:9) = interp1( table2array(a(6:end,1)), table2array(a(6:end,2:9)), b(:,1) ); %adjust the step size
 b(:,2:end) =  b(:,2:end).*fibre./P_r2; % correction for the non-ideality of the lambertian and the view angle of the fibre

ThetaExp= table2array( a(5,2:9) ); %convert to matrix

InpSpec(:,1)=b(:,1); InpSpecFlux(:,1)=b(:,1); 
d=  (b(:,2:9).*abs(sind(ThetaExp-85)))  ; %microwatts/s/cm2/nm. The Sine term is due to the solid angle definition in 3D
InpSpec(:,2) = 2*pi.*trapz( (ThetaExp-90)'.*(pi/180), d' ) ; %integration over angle. Careful about the angle definition
InpInt=2*pi.*(trapz( (ThetaExp-90)*(pi/180), trapz( b(:,1),b(:,2:9)) .*abs(sind(ThetaExp-85)) ) ); %Total input irradiance. integration over angle and wavelength.

%same thing but in flux 
b(:,2:end) = (b(:,2:end)).*((b(:,1).*10^-9)./(6.63*10^-34*3*10^8.*10^6)) ; %convert to flux
d=  (b(:,2:9).*abs(sind(ThetaExp-85)))  ; %microwatts/s/cm2/nm
InpSpecFlux(:,2) = 2*pi.*trapz( (ThetaExp-90)'.*(pi/180), d' ) ;
InpIntFlux = 2*pi.*(trapz( ((ThetaExp-90)*(pi/180)), trapz(b(:,1), b(:,2:9).*abs(sind(ThetaExp-85)) ) ) ); 
disp(strcat('Total Input Flux=',num2str(InpIntFlux) ) )
InpSpecFlux(:,2)= InpSpec(:,2) ./(h.*c.*10^6./(InpSpec(:,1).*10^-9) ) ; %micro/s/cm2/nm
clear a b

%% Load Luminophore Data
load('Dye_Measured.mat')
b = a.Final(:,[1 3]);
LumInt = trapz(b(:,1),b(:,2)) ; %Check if dye measurements follow conservation laws. Gives out less photons than input.
clear b

%Now in flux
b = a.FinalFlux(:,[1 3]);
LumIntFlux = trapz(b(:,1), b(:,2) ) ;
disp(strcat('At normal incidence, Flux by the dye/Input=',num2str(LumIntFlux./InpIntFlux) ) ) %the ratio of flux-due-to-dye to total input flux 
%%

% constant parameters
lambdas = linspace(510,800,151); % wavelength from 510nm to 800nm divided into 151 values
% spectrum_fun = @(x) exp(-(x-640).^2/((50/2.3555555)^2)); %GEEN IDEE

LR(:,1)=lambdas';
LR(:,2:size(b,2)) = interp1(b(:,1),b(:,2:end),LR(:,1)) ;
LR(1:20,2:size(b,2))=zeros(20,(size(b,2)-1));

% [spectrum1, l0,e0] = load_emission(lambdas); %converts the downshifted emission of the dye from the original step size to the step size of lamdas
spectrum = LR(:,end)';
spectrum=spectrum./trapz(lambdas, spectrum);
n_substrate = 1.52; %Refractive index of the material

x = linspace(520,800,201); % wavelength from 520nm to 800nm divided into 201 values
[spectrumfull, l1,e1] = load_emission(x); %converts the downshifted emission of the dye from the original step size to the step size of x

figure()
% plot(lambdas, spectrum1, 'LineWidth', 3); hold on
plot(lambdas, spectrum, 'LineWidth', 3); hold on
legend( 'Luminophore Emission'); set(gca, 'FontSize', 20) %in irradiance units

norm_emission = trapz(x,spectrumfull); %Normalising factor
spectrum = spectrum/norm_emission; %Normalise
spectrum = max(spectrum,0); %eliminating the negative values because negative emission makes no sense
[g_func,q_func] = setup_matrix_abstract_v3(); % intensity at every node and total efficiency of the system


%% load reflectance data- Average over s- and p- polarisation
load('RT_OG.mat'); % load the spectro-angular reflectance of the coating

R_in_full = squeeze(result.R_in.Rp)/2 + squeeze(result.R_in.Rs)/2; %Reflectance inside the waveguide that is distort according to Snell's law,
R_out_full = squeeze(result.R_out.Rp)/2 + squeeze(result.R_out.Rs)/2;%%Reflectance in free-sapce 

lambda_full = result.R_in.lambda*1e9; %lambda set by lumerical
thetad = result.R_in.theta; %viewing theta set by lumerical
if thetad(end) >= 90 %eliminate 90degree simulation or edge simulation
    thetad = thetad(1:end-1);
    R_in_full = R_in_full(:,1:end-1);
    R_out_full = R_out_full(:,1:end-1);
end
theta = deg2rad(thetad); %convert to radians
R_in = interp1(lambda_full,R_in_full,lambdas*1.0); %adjust the lambda step size
R_out = interp1(lambda_full,R_out_full,lambdas*1.0); %adjust the lambda step size

T_in = 1 - R_in; %Transmission inside the waveguide. Absorption is assumed to be zero
T_out = 1 - R_out; %Transmission in free-space. Absorption is assumed to be zero
P_r1s = trapz(theta,R_in.*sin(theta'), 2) / trapz(theta,sin(theta)); %ANGULAR NORMALISATION
dthi_dtho = n_substrate^(-2)*cos(theta)./sqrt(1-(sin(theta)/n_substrate).^2); % compensation for refraction
styles = {'-','-','--'}; %this is a styling choice for the angular emission plot

%% Set paramters for calculation and start calculating
figure
for ideal_case = [1,0];
if ideal_case
     P_r2 =1;%.99;
    eta_q = 1;
Q_M_inv = 0;
    A1 = 0; % get absorption before a reflection
    A2 = 0; % get absorption from top to bottom
    eta_abs = 1; % get blue light absorption
else
    P_r2 =.98;%.99;
    eta_q = QYdye;
    T_1 = 1;
    a1d = 1.5;

    Q = 20;
    Q_inv = 1/Q;
    a2d = a1d*Q_inv;
 Q_M_inv = 1/100;
    A1 = get_A1(a2d); % get absorption before a reflection
    A2 = get_A2(a2d); % get absorption from top to bottom
    eta_abs = get_eta_abs(a1d,n_substrate,P_r2); % get blue light absorption
end

%% Calculate q
for i=1:length(lambdas) %wavelength resolved calculation
    q = q_func(A1,A2,P_r1s(i),P_r2,Q_M_inv,eta_q,eta_abs); % actual calculation at nodes
    qs(i,:) = q;
end
L_top = qs(:,4)'; %lluminosity at node 4
L_air = L_top'.*T_out.*dthi_dtho';%luminosity in free space after refraction of the cone

lightout = trapz(theta,L_air.*sin(theta'), 2); %total no of photons out integarted over the angles

%% plot integrated over spectrum

cone_angle = 30; %confinement angle
in_cone = (thetad<=cone_angle)*1.0;  %calculating for photons only within the cone

total_L_air(:,ideal_case+1) = trapz(lambdas,spectrum.*L_air',2); %total photons ideal integrated over the walengths
eta_tot2(ideal_case+1) = trapz(theta,total_L_air(:,ideal_case+1).*sin(theta)); %total photons integrated over the all angles
eta_tot(ideal_case+1) = trapz(lambdas',spectrum'.*qs(:,5)); %total photons integrated over the walengths with all the paarameters of the ocating and dye

eta_em(ideal_case+1)  = trapz(theta,total_L_air(:,ideal_case+1).*sin(theta).*in_cone)/eta_tot(ideal_case+1); %inetgrate only if within the cone and then normalised
eta_to_cone(ideal_case+1)  = trapz(theta,total_L_air(:,ideal_case+1).*sin(theta).*in_cone); %inetgrate only if within the cone
eta_to_cone_lam  = trapz(theta,2*cos(theta).*sin(theta).*in_cone); %efficiecny for a lambertian

plot(thetad, total_L_air(:,ideal_case+1),styles{ideal_case+1},'linewidth',2.5)
hold on
% plot(thetad, total_L_air(:,ideal_case+1)./cos(theta)/2,styles{ideal_case+1},'linewidth',2)

end

plot(thetad,1.96*cos(theta),styles{3},'linewidth',2.5) %plot lambertian

% plot(thetad,theta./theta)
% plot(thetad,theta./theta*1.5)

legend('Realistic','Ideal', '98% Lambertian')
xlabel('Emission Angle [deg]')
ylabel('Emission [sr^{-1}]')
set(gca,'fontsize',fontsize)
grid on
xlim([0,90])
% ylim([0,4])
%% Final

Final =spectrum.*L_air';
Final=Final.*InpIntFlux;

a=trapz(lambdas, Final(:,:)' ); %integrate the FSLSC output over wavelength
FSLSCphotons= trapz(theta,a'.*sin(theta)) ; %integrate the FSLSC output over solid angle
trapz(theta,InpIntFlux.*sin(theta).*total_L_air); %Flux by the FSLSC real and ideal FSLSC
disp('ideal and real FSLSC Output/Input=') 
trapz(theta,InpIntFlux.*sin(theta).*total_L_air)./InpIntFlux

Lambphotons= InpIntFlux.*trapz(theta,(1.96*cos(theta)).*sin(theta)) ; %Analytical caluclation of flux by 98% ideal lambertian
clear a
%% Plot Spectro-angular Reflectance inside the waveguide

figure()
fontsize=15;
imagesc(lambdas,rad2deg(theta),flipud(min(R_in',10)))
yticklabels([80:-10:0])
xlabel('Wavelength (nm)')
ylabel('Emission Angle [deg]')
improve_figure;
ylabel(h,'Reflectance')
set(gca,'LineWidth',1,'TickLength',[0.015 0.015]);

%% Plot Spectro-angular Reflectance free-space
figure()
fontsize=15;
imagesc(lambdas,rad2deg(theta),flipud(min(R_out',10)))
yticklabels([80:-10:0])
xlabel('Wavelength (nm)')
ylabel('Emission Angle [deg]')
improve_figure;
ylabel(h,'Reflectance')
set(gca,'LineWidth',1,'TickLength',[0.015 0.015]);
%% Plot spectro-angular emission
figure()
fontsize=15;
imagesc(lambdas,rad2deg(theta),flipud(Final./(2*pi)))
yticklabels([80:-10:0])
xlabel('Wavelength (nm)')
ylabel('Emission Angle [deg]')
improve_figure;
ylabel(h,'Emission [sr^{-1}nm^{-1}]')
set(gca,'LineWidth',1,'TickLength',[0.015 0.015]);
