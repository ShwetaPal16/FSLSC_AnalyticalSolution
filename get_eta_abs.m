function eta_abs = get_eta_abs(a1d,n,P_r2)
% clear all; close all;
% n=1.5;
% P_r2 = .98; % reflection at the bottom surface
% a1d = 1;
% theta = linspace(0,pi/2*.7,181);

N_points = 1801;
theta_esc = asin(1/n);
theta_out = linspace(0,pi/2,N_points); % angle of the incident radiation n=1
theta_in = asin(1/n*sin(theta_out));% angle of the radiation inside the glass n=n
% L_out =2*cos(theta_out); % incoming radiance (function of theta_out)
L_out = 1;
% dthi_dtho = 1/n*cos(theta_out)./sqrt(1-(sin(theta_out)/n).^2);
% L_in= L_out./dthi_dtho*n;

atten = exp(-(a1d)./cos(theta_in));
% A1 = 1-trapz(theta_in,L_in.*atten.*sin(theta_in))/trapz(theta_in,L_in.*sin(theta_in)) % absorption when going down
A2 = get_A2(a1d); %absorption when moving back up;
A1 = 1-trapz(theta_out,L_out.*atten.*sin(theta_out))/trapz(theta_out,L_out.*sin(theta_out)); % absorption when going down

% plot(L_out)
% hold on
% plot(L_in)

% plot(rad2deg(theta_in),L_in)
% trapz(theta_in,sin(theta_in).*L_in)

eta_abs = (A1 + P_r2 * (1 - A1) * A2);
end