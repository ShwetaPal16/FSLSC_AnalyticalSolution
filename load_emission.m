function [em, lambdar, emn] = load_emission(lambda0)
load('emn.mat')
em = interp1(lambdar,emn,lambda0); %adjusting teh step size
em = em/trapz(lambda0,em); % WHY NORMALISE? FYI: the integral is 1
end

