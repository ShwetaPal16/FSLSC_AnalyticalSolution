function A2 = get_A2(a2d)
        theta = linspace(0,pi/2,1001);
        atten = exp(-(a2d)./cos(theta));
        A2 = 1-2*trapz(theta,atten.*cos(theta).*sin(theta));
end