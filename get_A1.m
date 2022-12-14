function A1 = get_A1(a2d)
% clear all
% a2d=.00001;
% 
n=0;
NN = 10;
for N = 2.^[NN,NN+1]
    n=n+1;
%     a2x = linspace(0,a2d,N)';
%     theta = linspace(0,pi/2,N);
x = linspace(1/2/N,1,N);
%     atten = exp(-(a2d-a2x)./abs(cos(theta)));
    
%     Ax = 1-trapz(theta,atten.*sin(theta),2)/trapz(theta,sin(theta),2);
%     A1 = trapz(a2x,Ax)/a2d
%     A_temp = sin(theta).*cos(theta).*exp(-a2d./cos(theta));
%     A1 = 1-1/a2d*(1/2-trapz(theta,A_temp))
    A2(n) = 1+1/a2d*(trapz(x,x.*(exp(-a2d./x)-1)));
end
% integral = igamma(-2,a2d)*a2d;
% A1 = 1-1/2/a2d+integral;% A2 = extrapolate(fliplr(A1));


p=2;
A1=2^p/(2^p-1)*A2(2)-1/(2^p-1)*A2(1);
dA = A1-A2(2);
if dA/A1>.01
    disp('A1 is not accurate')
end
end