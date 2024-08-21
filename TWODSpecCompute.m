function out = TWODSpecCompute(angle, fc, dist, cSpeed, Mx, My, noiseSub)
% This function used to compute 2D-Spectrum base on given angle
% angle must be provided in [azimuth; elevation] format

ax = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(0:Mx-1)' ...
            *cosd(angle(1, 1))*sind(angle(2, 1)));
ay = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(1:My-1)' ...
            *sind(angle(1, 1))*sind(angle(2, 1)));
a = [ax; ay]; % Steering vector with variable phi and theta
out = (norm(a'*noiseSub).^2);


