clc;
clear;
delete(findall(0, 'Type', 'figure'));

% STEP a: Simulating the Narrowband Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 100; % Number of time snapshots
fs = 10^7; % Sampling frequency
fc = 10^6; % Center frequency of narrowband sources
Mx = 12; % Number of array elements on x axis
My = 8; % Number of array elements on y axis
N = 4; % Number of sources
sVar = 1; % Variance of the amplitude of the sources

% p snapshots of N narrowband sources with random amplitude of mean zero
% and covariance 1: (N x P)
s = sqrt(sVar)*randn(N, p).*exp(1i*(2*pi*fc*repmat((1:p)/fs, N, 1)));
% STEP a: Simulating the Narrowband Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% STEP b: Mixing the sources and getting the sensor signals %%%%%%%%%%%%%%
doa = [34,66;
       12,45;
       78,20;
       65,32]; %DOAs, azimuth and elevation in first and second column

cSpeed = 3*10^8 ; % Speed of light
dist = 180; % Sensors (i.e., antennas) spacing in meters in both x and y axes

% Constructing the Steering matrix
Ax = zeros(Mx, N);
Ay = zeros(My - 1, N);
for k = 1:N
    Ax(:, k) = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(0:Mx-1)' ...
                    *cosd(doa(k, 1))*sind(doa(k, 2))); 
    Ay(:, k) = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(1:My-1)' ...
                    *sind(doa(k, 1))*sind(doa(k, 2)));
end
A = [Ax; Ay]; % Steering matrix of all sensors: ((Mx + My - 1) x N))

noiseCoeff = 1; % Variance of added noise
x = A*s + sqrt(noiseCoeff)*randn(Mx + My - 1, p); 
% Sensor signals: ((Mx + My - 1) x P))
% STEP b: Mixing the sources and getting the sensor signals %%%%%%%%%%%%%%


% STEP c: Estimating the covariance matrix of the sensor array %%%%%%%%%%%
R = (x*x')/p; % Empirical covariance of the antenna data
% STEP c: Estimating the covariance matrix of the sensor array %%%%%%%%%%%


% STEP d: Finding the noise subspace and estimating the DOAs %%%%%%%%%%%%%
[V, D] = eig(R);
noiseSub = V(:, 1:Mx + My - 1 - N); % Noise subspace of R
% dimensionality: ((Mx + My - 1) x (Mx + My - 1 -N))

phi = 0:0.1:90; % searching scope of azimuth
theta = 0:0.1:90; % searching scope of elevation
% assuming incident sources are restricted in first octant

res = zeros(length(phi),length(theta));
% spectrum value with respect to phi and theta being searched
for n = 1:length(phi)
    for m = 1:length(theta)
        ax = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(0:Mx-1)' ...
                    *cosd(phi(n))*sind(theta(m)));
        ay = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(1:My-1)' ...
                    *sind(phi(n))*sind(theta(m)));
        a = [ax; ay]; % Steering vector with variable phi and theta
        res(n, m) = (norm(a'*noiseSub).^2);
    end
end

% [resSorted, orgInd] = sort(res, 'descend');
% DOAs = theta(orgInd(1:N, 1));
% STEP d: Finding the noise subspace and estimating the DOAs %%%%%%%%%%%%%

% STEP e: Plotting spectrum vs phi & theta (deg)
figure;
mesh(theta, phi, res);
xlabel('Elevation(°)');
ylabel('Azimuth(°)');
zlabel('Norm^2');
title('spectrum (Norm^2) vs phi & theta (deg)');
% STEP e: Plotting spectrum vs theta(deg)

% STEP f: Use ULA on X-axis/Y-axis to calculate the angle of sources 
% with X-axis/Y-axis, denoted by Xdoa/Ydoa
RMx = R(1:Mx, 1:Mx); % covariance matrix of signals received by ULA on X-axis
[VMx, DMx] = eig(RMx);
noiseSubMx = VMx(:, 1:Mx - N); % Noise subspace of RMx: (Mx x (Mx - N))
% Note: relationship (N < Mx) is required
Xdoa = 0:0.1:90;
resMx = zeros(length(Xdoa), 1);
for n = 1:length(Xdoa)
    ax = exp(-1i*2*pi*fc*dist*(1/cSpeed)*cosd(Xdoa(n))*(0:Mx-1)');
    % This expression is used for 1D-Searching
    resMx(n, 1) = (norm(ax'*noiseSubMx).^2);
end

RMy = R(Mx+1:Mx+My-1, Mx+1:Mx+My-1); % Original sensor is omitted
RMy = [R(Mx+1:Mx+My-1, 1), RMy];
RMy = [[R(1, 1), R(1, Mx+1:Mx+My-1)]; RMy]; % (My x My)
% covariance matrix of signals received by ULA on Y-axis
[VMy, DMy] = eig(RMy);
noiseSubMy = VMy(:, 1:My - N); % Noise subspace of RMx: (My x (My - N))
% Note: relationship (N < My) is required
Ydoa = 0:0.1:90;
resMy = zeros(length(Ydoa), 1);
for n = 1:length(Ydoa)
    ay = exp(-1i*2*pi*fc*dist*(1/cSpeed)*cosd(Ydoa(n))*(0:My-1)');
    % This expression is used for 1D-Searching
    resMy(n, 1) = (norm(ay'*noiseSubMy).^2);
end
% Plot spectrum vs Xdoa/Ydoa
figure;
plot(Xdoa, resMx, 'b-', Ydoa, resMy, 'r-');
hold on;
alpha = acosd(cosd(doa(:, 1)) .* sind(doa(:, 2)));
beta = acosd(sind(doa(:, 1)) .* sind(doa(:, 2)));
% alpha and beta are angles of sources with X-axis and Y-axis
for i = 1:size(doa, 1)
    xline(alpha(i), 'b--', 'Label', ['alpha', ...
                    num2str(i), '=', num2str(alpha(i))]);
    xline(beta(i),  'r--', 'Label', ['beta', ...
                    num2str(i), '=', num2str(beta(i))]);
end
% Annotate angles of sources with X-axis and Y-axis in figure
hold off;
% STEP f

% STEP g: Estimate angles of sources with X-axis and Y-axis according to
% curves "spectrum vs Xdoa/Ydoa"
estN = numSourcesa(D); % Estimate number of sources based on D matrix
estalpha = ONEDSrch([Xdoa', resMx], estN);
estbeta = ONEDSrch([Ydoa', resMy], estN);
% Estimate angles of sources with X-axis and Y-axis, denoted as 
% estalpha and estbeta. Number of Elements contained in 
% estalpha and estbeta is equal or less than 'estN'.
estdirec = [];
for i=1:size(estalpha, 1)
    for j=1:size(estbeta, 1)
        if cosd(estalpha(i, 1))^2 + cosd(estbeta(j, 1))^2 <= 1
            estphi = atand(cosd(estbeta(j, 1))/cosd(estalpha(i, 1)));
            esttheta = acosd(sqrt(1 - ...
            cosd(estalpha(i, 1))^2 - cosd(estbeta(j, 1))^2));
            estdirec = [estdirec; [estphi, esttheta]];
        end
    end
end
% Estimate azimuth and elevation. There's at most 'estN^2' number of 
% elements contained in estdirec. They will be further filtered to 'estN'
% elements.
temp = zeros(size(estdirec, 1), 1);
for i = 1:length(estdirec)
    ax = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(0:Mx-1)' ...
                *cosd(estdirec(i, 1))*sind(estdirec(i, 2)));
    ay = exp(-1i*2*pi*fc*dist*(1/cSpeed)*(1:My-1)' ...
                *sind(estdirec(i, 1))*sind(estdirec(i, 2)));
    a = [ax; ay]; % Steering vector with variable phi and theta
    temp(i, 1) = (norm(a'*noiseSub).^2);
end
estdirec = [estdirec, temp];
[~, idx] = sort(estdirec(:, 3));
estdirec = estdirec(idx, :);
LocSrchP = estdirec(1:estN, :);
% Sort 'estdirec' according to spectral value, then retain 'estN' 
% elements for following Local Fine 2D-Searching.
% STEP g

% STEP h: Local Fine 2D-Searching. Use filtered 'LocSrchP' as starting 
% point, perform four-directions searching to find minimum points
mov = [[0, 1]; [-1, 0]; [0, -1]; [1, 0]]; % right, up, left, down
SrchedPnt = zeros(size(LocSrchP, 1), 2);
for i = 1:size(LocSrchP, 1)
    SrchPnt = LocSrchP(i, 1:end - 1);
    cnt = 1;
    loop = 0;
    temp = 0;
    endflag = cnt;
    crtPnt = SrchPnt;
    crtspec = TWODSpecCompute(crtPnt', fc, dist, cSpeed, ...
    Mx, My, noiseSub);
    movdirec = mov(cnt, :);
    nxtPnt = crtPnt + movdirec * 0.1;
    nxtspec = TWODSpecCompute(nxtPnt', fc, dist, cSpeed, ...
    Mx, My, noiseSub);
    for j = 1:100
    % This loop will perform four-directions searching for a specific 
    % starting point. Assuming the starting point is close enough to 
    % minimum point, entire searching contains no more than 100 times of 
    % 2D spectrum computation. Right will be the initial direction, then 
    % searching direction will change in a counterclockwise manner.
    % If spectrum of current point is larger than that of next point along 
    % a certain direction, then current point will move along this 
    % direction, and this same direction will be next first searching 
    % direction. If spectrum of current point is less than that of points 
    % along all four-directions, this loop will 'break' and current point
    % will be treated as minimum point, then this loop will be performed
    % on other starting points (stored in 'LocSrchP'). Note that some
    % compares can be omitted during searching process, that's why variable
    % 'temp' is set.
        if crtspec <= nxtspec
            if cnt == 4, cnt = 1; loop = 1; else, cnt = cnt + 1; end
            if cnt == temp, cnt = cnt + 1; end
            if cnt == endflag && loop == 1, break; end
            movdirec = mov(cnt, :);
            nxtPnt = crtPnt + movdirec * 0.1;
            nxtspec = TWODSpecCompute(nxtPnt', fc, dist, cSpeed, ...
                                    Mx, My, noiseSub);
        else
            crtPnt = crtPnt + movdirec * 0.1;
            nxtPnt = crtPnt + movdirec * 0.1;
            endflag = cnt;
            loop = 0;
            temp = mod(cnt + 2, 4);
            crtspec = TWODSpecCompute(crtPnt', fc, dist, cSpeed, ...
                Mx, My, noiseSub);
            nxtspec = TWODSpecCompute(nxtPnt', fc, dist, cSpeed, ...
                                    Mx, My, noiseSub);
        end
    end
    SrchedPnt(i, :) = crtPnt;
end
% Step h


