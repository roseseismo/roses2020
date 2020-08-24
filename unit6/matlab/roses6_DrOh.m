% Authors: Tolulope Olugboji and Liam Moser
% Key to exercises --- 
%% ExerciseI -- true back azimuth , 280 degrees
disp('EX I Begin Here'); 
%dat = load('NOMELT13_s.txt'); % Z 1 2;
datZ = load('./SampleSeismogram/sample_BHZ_long.txt'); % Z 1 2;
dat1 = load('./SampleSeismogram/sample_BH1_long.txt'); % Z 1 2;
dat2 = load('./SampleSeismogram/sample_BH2_long.txt'); % Z 1 2;

%dat = fliplr(dat);
dat = [datZ, dat1, dat2];  
maxdat = max(abs(dat), [], 'all');
normdat = dat ./ maxdat;            % normalize by largest amplitude

nlen = length(dat);

srate = 50;    % sample rate in Hz
delta = 1/srate;
tvec = 0:delta:(nlen -1)*delta;         % time vector

twin = [2.5,4];
nWin = twin * 60 ./ delta; n1 = nWin(1); n2 = nWin(2);


%
figure(1)
cols = {'k', 'b', 'r'};
clf
hold on

loco = [0.1 0.3];               % corner frequencies to filter waveform
fltdat = zeros(size(normdat));
for ichan = 1:3
    trace = normdat(:,ichan) ;
    
    % bandpass filter from Fred. Simmons
    flted = bpass(trace, srate, loco(1),loco(2), 4, 2, 'butter');
    plot(tvec ./ 60, (ichan-1)*0.1 + flted,  ...
        'linewidth', 2, 'color' , cols{ichan})
    fltdat(:,ichan) = flted;
    
end

legend({'Z', 'H1', 'H2'})
%ylim([-1 5])
%yticks([1 2 3])
%yticklabels({'Z', '1', '2'})
grid on
xlim([0.2 6])
xlabel('Minutes')

% -- find angle by line search for maximum rotation angle in alpha
alpha = 0:5:360;
tval = zeros(size(alpha));  % store transverse energy here

if 1
    
    for i = 1:length(alpha)
        a = alpha(i);
        
        % 2 x 2 rotation matrix in the horizontal plane
        ROTM = [cosd(a), sind(a); -sind(a), cosd(a)];
        
        H12 = [fltdat(:,2), fltdat(:,3)]';
        RT_dat = ROTM * H12;
        tval(i) = RT_dat(2,n1:n2) * RT_dat(2,n1:n2)'; 
        % above, length of transverse channel using inner product
    end
end


 
figure(2)
plot(alpha, tval, 'r-', 'linewidth', 2);
xlabel('Search angle, \theta', 'fontsize', 20)
ylabel(' E(\theta)', 'fontsize', 20)

grid on

[val, ind] = min(tval);
H1rot = alpha(ind); thetat = num2str(H1rot);
H1ang = wrapTo360(280+H1rot);
vv = vline(H1rot); vv.LineWidth = 2;
text(H1rot+2, max(val), ['\theta_t = ' thetat '^\circ'], 'Fontsize', 20)



%% ExerciseII
disp('Ex. II here');
syntrace = load('Rondenay2019_syntrace.txt'); % N E Z - Time;
traces = syntrace(:,1:3);
tvec = syntrace(:,4);

figure(3)
clf
%subplot(131)
%hold on
%for ichn = 1: 3
%    plot(tvec, ichn + traces(:,ichn), 'linewidth', 2)
%end

ylim([-1 5])
yticks([1 2 3])
yticklabels({'N', 'E', 'Z'})
grid on

% rotate into Z R T here 
a = 210;
ROTM1 = [cosd(a), sind(a), 0; -sind(a), cosd(a), 0; 0, 0, 1];

% LQT rotation matrix here ..
alpha = 6.0 ;
rayp = 0.08 ;
a= asind(rayp*alpha); % apparent incidence angle of p
ROTM2 = [cosd(a), sind(a), 0; -sind(a), cosd(a), 0; 0, 0, 1];

% PSV rotation matrix here ..
beta = 3.4; 
q_a  = sqrt(alpha^(-2.0) - rayp^2);% vertical p slowness
q_b  = sqrt(beta^(-2.0) - rayp^2);% vertical s slowness

m11 = (1 - 2*beta^2*rayp^2)/(2*alpha*q_a); m12= (rayp*beta^2)/(alpha);
m21 = (beta*rayp); m22= (2*rayp^2*beta^2  - 1)/(2*beta*q_b);

ROTM3 = [m11, m12, 0; m21, m22, 0; 0, 0, 0.5];

traces_r = ROTM1*traces';
traces_rr = [traces_r(3,:); traces_r(1:2,:)];

subplot(131)
hold on
for ichn = 1: 3
    plot(tvec, ichn + traces_r(ichn,:), 'linewidth', 2)
end

ylim([-1 5])
yticks([1 2 3])
yticklabels({'R', 'T', 'Z'})
grid on


%traces_l = ROTM2* ROTM1*traces';
traces_l = ROTM2* traces_rr;


subplot(132)
hold on
for ichn = 1: 3
    plot(tvec, ichn + traces_l(ichn,:), 'linewidth', 2)
end

ylim([-1 5])
yticks([1 2 3])
yticklabels({'L', 'Q', 'T'})
grid on

subplot(133)
hold on
traces_l = ROTM3* traces_rr;

for ichn = 1: 3
    plot(tvec, ichn + traces_l(ichn,:), 'linewidth', 2)
end

ylim([-1 5])
yticks([1 2 3])
yticklabels({'P', 'SV', 'SH'})
grid on


%% ExerciseIII
clc

disp('Shear Wave Splliting and Anisotropy - III');

dat = importdata('ILON_data.txt', ' ', 4);
tEN = dat.data;
tvec = tEN(:,1);
E = tEN(:,2) - mean(tEN(:,2));
N = tEN(:,3) - mean(tEN(:,3));
nLen = length(E);

dift = diff(tvec);
deltaT = dift(1);   % delta t -- sampling interval in secoNds/sample

figure(4)
clf
subplot(121)
plot(tvec, E, 'k-', 'linewidth', 2)
hold on
plot(tvec, N, 'k--', 'linewidth', 2); text(5, 8, 'North'); text(10, 2, 'East');
grid on

subplot(122)
plot(E,N, 'k-', 'linewidth', 2)
xlabel('East'); ylabel('North');
grid on


% start search through fast axis here ...
phi = 0:1:360; nphi = length(phi);
dt =  deltaT:deltaT:30*deltaT; ndt = length(dt);

[phig, dtg] = meshgrid(phi, dt);

minsurf = zeros(size(phig));  % ratio of eigen values 

%
for i = 1:ndt
    for j = 1:nphi
        
        p = phi(j); 
        ROTM = [cosd(p), sind(p); -sind(p), cosd(p)];
        
        t =  dt(i);
        nshift = floor(t / deltaT);
        
        XY = [E'; N'];
        U_FS = ROTM * XY;   % rotated into fast slow
        ww = [ U_FS(2,nshift:end), U_FS(2,1:nshift-1)] ; % slow is delayed ...
        U_FS(2,:) = ww;
        
     
        %normdat(:,2) = RT_dat(1,:)';
        %normdat(:,3) = RT_dat(2,:)';
        
        C =  (U_FS * U_FS') ./ nLen;
        [V, D] = eig(C);
        [ds, indx] = sort(diag(D),1,'descend');
        
        minsurf(i,j) = ds(1) ./ ds(2);
    end
end

maxval = max(minsurf, [], 'all');
[i, j] =  find(minsurf == maxval, 1)
maxval = minsurf(i, j);
phimax = phig(i, j);
dtmax = dtg(i,j);

figure(5);
clf
[~, c] = polarPcolor(dt,phi,minsurf)
ylabel(c, '\lambda_1 / \lambda_2', 'fontsize', 20)

%pcolor(phig, dtg, minsurf);
%shading flat
%colorbar


%hold on;
%plot(phimax, dtmax, 'wo')


p = phimax; b = 166 - phimax;
ROTM = [cosd(p), sind(p); -sind(p), cosd(p)];
ROTR = [cosd(b), sind(b); -sind(b), cosd(b)];

t =  dtmax;
nshift = floor(t / deltaT);
XY = [E'; N'];

U_FS = ROTM * XY;   % rotated into fast slow and radial



%C =  (XY * XY') ./ nLen;
%[V, D] = eig(C);
%[ds, indx] = sort(diag(D),1,'descend');

figure(6);
clf
subplot(121)
plot(tvec, U_FS(1,:), 'k-', 'linewidth', 2)
hold on
%plot(tvec,U_FS(2,:), 'k--', 'linewidth', 2); text(5, 8, 'North'); text(10, 2, 'East');
plot(tvec,U_FS(2,:), 'r--', 'linewidth', 2); 
%text(5, 8, 'North'); text(10, 2, 'East');
legend({'u_{fast}', 'u_{slow}'})
grid on

subplot(122)
ww = [ U_FS(2,nshift:end), U_FS(2,1:nshift-1)] ; % correct for anisotropy
U_FS(2,:) = ww;   % save correction
XY = ROTR * U_FS;   % rotate into correct orientation

plot(XY(2,:), XY(1,:), 'r-', 'linewidth', 2)   
xlabel('East'); ylabel('North');
grid on

%% ExerciseIV
disp('IV');


