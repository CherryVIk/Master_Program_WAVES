clear all;
close all;
%% SONAR Parameters
% Speed of sound underwater
cWater = 1500.0;
% Number of transmitting projectors
NTx = 1; % SIMO: Single-input
centerTx = [0 0 0];
% Number of receiving hydrophones
NRx = 32; % SIMO: Multi-output
centerRx = [0 0 0];
% Beamforming angles
angles = -90:2:90;
% Number of beams
NBeams = length(angles);
time_end = 20;
bubbleVelocity = 0.8; % v = 1m / 0.1s;
for time_step = 1:time_end % can be assumed as 0.1s
%% Signal Parameters
% Sample frequency
fs = 192000;
% Signal bandwidth
fB = 20000;
% Center frequency
fC = 70000;
% Min & max frequencies
fMin = fC - fB/2;
fMax = fC + fB/2;
% Ping duration
tPing = 1.0; % in seconds
nPing = tPing * fs; % in samples
% Signal duration
tSig = 0.05; % in seconds
nSig = tSig * fs; % in samples
t = linspace(0, tSig, nSig);
% Signal types
eSignalTypes.CW = 'CW';
eSignalTypes.blNoise = 'blNoise';
eSignalTypes.HFM = 'HFM';
eSignalType = eSignalTypes.blNoise;
bRandomSeed = 42;%42
rng(bRandomSeed)

%% Adaption of SONAR-system wrt. signal design
% Set array element distance to half of wavelength of max freq.
dRx = cWater / (fMax) / 2; % Lambda half distance of array elements
dTx = dRx;

%% FFT-parameters
nNextPow2 = nextpow2(nSig*2); % find nearest x so 2^x = nSig
NFFT = 2^nNextPow2; % FFT-length as multiple of 2
% NFFT = nSig;
NBins = NFFT / 2 + 1; % FFT-bins of pos. frequencies
bins = 0:NBins-1; % Freq. bin support vector
fBinRes= fs / NFFT;
nfMin = floor(fMin/fBinRes);
nfMax = ceil(fMax/fBinRes);
f = linspace(-fs/2, fs/2, NFFT);%NFFT
% f = linspace(0, fs/2, NBins);%NFFT

%% Generate transmit sequence
% Bandpass filter design
Fstop1 = fMin-100;       % First Stopband Frequency
Fpass1 = fMin;       % First Passband Frequency
Fpass2 = fMax;       % Second Passband Frequency
Fstop2 = fMax+100;       % Second Stopband Frequency
Astop1 = 100;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 100;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, fs);
Hd = design(h, 'cheby2');



if strcmp(eSignalType, eSignalTypes.blNoise)
    % Generate Gaussian white noise
    tx = randn(nSig, NTx);
    tx = filter(Hd, tx);
    %tx = filtfilt(Hd.sosMatrix, Hd.ScaleValues, tx);
    % Transform time to freq. domain signal
    Tx = fft(tx, NFFT);%NFFT
    % Only save positive freq.
    Tx = Tx(1:NBins, :);
% The following has a commented out ideal (but impractical) bandpass filter
    % Bandpass
    %Tx(1:nfMin, :) = 0;
    %Tx(nfMax:end, :) = 0;
    % Freq. to time domain
    %tx = ifft(Tx, NFFT, 'symmetric');
    %tx = tx(1:nSig);
% End of ideal bandpass example
end

%% Plot transmit sequence

% fig=figure(5);
% subplot(211);
% grid on;
% best_plot_ever(t, tx(:,1),"Time domain signal", fig)
% subplot(212);
% logTx = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
% grid on;
% best_plot_ever(f(end/2:end), logTx,"Log. frequency spectrum ", fig)

%% Bubble environment settings
%  source location constrains a, b
x_lims=[50 51];
y_lims=[50 51];
z_lims=[0 0];
Nbubbles=10;
bubbleOsc_lims = [-1,1];
maxRadius = 1000e-6;
minAllowableDistance = max([585e-6, 2 * maxRadius]);
if time_step == 1
    % Generate bubbles in some constrained space
    posTar = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles);
else
    rng(time_step)
    % posTarNew = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles);
    posTar(:,3) = posTar(:,3) + bubbleVelocity*ones(NTargets, 1);
    bubbleOscillations = bubbleOsc_lims(1) + (bubbleOsc_lims(2) - bubbleOsc_lims(1))*rand(Nbubbles,2);
    posTar(:,1:2) = posTar(:,1:2) + bubbleOscillations;
    % posTar = [posTar; posTarNew];
end
NTargets = size(posTar, 1);
bDirectSound = 0;

x = posTar(:,1);
y = posTar(:,2);
z = posTar(:,3);
bubbles_mov = figure(10);
hold on
plot3(x,y,z, '-ok',"MarkerEdgeColor" ,	"#4DBEEE")
axis_x = x_lims+bubbleOsc_lims*time_end*bubbleVelocity;
axis([axis_x   axis_x     0 time_end*bubbleVelocity])
grid on
% axis square
view([10  20])
refreshdata
drawnow
Frame = getframe(bubbles_mov);
make_gif(Frame, time_step, "Bubble_mov3.gif");
% hold off
end

%% Functions
function best_plot_ever(x,y,titlename, fig)
    plot(x,y);
    title(titlename);
    font = 20;
    a = 36; % set this parameter and keep it forever
    b = 0.55; % feel free to play with this ratio
    set(gca,'FontSize',font)
    set(findall(fig,'-property','Box'),'Box','off') % optional
    set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
    set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(fig,'Units','centimeters','Position',[3 3 a b*a])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    set(gca, 'XDir', 'normal', 'YDir', 'normal');
    grid on
end
function posTar = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles)
    posTarX = x_lims(1) + (x_lims(2)-x_lims(1))*rand(Nbubbles,1);
    posTarY = y_lims(1) + (y_lims(2)-y_lims(1))*rand(Nbubbles,1);
    posTarZ = z_lims(1) + (z_lims(2)-z_lims(1))*rand(Nbubbles,1);
    posTar = [posTarX posTarY posTarZ];
end
function make_gif(Frame, ii, filename)
    im = frame2im(Frame); 
    [imind, CM] = rgb2ind(im,256); 
    % Write the animation to the gif File: MYGIF 
    if ii == 1 
      imwrite(imind, CM,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind, CM,filename,'gif','WriteMode','append'); 
    end 
end