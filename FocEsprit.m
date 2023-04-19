%Azimuth and Elevation Direction Finding
%Using Arbitrary Array Geometries
close all;
clear all;

IsSimulatedData = true;

%Create URA
c = 340;       %m/s
f_s = 150000; %sampling rate , Hz (or Sps)
f = 20000;  %Hz

%snap shot number.
%N = 400;
N = 64;

if 1
    %URA
    lambda = c/f;
    d = 0.5 * lambda;
    m = 8;
    mic_pos = -(m-1)/2:(m-1)/2;
    mic_pos = d * mic_pos;

    [X,Y] = meshgrid(mic_pos);

    mic_pos = [X(:), Y(:), zeros(length(X(:)),1)];
else
    createMicDrone;
    mic_pos = [xcfg, ycfg, zeros(length(xcfg),1)];
end

M = size(mic_pos,1);

%theta and phi range.
%thetalim         = [0 180];
thetalim         = [0 90];      %90~180 is just the mirror
%theta_separation = 1;
theta_separation = 5;
%theta_separation = 10;

% Bearing grid
theta_range  = (thetalim(1):theta_separation:thetalim(2))';
Ntheta = length(theta_range);

% range of angle space
philim         = [0 360];
%phi_separation = 10;
phi_separation = 5;
%phi_separation = 1;

% Bearing grid
phi_range  = (philim(1):phi_separation:philim(2))';
Nphi = length(phi_range);

NSpaceGrid = Ntheta * Nphi;

[Thetas, Phis] = meshgrid(theta_range, phi_range);

theta_phi_pairs = zeros(NSpaceGrid, 2);
theta_phi_pairs(:,1) = Thetas(:);
theta_phi_pairs(:,2) = Phis(:);

A = steerVector3D(c, f, mic_pos, theta_phi_pairs);

if IsSimulatedData
    %simulate signals
    %deg(theta, phi), amplitude.
    signals_of_interest = [ ...
        [15,15,1];...
        [20,20,1];...
        ];
    
    %deg(theta, phi), amplitude.
    INR = 50;   %dB compared to signal
    I_amp = 10 ^ (INR/20);
    %I_amp = 1;
    signals_interfering = [ ...
        [70,280,I_amp];...
        
        [60,60,I_amp];...
        ];
    
    thetas = [signals_of_interest(:,1); signals_interfering(:,1)];
    phis = [signals_of_interest(:,2); signals_interfering(:,2)];
    amps = [signals_of_interest(:,3); signals_interfering(:,3)];
    
    signal_theta_phi_pais = [thetas(:), phis(:)];
    A_signals = steerVector3D(c, f, mic_pos, signal_theta_phi_pais);
    
    SNR = 0;
    
    CSM = zeros(M,M);
    for n = 1:N
        %generate signal.
        %random phase
        rph = randn(length(amps),1);
        rph = exp(-j*rph);
        %rph = ones(length(amps),1);
        if 0
            pure_signals = A_signals * (rph .* amps);
            signals  = awgn(pure_signals,SNR);
        else
            signals = zeros(M,1);
            for i=1:length(amps)
                pure_sig_i = A_signals(:,i) * rph(i) * amps(i);
                sig_i = awgn(pure_sig_i,SNR);
                signals = signals + sig_i;
            end
        end
        
        CSM_n = signals * signals' ;
        CSM = CSM + CSM_n;
    end
    
    CSM = CSM / N;
else %IsSimulatedData
    %true signals
        %'/home/wywork/Downloads/DroneData/20221125/22-11-25-15-38-10-BigSoundStatic/Mobile20KStatic_0_29.bin',...
	[signals, CSM]=TrueSignal(...
        '/home/wywork/Downloads/DroneData/20221125/22-11-25-15-44-14-BigSoundCircle/BigSoundCircle_550_579.bin',...
            N);
end

if 0
    %do CSM normalization!
    MaxCsmAbs = max(max(abs(CSM)));
    CSM = CSM/MaxCsmAbs;

    figure;
    surf(abs(CSM));
    xlabel('x');
    ylabel('y');
    title('Element space CSM abs normalized!');
end

%do DAS beamforming to test the signal model.
%CBF(:,iF) = real(diag(conj(squeeze(A(:,:,iF)).') * squeeze(K(:,:,iF))* squeeze(A(:,:,iF))));   

CBF = real(sum(conj(A) .* (CSM* A),1));   

figure;
AbsXY = reshape(abs(CBF),Nphi,Ntheta);
surf(theta_range,phi_range,AbsXY);
hold on;
xlabel('\theta');
ylabel('\phi');

if IsSimulatedData
    text(thetas(1:end),phis(1:end),zeros(1,length(thetas)),'o','color','r');
    hold on;
end

title('Element Space DAS result');

    [V,D] = sorted_eig(CSM);
    if IsSimulatedData
        K = length(thetas);
        En = V(:,end - K + 1:end);
    else
        En = V(:,end-1:end);      %2 signals
    end
    P = En * En';
    CBF = 1./sum(real(conj(A) .* (P * A)),1);

    figure;
    AbsXY = reshape(abs(CBF),Nphi,Ntheta);
    surf(theta_range,phi_range,AbsXY);
    hold on;
    xlabel('\theta');
    ylabel('\phi');

if IsSimulatedData
    text(thetas(1:end),phis(1:end),zeros(1,length(thetas)),'o','color','r');
    hold on;
end

title('MUSIC result in element space');


[V,D] = sorted_eig(CSM);

	tic;
    All_CBFs = IpnsDLBatch(A, CSM, V, D);   %对应所有信号数目的值都在返回结果中了。
    toc;

    figure;
    AbsXY = reshape(abs(All_CBFs(:,M-K)),Nphi,Ntheta);  %对应的信号数目为K
    surf(theta_range,phi_range,AbsXY);
    hold on;
    xlabel('\theta');
    ylabel('\phi');

if IsSimulatedData
    text(thetas(1:end),phis(1:end),zeros(1,length(thetas)),'o','color','r');
    hold on;
end

title('ipnl result in element space');


