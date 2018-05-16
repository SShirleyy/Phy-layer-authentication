clc;
clear;
%% Initialization
EbN0_db = 15;                     % Eb/N0 values to simulate (in dB)
nr_bits_per_symbol = 4;             % Corresponds to k in the report
nr_cyclic_bits = 64;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 224;                % Size of each data sequence (in nr bits)
nr_blocks = 10000;                     % The number of blocks to simulate
nr_pilots = 8;                       % Number of pilot carriers

% sps = 16;



snr_point = 20;
%% Loop over different values of Eb/No.´
 H = zeros(1, length(EbN0_db)*nr_blocks);
 h = zeros(1, length(EbN0_db));
 uu = 0;
 kk= 1;
 for snr_point = 1:length(EbN0_db)
 for blk = 1:nr_blocks
    
    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);
 
    % Multiplex training and data into one sequence.
    b = b_data;

    %Interleaving coded data
    s2=size(b,2);
    q=s2/4;
    matrix=reshape(b,4,q);
    b_i = matrix';
%     b_i = matintrlv(double(matrix'),2,2)'; % Interleave.
    
    % Binary to decimal conversion
    dec=bi2de(b_i,'left-msb');
    dec = dec';
    %16-QAM Modulation
    M=16;
    y = qammod(dec,M, 'UnitAveragePower',true);
    y=y/sqrt(var(y));
    %y = qammod(dec,M);
    pilt=min(y);
    % Pilot insertion
    %lendata=length(y);
    %y1 = zeros(1,64);
    pp = linspace(1,64,nr_pilots);
%     pilt = y(pp);
    y1(pp) = pilt;
    d = pp(2)-pp(1)-1;
    y1(pp(1)+1:pp(2)-1) = y(1:d);
    for u = 2:length(pp)-1
    y1(pp(u)+1:pp(u+1)-1) = y((u-1)*d+1:u*d);
    end
    y1(pp(end)+1:63) = y((length(pp)-1)*d+1:end);
    newvar = var(y1);
    % IFFT
    ifft_sig = ifft(y1,64);
    
    % Add Cyclic Extension
    cext_data = zeros(1,80);
    cext_data(1:16) = ifft_sig(49:64);
    cext_data(17:80) = ifft_sig(1:64);


    %%%
    %%% Multipath Channel
    %%%
    
    ch = [sqrt(3/4), sqrt(1/4)];
    rx = multipath(cext_data, 1, ch);
    % tx = cext_data;

    % Received signal.
   % Compute variance of complex noise according to report.
    %sigma_sqr =  1 / 10^(EbN0_db(snr_point)/10);
    sigma_sqr =  newvar/64 / 10^(EbN0_db(snr_point)/10);
    %rx=rx/sqrt(var(rx));
    
    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(rx))+j*randn(size(rx)));
    % Received signal.
    rx = rx + n;
    
    %%%
    %%% Receiver
    %%%
    
    % Removing Cyclic Extension
    % Which is the first 16 symbols
    rxed_sig = rx(17:end);

    % FFT
    ff_sig=fft(rxed_sig,64);
    
    % Pylot Extraction
    r = ff_sig(pp(1)+1:pp(2)-1);
    for u = 2:length(pp)-1
        r = [r, ff_sig(pp(u)+1:pp(u+1)-1)];
    end
    r = [r, ff_sig(pp(end)+1:63)];
    
    r_pilots = ff_sig(pp);
    %r_pilots(nr_pilots) = ff_sig(64);
    % Phase Estimation
    rt = r_pilots.*conj(pilt);
    phihat = 1/nr_pilots*sum(atan(imag(rt)./real(rt)));
    r = r * exp(-1i*phihat);
    
    % FIR parameter calculation and filitering
    H_t = r_pilots./pilt;
   %H_t = interp1(pp, H_t, 1:64);
    %H_t = spline(pp, H_t, 1:64);
    % Authentication
    theta_A = [sqrt(3/4); sqrt(1/4)]; 
    H_AB = fft(theta_A, 64);
    H_AB = H_AB(pp);
    % Minimizing phase of Alice   
    phi = angle(H_AB'*H_t.');
    H_AB = H_AB.';
    z=H_t-H_AB*exp(1i*phi);
    H(kk) = ((z*z')/(64*newvar*sigma_sqr/(2*(pilt*pilt'))));
    if H(kk)>34.267
    h(snr_point) = h(snr_point) + 1;
    end
    kk = kk+1;
   % Next block.
  end
 % Next Eb/No value.
 end
%  p_FA = h/nr_blocks;
%  disp(mean(p_FA));
p_MD = (nr_blocks-h)/nr_blocks;
disp(p_MD);
%   plot(abs(H_t)); hold on; plot(abs(H_AB));
 % histfit(real(z));
% histfit(imag(z));
%histfit(H,20);
 histogram(H,'Normalization', 'pdf')
 hold on;
 x = 1:50;
 plot(x, chi2pdf(x,nr_pilots*2))