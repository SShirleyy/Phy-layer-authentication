clc;
clear;

%% Initialization
EbN0_db = 10;                     % Eb/N0 values to simulate (in dB)
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 10;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 1000;                % Size of each data sequence (in nr bits)
nr_training_bits = 100;             % Size of training sequence (in nr bits)
nr_blocks = 10000;                     % The number of blocks to simulate
Q = 1;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.
pulse_shape = ones(1, Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);
sps = 16;

%% Compute f0
 % Generate training sequence.
    b_train = training_sequence(nr_training_bits);
    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);

    % Generate guard sequence.
    b_guard = random_data(nr_guard_bits);
 
    % Multiplex training and data into one sequence.
    b = [b_guard b_train b_data b_guard];
    
    % Map bits into complex-valued QPSK symbols.
    d = qpsk(b);

    % Upsample the signal, apply pulse shaping.
    tx = upfirdn(d, pulse_shape, Q, 1);
%    
%     % Baseband to passband with imbalance
     alpha_tx = 0.03;
     theta_tx = deg2rad(5);
     tx = real(tx)+1i*imag(tx)*(1+alpha_tx)*exp(1i*theta_tx);
    
    rx = tx;
  
    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
   
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start=1+Q*(nr_guard_bits/2); 
    t_end=t_start+5*Q; 
    t_samp = sync(mf, b_train, Q, t_start, t_end);
   
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);
%     r = r/Q;
    %Authentication
    r_training = r(1:nr_training_bits/nr_bits_per_symbol);
    f0 = r_training-qpsk(b_train);


%% Loop over different values of Eb/No.
rp = 1;
H = zeros(1, length(EbN0_db));
h = zeros(1, length(EbN0_db));
for snr_point = 1:length(EbN0_db)
for i=1:rp
  % Loop over several blocks to get sufficient statistics.
  for blk = 1:nr_blocks
    H(snr_point) = 0;
    %%%
    %%% Transmitter
    %%%

    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);

    % Generate guard sequence.
    b_guard = random_data(nr_guard_bits);
 
    % Multiplex training and data into one sequence.
    b = [b_guard b_train b_data b_guard];
    
    % Map bits into complex-valued QPSK symbols.
    d = qpsk(b);

    % Upsample the signal, apply pulse shaping.
    tx = upfirdn(d, pulse_shape, Q, 1);
   
    % Baseband to passband with imbalance
     alpha_tx = -0.05;
     theta_tx = deg2rad(-10);
     tx = real(tx)+1i*imag(tx)*(1+alpha_tx)*exp(1i*theta_tx);
    
    %%%
    %%% AWGN Channel
    %%%
    
    % Compute variance of complex noise according to report.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);

    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(tx))+j*randn(size(tx)));
    % Received signal.
    rx = tx + n;


    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
   
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start=1+Q*(nr_guard_bits/2); 
    t_end=t_start+5*Q; 
    t_samp = sync(mf, b_train, Q, t_start, t_end);
   
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);
%     r = r/Q;
    %Authentication
    r_training = r(1:nr_training_bits/nr_bits_per_symbol);
    f = r_training-qpsk(b_train);
    a = f-f0;
    H(snr_point) = (a*a')/(sigma_sqr/2);
    if H(snr_point)>118.498
    h(snr_point) = h(snr_point) + 1;
    end

    % Phase estimation and correction.
    phihat = phase_estimation(r, b_train);
    r = r * exp(-1i*phihat);
       
    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat = detect(r);
    
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
%     temp=bhat(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
%     nr_errors(snr_point) = nr_errors(snr_point) + sum(temp);

    % Next block.
  end


    % Next Eb/No value.
end
end
% p_fa = h/(rp*50);
p_md = (nr_blocks-h)/(nr_blocks);
disp(p_md);
% plot(EbN0_db, p_fa);
%hold on;plot(EbN0_db, p_md);
% Compute the BER. 
% BER = nr_errors / nr_data_bits / nr_blocks;

 
