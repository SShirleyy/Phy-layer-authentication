clc;
clear;
addpath('simulations');
%% Initialization
 EbN0_db = linspace(0,30,10);         % Eb/N0 values to simulate (in dB)
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 10;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 2000;                % Size of each data sequence (in nr bits)
nr_training =[128 192 256 512];     % Size of training sequence (in nr bits)
nr_blocks = 500;                     % The number of blocks to simulate
Q = 1;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.
pulse_shape = ones(1, Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);

%% Compute f0
 % for t_point = 1:4;
t_point = 4;
nr_training_bits = nr_training(t_point);
% Generate training sequence.
u_train = training_sequence(nr_training_bits/2);

M = zeros(1, length(EbN0_db));
T = zeros(1, length(EbN0_db));

 P = linspace(0.01,0.1,10);
 h = zeros(1, length(P));
 p_md = zeros(1, length(P));
   for i = 1:length(P);
       p = P(i);
       m = zeros(1, length(EbN0_db));
for snr_point = 1:length(EbN0_db)
        for blk = 1:nr_blocks
    %%%
    %%% Transmitter
    %%%
    
    %repeat training sequence
    b_train = repmat(u_train, 1, 2);
    
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
 
    rx = tx;
    
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
    
    r_training = r(1:nr_training_bits/nr_bits_per_symbol);
   
    % epsilon between -1/4Lt, 1/4Lt
     epsilon_A = 0.0007;  
     epsilon = 0.0008;
     L_t = (nr_training_bits/nr_bits_per_symbol)/2;
     y = exp(1i*2*pi*epsilon*(1:(2*L_t))).*r_training;
     
       
 
    % Compute variance of complex noise according to report.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(10/10);
%     disp(snr_point);
    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(y))+j*randn(size(y)));
    y = y + n;
    % computer epsilon_hat
     yt= y(1+L_t:end)*ctranspose(y(1:L_t));
     epsilon_h = (1/(2*pi*L_t))*(atan(imag(yt)/real(yt)));

    M(snr_point)= epsilon_h-epsilon_A;
    snr_t = var(qpsk(b_train))/sigma_sqr;
    sigma_t = 1/(4*pi^2*L_t^3*snr_t);
%     p=0.01;
     T(snr_point) = qfuncinv(p/2)*sqrt(sigma_t);
%      disp(T(snr_point))
    if abs(M(snr_point))>T(snr_point);
    m(snr_point) = m(snr_point) + 1;
    end

        % Next block.
        end
   % Next Eb/No value.
end

    % Next p

%     p_fa = m/nr_blocks;
%     plot(EbN0_db, p_fa ,'o-');
%     hold on;
%     p_t = p*ones(1,25); 
%     plot(EbN0_db,p_t);
p_md(i) = (nr_blocks-mean(m))/(nr_blocks);
   end
  plot(P,p_md);
hold on;
