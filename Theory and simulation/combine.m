% Performance enhanced wireless device authentication using multiple
% weighted device-specific characteristics SIMULATION.
clear; 
clc;
%%
P_fa = linspace(0.001,0.1,10);
% N features, M realizations
N = 4; 
M = 100000;

p_rfa = zeros(1,length(P_fa));
p_rd = zeros(1,length(P_fa));
p_d = zeros(1, length(P_fa));
T = zeros(1, length(P_fa));
A = zeros(1,N);
a = zeros(1,N);
w = zeros(1,N);
k = 1;

for  j = 1:length(P_fa)
    
                n_a = 0;
                n_e = 0;

            % Create vectors
            X_v = ones(N, M);
            X_v(1,:)= 20.6;
            X_v(2,:)= 50;
            X_v(3,:)= 71.6;
            X_v(4,:)= 49.7;
%             X_v(5,:)= 59.7;
            
            sigmas = 0.01*ones(N,1);
            
%             for i = 2:N
%                 sigmas(i) = sigmas(i-1) * 0.0008;
%             end
            
            % Creat Noise (every row(i) have var_i= sigmas(i)) 
            sigmas = sigmas';
            R = randn(N,M);
            N_w = bsxfun(@times,R',sqrt(sigmas));
            N_n = N_w';
            % Alice or Eve?
%              X_r = 1.0003*X_v;  %Alice
              X_r = 1.055*X_v;  %Eve
            
           
            % Y estimated features
            Y = X_r + N_n;

            Delta_X = Y - X_v;
            X = X_r - X_v;
            
            % Compute weights
            for i = 1:N
                C = 50;
                m = mean(X,2);
                w(i)= m(i)/(sigmas(i)*C);
            end
            
            % Normalization weight
            if k == 1
            w = w/sum(w);
            else
            w = ones(1,N)/N;
            end   


            % Variable decision
              a = ones(N,M);
              for i = 1:N
              a(i,:) = m(i);
              end 
              A = a + R;
              S = (bsxfun(@times,A',w))';
              s = sum(S);
              
             % Thresholds
            T(j) = sqrt(w*w')*qfuncinv(P_fa(j));
            p_d(j) = qfunc(qfuncinv(P_fa(j))-mean(s)/sqrt(w*w'));
            
       for i = 1:M
            if (s(i) < T(j))
                n_a = n_a + 1;
            else
                n_e = n_e + 1;
            end
            
            % Next s 
       end
        p_rd(j) = n_e/(n_a+n_e);
end
        
%        if s == 0
%             p_rfa = n_e/(n_a+n_e);
%         else
%             p_rd = n_e/(n_a+n_e);
%        end
%        
        
        

% plot(P_fa, p_d,'o');

plot(P_fa, p_rd);
hold on;