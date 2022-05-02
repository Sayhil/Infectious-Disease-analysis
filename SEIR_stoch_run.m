clc
clearvars
close all

params.m = 1/(70*365); % host death rate
params.r = 1/14;  % rate of recovery
params.s = 1/6; %rate of incubation
params.alpha = 0.5; %percent of V population susceptible
k = 0.2; %percent of population vaccinated
aaa = [100	87	33.625	121.5	116.3333333	17.66666667	53.66666667	68.83333333	97.33333333	96.66666667	69.2	50.42857143	118	84.66666667	39	98.14285714	79.85714286	40.33333333	124	2];
bbb = 20;
for kk = 1:1:length(aaa)
    initial.V = round(aaa(kk)*0.12);  % number of Vaccinated individuals
    initial.I = 1;   % number of infected individuals
    initial.S = -initial.V-initial.I+aaa(kk);  % number of susceptible individuals
    initial.E = 0;   % number of infected individuals with no syptoms
    initial.R = 0;   % number of recovered individuals
    N = initial.S+initial.I+initial.E+initial.V+initial.R;  %total population
    R0 = 0.8;  %secondary infection
    % params.b = R0*(params.m+params.r)*(params.m+params.s)/N; % infection rate
    params.b = R0*params.r/((1-k)*N+params.alpha*k*N); %infection rate
    
    run_count = bbb; % number of runs
    
    result.time = []; % collects the time results
    result.S = [];    % collects the S results
    result.V = [];    % collects the V results
    result.I = [];    % collects the I results
    result.E = [];    % collects the E results
    result.R = [];    % collects the R results
    
    % simulate several stochastic SIR models and collect data
    for i=1:run_count
        out = SEIR_stoch (params, initial);
        result.time = [result.time out.time];
        result.S = [result.S out.S];
        result.V = [result.V out.V];
        result.I = [result.I out.I];
        result.E = [result.E out.E];
        result.R = [result.R out.R];
        
    % end
        [time, m, n] = unique(result.time);
        S = result.S(m);
        I = result.I(m);
        V = result.V(m);
        E = result.E(m);
        R = result.R(m);
    
    %     plot(time,I)
    %     hold on
    %     yyaxis right
    %     plot(time,V)
    %     hold on
    %     plot(time,S)
    %     hold on
    
        q(i) = S(length(S)); % total suceptible after infection dies out
    end
    
    % xlabel('time (days)')
    % ylabel('I')
    
    for i = 1:1:run_count
        percent_infection(i) = 100*(1-q(i)/initial.S);
        total_infection(i) = initial.S - q(i);
    end
    
    avg_infection(kk) = mean(percent_infection);
    % plot(percent_infection)
    % xlabel('Class')
    % ylabel('Percentage of infection')
end