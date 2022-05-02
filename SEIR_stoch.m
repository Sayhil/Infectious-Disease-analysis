function result = SEIR_stoch (params, initial)
% Simulation of the stochastic SVEIR model

state = initial; % holds the state variables S, I, E, V and R

result.time = []; % receives the time results
result.S = [];    % receives the S results
result.I = [];    % receives the I results
result.V = [];    % receives the V results
result.E = [];    % recieves the E results
result.R = [];    % receives the R results

time = 0;
while ( state.I > 0)
        
    % calculate process probabilities for current state
    probs = probablities(state, params);

    % WHEN does the next process happen?
    tau = exprnd(1/sum(probs));
    
    %update time
    time = time + tau;
    
    % determine WHICH process happens after tau
    which = process(probs);
    
    % update state
    switch which
        case 1
            state.E = state.E - 1;
            state.I = state.I + 1;
        case 2
            state.E = state.E + 1;
            state.S = state.S - 1;
        case 3
            state.I = state.I - 1;
            state.R = state.R + 1;
        case 4
            state.V = state.V - 1;
            state.E = state.E + 1;
    end
    
    % store results
    result.time = [result.time time];
    result.S = [result.S state.S];
    result.V = [result.V state.V];
    result.I = [result.I state.I];
    result.E = [result.E state.E];
    result.R = [result.R state.R];
    
end

function which = process (probs)
% Determines which process happens

r = rand * sum(probs);
which = 1;
s = probs(1);
while (r > s)
    which = which + 1;
    s = s + probs(which);
end


function a = probablities (state, params)
% Calculates process probabilities for given state and parameters

a(1) = params.s * state.E; % sypmtoms
a(2) = params.b * state.S * state.I; % infection
a(3) = params.r * state.I; % recovery
a(4) = params.alpha * params.b * state.V * state.I; %infection of V


