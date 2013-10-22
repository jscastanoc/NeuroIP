% Script for "perfect" dictionary

%% Init
clear; close all; clc;

verbose = true;


% Load data, define sample rate, number of eeg channels, etc...
load_data;

% Simulate brain activity and generate pseudo EEG
gen_eeg;

% Show simulated activity
if verbose
    nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),[]);
    pause(0.01)
end

% Depth bias compensation
model.L = nip_depthcomp(model.L,0.4);

%% Perfectly located dictionaries

fuzzy = nip_fuzzy_sources(model.cortex,3);
Nd = size(fuzzy,1);
for i = 1:length(actidx)
    for j = 1:3
        for k = 1:3
            if k == j
                D(:,k,i*j) = fuzzy(:,actidx(i)); % Gaussian blobs in the simualted active dipoles
            else
                D(:,k,i*j) = zeros(Nd,1);
            end
        end
    end
end
D = nip_trans_solution(D);

% Matrix with the perfectly located dictionaries mapped in to the lead
% field
L = model.L*D;
%% Compute hyperparameters and show results

h = nip_kalman_hyper(model.y,L);

if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    plot(model.t,h')
    title('Temporal evolution of the hyperparameters')
    pause(0.01)
end


%% Solve


