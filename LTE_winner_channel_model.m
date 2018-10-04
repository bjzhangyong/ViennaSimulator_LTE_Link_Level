function [channel_out, delays, out] = LTE_winner_channel_model(N_subframes,Arrays,varargin)
% LTE winner channel model - to generate channel realization using Winner
% Model II [1]
%
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% [1]   IST-WINNER D1.1.2 P. Kyösti, et al., "WINNER II Channel Models", ver 1.1, Sept. 2007. 
%       Available: https://www.ist-winner.org/WINNER2-Deliverables/D1.1.2v1.1.pdf
%
% input :   N_subframes                 ... [1x1]   number of channel realization
%           Arrays                      ... struct  -> antenna specification
%           out                         ... struct  -> state of previous channel generation using winner model
% output:   channel                     ... [ x ] channel matrix
%           delays                      ... [ x ] delays matrix
%           out                         ... struct  -> output state of the winner model
%
% date of creation: 2009/10/12
% last changes: 2008/12/12  Simko

global LTE_params;

% add Winner Model code path
path(path,'.\Winner Channel Model');

%% winner II channel model

% set parameters {default option}
wimpar=wimparset;
wimpar.Scenario                 = ScenarioMapping(LTE_params.ChanMod_config.winner_settings.Scenario);
wimpar.PropagCondition          = LTE_params.ChanMod_config.winner_settings.PropagCondition;                % [LOS,{NLOS}]
switch LTE_params.ChanMod_config.filtering
    case 'BlockFading'
        Channel_Sampling_Time = LTE_params.Tsubframe;
    case 'FastFading'
        Channel_Sampling_Time = LTE_params.SamplingTime;
end
wimpar.SampleDensity            = LTE_params.speed_of_light/(2*LTE_params.carrier_freq*Channel_Sampling_Time*LTE_params.UE_config.user_speed); 
                                                                                                            % number of time samples per half wavelength [ {2} ]
wimpar.NumTimeSamples           = N_subframes;                                                              % number of time samples [ {100} ]
wimpar.UniformTimeSampling      = LTE_params.ChanMod_config.winner_settings.UniformTimeSampling;            % Use same time sampling grid for all links [ yes | {no} ] 
wimpar.FixedPdpUsed             = LTE_params.ChanMod_config.winner_settings.FixedPdpUsed;            	    % nonrandom path delays and powers [ yes | {no}]
wimpar.FixedAnglesUsed          = LTE_params.ChanMod_config.winner_settings.FixedAnglesUsed;                % nonrandom AoD/AoAs [ yes | {no} ]
wimpar.PolarisedArrays          = LTE_params.ChanMod_config.winner_settings.PolarisedArrays;                % usage of dual polarised arrays [ {yes} | no ]
wimpar.TimeEvolution            = LTE_params.ChanMod_config.winner_settings.TimeEvolution;                  % usage of time evolution  [ yes | {no} ]
wimpar.CenterFrequency          = LTE_params.carrier_freq;                                                  % carrier frequency in Herz [ {5.25e9} ]
wimpar.DelaySamplingInterval    = LTE_params.SamplingTime;                                                  % delay sampling grid [ {5e-9} ]LTE_params.SamplingTime
wimpar.PathLossModelUsed        = LTE_params.ChanMod_config.winner_settings.PathLossModelUsed;              % usage of path loss model [ yes | {no} ]
wimpar.ShadowingModelUsed       = LTE_params.ChanMod_config.winner_settings.ShadowingModelUsed;             % usage of shadow fading model [ yes | {no} ]
wimpar.PathLossModel            = LTE_params.ChanMod_config.winner_settings.PathLossModel;      	        % path loss model function name [ {pathloss} ]
wimpar.PathLossOption           = LTE_params.ChanMod_config.winner_settings.PathLossOption;                 % ['{CR_light}' or 'CR_heavy' or 'RR_light' or 'RR_heavy', CR = Corridor-Room, RR = Room-Room nlos} 
wimpar.RandomSeed               = LTE_params.ChanMod_config.winner_settings.RandomSeed;                     % sets random seed [ {[empty]} ]
wimpar.UseManualPropCondition   = LTE_params.ChanMod_config.winner_settings.UseManualPropCondition;         % whether to use manual propagation condition (los/nlos) setting or not. If not, the propagation condition is drawn from probabilities.  [ {yes} | no] 



% MsAAIdx = LTE_params.UE_config.nRX * ones(1,LTE_params.nUE);
% BsAAIdxCell = {[LTE_params.BS_config.nTx + 4]};
MsAAIdx = ones(1,LTE_params.nUE);   %every user is using antenna defined Arrays(1)
BsAAIdxCell = {[2]};    %   base station is using antenna defined Arrays(2)

layoutpar=layoutparset(MsAAIdx,BsAAIdxCell,LTE_params.nUE,Arrays);

layoutpar.ScenarioVector = LTE_params.ChanMod_config.winner_settings.Scenario*ones(1,LTE_params.nUE);       %   1=A1, 2=A2, 3=B1, 4=B2, 5=B3, 6=B4, 7=B5a, 8=B5c, 9=B5f, 10=C1,
                                                                                                            %   11=C2, 12=C3, 13=C4, 14=D1 and 15=D2a
                                                                                                            %   for more details look in  ScenarioMapping.mat 
switch LTE_params.ChanMod_config.winner_settings.PropagCondition
    case 'LOS'
        layoutpar.PropagConditionVector =1*ones(1,LTE_params.nUE);  %   (NLOS=0/LOS=1)
    case 'NLOS'
        layoutpar.PropagConditionVector =0*ones(1,LTE_params.nUE);  %   (NLOS=0/LOS=1)
end
for uu = 1:LTE_params.nUE
    layoutpar.Stations(1,uu+1).Velocity = [LTE_params.UE_config.user_speed;0;0];
end

optargin = size(varargin,2);
if optargin==1
    out = varargin{1};
    [channel, delays, out] = wim(wimpar,layoutpar,out);
elseif optargin==0
    [channel, delays, out] = wim(wimpar,layoutpar);
else
    error('Wrong number of input variables');
end
delays = round(delays/LTE_params.SamplingTime);             %   correct sampling in the delay domain


for user_i = 1:LTE_params.nUE
    channel{user_i}(isnan(channel{user_i})) = 0;
    channel_matrix_size = size(channel{user_i});
    channel_matrix_size(3) = max(delays(user_i,:))+1;
    channel_out{user_i} = zeros(channel_matrix_size);
    
    for tap_i = 1:channel_matrix_size(3)
        tap_positions = find(delays(user_i,:) == tap_i-1);
        if sum(tap_positions)>0
            channel_out{user_i}(:,:,tap_i,:) = sum(channel{user_i}(:,:,tap_positions,:),3);
        end
    end
    % channel normalization
    channel_energy = mean(sum(sum(sum(abs(channel_out{user_i}).^2,3),2),1),4);              %   mean channel energy = sum over transmit and receive antennas and over taps energies averaged over all realizations    
    channel_out{user_i} = sqrt(LTE_params.UE_config.nRX * LTE_params.BS_config.nTx) * channel_out{user_i} / sqrt(channel_energy);       % channel is normalized to have mean energy Nt*Nr
end
