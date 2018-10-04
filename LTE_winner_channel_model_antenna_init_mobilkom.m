function Arrays = LTE_winner_channel_model_antenna_init_mobilkom(spacing,TX_type,RX_type)
% LTE winner channel model - to generate channel realization using Winner
% Model II [1]
%
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% [1]   IST-WINNER D1.1.2 P. Kyösti, et al., "WINNER II Channel Models", ver 1.1, Sept. 2007.
%       Available: https://www.ist-winner.org/WINNER2-Deliverables/D1.1.2v1.1.pdf
% [2]   TSG-RAN Working Group 4 (Radio) meeting #38 R4-060334
%
% input :
% output:   Arrays                      ... struct antenna specification using Winner II channel model
%
% date of creation: 2009/10/13
% last changes: 2009/10/13  Simko

global LTE_params;

% add Winner Model code path
path(path,'.\Winner Channel Model');
if exist('.\Winner Channel Model\dipole.m')==0
    error('Please download the Winner+ II channel model from http://projects.celtic-initiative.org/winner+/Documents/Publications/WIM2_3D_ant_ver064_220908.zip and extract it into this folder');
end


%% Base Station antennas

NAz=3*120; %3 degree sampling interval
Az=linspace(-180,180-1/NAz,NAz);
lambda = LTE_params.speed_of_light/LTE_params.carrier_freq;
% spacing = 1;  %in lambda

switch LTE_params.BS_config.nTx
    case 2
        switch TX_type
            case 'linear'
                 pattern(1,:,1,:) = dipole(Az,0);
                 BsArrays = AntennaArray('ULA',2,spacing*lambda,'FP-ECS',pattern,'Azimuth',Az); %ULA-2 
            case 'cross'
                Position = [0 0 0;0 0 0];
                Rotation = zeros(2,3);
                pattern = zeros(2,2,1,length(Az));
                pattern(1,:,1,:) = dipole(Az,45);
                pattern(2,:,1,:) = dipole(Az,-45);
                BsArrays = AntennaArray('Pos',Position,'Rot',Rotation,'FP-ECS',pattern);
            otherwise
                error('Something wrong')
        end
              

    case 4
        switch TX_type
            case 'linear'
                pattern(1,:,1,:) = dipole(Az,0);
                BsArrays = AntennaArray('ULA',4,spacing*lambda,'FP-ECS',pattern,'Azimuth',Az); %ULA-4 
            case 'cross'
                Position = [-1/2*spacing*lambda 0 0;-1/2*spacing*lambda 0 0;1/2*spacing*lambda 0 0;1/2*spacing*lambda 0 0];
                Rotation = zeros(4,3);
                pattern = zeros(4,2,1,length(Az));
                pattern(1,:,1,:) = dipole(Az,45);
                pattern(2,:,1,:) = dipole(Az,-45);
                pattern(3,:,1,:) = dipole(Az,45);
                pattern(4,:,1,:) = dipole(Az,-45);
                BsArrays = AntennaArray('Pos',Position,'Rot',Rotation,'FP-ECS',pattern);
            otherwise
                error('Something wrong')
        end
        % double cross-dipole with 45/-45 poralized dipols with lambda/2 spacing
        %BsArrays(4)=AntennaArray('ULA',4,dist,'FP-ECS',pattern,'Azimuth',Az); %UCA-4 1cm spacing
        
    otherwise
        error('Invalid number ot transmit antennas');
end


%% User antennas
clear pattern Rotation
NAz=120; %3 degree sampling interval
Az=linspace(-180,180-1/NAz,NAz);
switch LTE_params.UE_config.nRX
    case 2
        switch RX_type
            case 'linear'
                pattern(1,:,1,:) = dipole(Az,0);
                UserArrays = AntennaArray('ULA',2,1/10*lambda,'FP-ECS',pattern,'Azimuth',Az); %ULA-2 lambda/10 spaced
            case 'cross'
                Position = [0 0 0;0 0 0];
                Rotation = [0 0 0;0 0 0];
                pattern = zeros(2,2,1,length(Az));
                pattern(1,:,1,:) = dipole(Az,45);
                pattern(2,:,1,:) = dipole(Az,-45);
                UserArrays = AntennaArray('Pos',Position,'Rot',Rotation,'FP-ECS',pattern,'Azimuth',Az);
            otherwise
                error('Something wrong')
        end
    otherwise
        error('Invalid number ot receive antennas');
end

Arrays = [UserArrays,BsArrays];