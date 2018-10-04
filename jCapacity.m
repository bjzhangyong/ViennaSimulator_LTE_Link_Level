function varargout = jCapacity(Hf, varargin)
% 
% [CTxWF, CNoTx, PfTxWF, P, N] = jCapacity(Hf, P, N)
%
% Computes the capacity in bits/s/Hz (bps/Hz) for frequency-selective 
% wireless channels defined as a set of parallel channels in the frequency
% domain. The capacity is computed with channel knowledge at the
% transmitter side by using water filling power allocation over frequency
% (CTxWF) and without channel knowledge at the transmitter side (except 
% the knowledge of the SNR), thus the power is equally allocated over all
% frequencies (CNoTx).
%
% Input parameters:
%   Hf : attenuation values (channel coefficients). The dimension of Hf is
%     assumed like this:
%       size(hn) = [L, nRA, nTA]
%     where:
%       L is the number of frequency bins (number of parallel channels in the
%         frequency domain.
%       nRA = number of receive antennas.
%       nTA = number of transmit antennas.
%   P : sum power in linear units. The total transmit power. 
%     It is very important to note that in the case of capacity withouth
%     channel knowledge at the transmitter side all frequencies allocate
%     the same amount of power equal to: P / size(Hf, 1)
%   N : noise power in linear units.
%
% Output parameters:
%   CTxWF : resulting capacity in bps/Hz with channel knowledge at the
%     transmitter side (water filling capacity).
%   CNoTx : resulting capacity in bps/Hz without channel knowledge at the
%     transmitter side.
%   PfTxWF: transmit power per frequency bin, i.e., the water filling
%     coefficients. Restriction: sum(PfTxWF) <= P. The transmit power
%     coefficients for the case when no channel knowledge at the
%     transmitter is available are always equal to P / size(Hf, 1).
%   P,N can be recovered again.
%
% [CTxWF, CNoTx, PfTxWF, P, N] = jCapacity(Hf, SNR)
%
% The same as before but specifying the SNR instead of the power of
% the noise and the total transmit power.
%
% [CTxWF, CNoTx, PfTxWF, P, N] = jCapacity(Hf)
%
% Assumes SNR = 1, which means that the SNR is somehow included in the
% attenuation values instead of being normalized.
%

% vargargout initialisation
for index=1:nargout
    varargout{index} = -1;
end

if (nargin == 1)
    % No input parameter for the SNR. It means SNR == 1.
    P = 1;
    N = 1;
elseif (nargin == 2)
    % The input parameter is the SNR value.
    P = varargin{1};
    N = 1;
elseif (nargin == 3)
    % The input parameters are (by this order) the total transmit power and
    % the noise power, both in linear units.
    P = varargin{1};
    N = varargin{2};
else
    error(' Wrong input parameters. ');
end

L = size(Hf, 1);
nRA = size(Hf, 2);
nTA = size(Hf, 3);

% MIMO => SVD => spatial water filling
gs = zeros(L,min(nRA,nTA));

for frqIndex = 1:L,
    G = squeeze(Hf(frqIndex, :, :));
    singValues = svd(G);
    gs(frqIndex,:) = singValues;
end

Hf = reshape(gs,1,[]);



% The following code is valid for both SISO and MIMO calculations.
epsilon = 0.1;
Hf = abs(Hf).^2;
lr = min(N./Hf); % Cut-off value
Psum = 0;
while (Psum < P)
    % Cut-off value calculation
    lr = lr * (1 + epsilon);
    
    % Water filling coefficients calculation
    Pf = lr - (N ./ Hf);
    Pf(Pf < 0) = 0;

    % Total allocated power
    Psum = sum(Pf);
    
    % Unassigned power
    epsilon = abs(0.05*(P - Psum)/P);
    
% % % Code for debugging
% %     fprintf(' P = %2.2f, PTotal = %2.2f, epsilon = %2.5f\n', P, PTotal, epsilon);
    if (epsilon < 1e-5) 
        Psum = P; 
    end;
end

Pf_len = length(Pf);
varargout{1} = sum(log2(1 + ((Pf .* Hf) ./ N))) / L;
if (nargout >= 2)
    % Important!!!!!!
    % In order to not produce mutual information values greather than
    % capacity ones and given that the algorithm don't use some residual
    % power (that is not allocated using water filling), then the total
    % allocated power with water filling is used for the mutual information
    % calculation. 
    PfNoTx = sum(Pf) ./ Pf_len;
    varargout{2} = sum(log2(1 + (PfNoTx .* Hf) ./ N)) / L;
end
if (nargout >= 3)
    varargout{3} = Pf;
end
if (nargout >= 4)
    varargout{4} = P;
end
if (nargout >= 5)
    varargout{5} = N;
end
if (nargout > 5)
    error(' Wrong output parameters. ');
end

