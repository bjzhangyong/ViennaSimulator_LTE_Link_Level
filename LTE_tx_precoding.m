function [precode_y, UE, D, W, U] = LTE_tx_precoding(LTE_params, layer_x, UE, AtPorts, codebook_index, RI, CDD,rb_numbers)
% Precoding according to 3GPP TS 36.211-820, Section 6.3.4 page 37
% TODO: reorganize this huge tree so it more readable. Check whether
% something can be precalculated and stored 
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
d = 0;

if (AtPorts ~= 1)   % for 1 transmit antenna no precoding is used (precoding matrix = 1)
    if (UE.mode == 2) % transmit diversity, Section 6.3.4.3
        if (AtPorts == 2)
            c = length(layer_x(1,:));
            precode_y = zeros(2,2*c);   
            X =1/sqrt(2)*[1,0,1i,0;0,-1,0,1i;0,1,0,1i;1,0,-1i,0]*[real(layer_x(1,:));real(layer_x(2,:));...
                imag(layer_x(1,:));imag(layer_x(2,:))];
            precode_y(1,1:2:2*c-1)=X(1,:); 
            precode_y(2,1:2:2*c-1)=X(2,:);
            precode_y(1,2:2:2*c)=X(3,:); 
            precode_y(2,2:2:2*c)=X(4,:);
            D = 0;  % these matrices are needed in the receiver for undoing the precoding 
            W = 0;
            U = 0;
        else
            c = length(layer_x(1,:));
            precode_y = zeros(4,4*c); % NOTE: cleanup
            Z = [1,0,0,0,1i,0,0,0;0,0,0,0,0,0,0,0;0,-1,0,0,0,1i,0,0;0,0,0,0,0,0,0,0;0,1,0,0,0,1i,0,0;...
                0,0,0,0,0,0,0,0;1,0,0,0,-1i,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,1,0,0,0,1i,0;...
                0,0,0,0,0,0,0,0;0,0,0,-1,0,0,0,1i;0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,1i;0,0,0,0,0,0,0,0;...
                0,0,1,0,0,0,-1i,0];
            X =1/sqrt(2)*Z*[real(layer_x(1,:));real(layer_x(2,:));real(layer_x(3,:));real(layer_x(4,:));...
                imag(layer_x(1,:));imag(layer_x(2,:));imag(layer_x(3,:));imag(layer_x(4,:))];
            precode_y(1,1:4:4*c-3)=X(1,:); 
            precode_y(2,1:4:4*c-3)=X(2,:);
            precode_y(3,1:4:4*c-3)=X(3,:); 
            precode_y(4,1:4:4*c-3)=X(4,:);
            precode_y(1,2:4:4*c-2)=X(5,:); 
            precode_y(2,2:4:4*c-2)=X(6,:);
            precode_y(3,2:4:4*c-2)=X(7,:); 
            precode_y(4,2:4:4*c-2)=X(8,:);
            precode_y(1,3:4:4*c-1)=X(9,:); 
            precode_y(2,3:4:4*c-1)=X(10,:);
            precode_y(3,3:4:4*c-1)=X(11,:); 
            precode_y(4,3:4:4*c-1)=X(12,:);
            precode_y(1,4:4:4*c)=X(13,:); 
            precode_y(2,4:4:4*c)=X(14,:);
            precode_y(3,4:4:4*c)=X(15,:); 
            precode_y(4,4:4:4*c)=X(16,:);
            D = 0;
            W = 0;
            U = 0;
        end
    else    % spatial multiplexing
        
        if (AtPorts == 2)   % precoding matrix according to Table 6.3.4.2.3-1
            switch codebook_index
                case 0
                    if (RI == 2)
                        W=1/sqrt(2)*[1,0;0,1];
                    else
                        W = [1;0];
                    end
                case 1
                    if (RI == 2)
                        W=1/(2)*[1,1;1,-1];
                    else
                        W = [0;1];
                    end
                case 2
                    if (RI == 2)
                        W=1/(2)*[1,1;1i,-1i];
                    else
                        W = 1/sqrt(2)*[1;1];
                    end
                case 3
                    if (RI == 2)
                        error('codebook index not supported');
                    else
                        W = 1/sqrt(2)*[1;-1];
                    end
                case 4
                    if (RI == 2)
                        error('codebook index not supported');
                    else
                        W = 1/sqrt(2)*[1;1i];
                    end
                case 5
                    if (RI == 2)
                        error('codebook index not supported');
                    else
                        W = 1/sqrt(2)*[1;-1i];
                    end
                otherwise
                    error('codebook index not supported');
            end
        else % precoding matrix according to Table 6.3.4.2.3-2
            for i = 1:length(codebook_index)    % for open loop spatial multiplexing 
                                                % codebook_index = [12,13,14,15]
                                                % see 3GPP TS 36.213-820, Section 7.1.3 
                                                % NOTE: store this in some pre-stored variable
                switch codebook_index(i)
                    case 0
                        u= [1;-1;-1;-1];
                    case 1
                        u= [1;-1i;1;1i];
                    case 2
                        u= [1;1;-1;1];
                    case 3
                        u= [1;1i;1;-1i];
                    case 4
                        u= [1;(-1-1i)/sqrt(2); -1i;(1-1i)/sqrt(2)];
                    case 5
                        u= [1;(1-1i)/sqrt(2); 1i;(-1-1i)/sqrt(2)];
                    case 6
                        u= [1;(1+1i)/sqrt(2); -1i;(-1+1i)/sqrt(2)];
                    case 7
                        u= [1;(-1+1i)/sqrt(2); 1i;(1+1i)/sqrt(2)];
                    case 8
                        u= [1;-1;1;1];
                    case 9
                        u= [1;-1i;-1;-1i];
                    case 10
                        u= [1;1;1;-1];
                    case 11
                        u= [1;1i;-1;1i];
                    case 12
                        u= [1;-1;-1;1];
                    case 13
                        u= [1;-1;1;-1];
                    case 14
                        u= [1;1;-1;-1];
                    case 15
                        u= [1;1;1;1];
                    otherwise
                        error('codebook index not supported');
                end
                W1 = diag(ones(1,4))-2*u*u'/(u'*u);
                switch RI  
                    case 1
                        W = W1(:,1);
                    case 2
                        switch codebook_index(i)
                            case {0,4,5,9}
                                W(:,:,i) = W1(:,[1 4])/sqrt(2);
                            case {1,2,3,8,12,15}
                                W(:,:,i) = W1(:,[1 2])/sqrt(2);
                            otherwise
                                W(:,:,i) = W1(:,[1 3])/sqrt(2);
                        end
                    case 3
                        switch codebook_index(i)
                            case {0,4,5,8}
                                W(:,:,i) = W1(:,[1 2 4])/sqrt(3);
                            case {1,2,3,10,12,13,14,15}
                                W(:,:,i) = W1(:,[1 2 3])/sqrt(3);
                            otherwise
                                W(:,:,i) = W1(:,[1 3 4])/sqrt(3);
                        end
                    case 4
                        switch codebook_index(i)
                            case {0,1,4,5,8,9,12,15}
                                W(:,:,i) = W1(:,[1 2 3 4])/2;
                            case {6,7,10,11,13}
                                W(:,:,i) = W1(:,[1 3 2 4])/2;
                            otherwise
                                W(:,:,i) = W1(:,[3 2 1 4])/2;
                        end
                    otherwise
                        error('RI not supported');
                end
            end
        end
        switch CDD  % cyclic delay diversity, precoding according to Section 6.3.4.2
            case 0  % zero delay CDD
                D = diag(ones(1,AtPorts));  
                precode_y = D*W*layer_x;
                U = 0;
            case 1  % small delay CDD 
                eta = [128 256 512 1024 2048];
                index = find((eta-LTE_params.Nrb*LTE_params.Nsc*ones(1,5))>=0);
                n = eta(index(1));
                D = zeros(LTE_params.Nsc*length(rb_numbers),AtPorts);
                vec = [];
                for ii = 1:length(rb_numbers)
                    if(rb_numbers(ii) > LTE_params.Nrb)
                        vec =[vec;((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc).'];
                    else
                        vec =[vec;((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc).']; 
                    end
                end
                switch AtPorts 
                    case 2
                        d = 2/n;
                        D(:,1)=1;
                        D(:,2) = exp(-1i*2*pi*d*vec);
                    case 4
                        d = 1/n;
                        D(:,1)=1;
                        D(:,2) = exp(-1i*2*pi*d*vec);
                        D(:,3) = exp(-1i*2*pi*d*vec*2);
                        D(:,4) = exp(-1i*2*pi*d*vec*3);
                end
                precode_y = W*layer_x;
                U = 0;
            case 2  % large delay CDD
                switch RI
                    case 1
                        precode_y = W*layer_x;
                        D =0;
                    case 2
                        U = 1/sqrt(2)*[1,1;1,exp(-1i*pi)];
                        D = [1,0;0,exp(-1i*pi)];
                           if (AtPorts ==2)    % faster encoding (for 2 Antennaports W is not dependent on the time index)
                            precodey1= W*D*U*layer_x(:,1:2:end);
                            precodey2= W*D^2*U*layer_x(:,2:2:end);
                            ende = floor(length(layer_x)/2);
                            precode_y = reshape([precodey1(:,1:ende);precodey2(:,1:ende)],AtPorts,[]);
                            if (length(precodey1)>ende)
                                precode_y=[precode_y,precodey1(:,end)]; 
                            end
                           end
                    case 3
                        U = 1/sqrt(3)*[1,1,1;1,exp(-1i*2*pi/3),exp(-1i*4*pi/3);1,exp(-1i*4*pi/3),exp(-1i*8*pi/3)];
                        D = [1,0,0;0,exp(-1i*2*pi/3),0;0,0,exp(-1i*4*pi/3)];
                    case 4
                        U = 1/2*[1,1,1,1;1,exp(-1i*2*pi/4),exp(-1i*4*pi/4),exp(-1i*6*pi/4);...
                            1,exp(-1i*4*pi/4),exp(-1i*8*pi/4),exp(-1i*12*pi/4);...
                            1,exp(-1i*6*pi/4),exp(-1i*12*pi/4),exp(-1i*18*pi/4)];
                        D = [1,0,0,0;0,exp(-1i*2*pi/4),0,0;0,0,exp(-1i*4*pi/4),0;0,0,0,exp(-1i*6*pi/4)];
                    otherwise
                        error('RI not supported');
                end
                if (RI > 1 && AtPorts == 4) % slower encoding version (if W is time dependent)
                    l = 1:length(layer_x);
                    k = mod(floor((l-1)/RI),4)+1;
                    p = mod(l-1,RI)+1;
                    for i = 1:l(end)
                       precode_y(:,i) = W(:,:,k(i))*D^(p(i))*U*layer_x(:,i);
                    end
                end
            otherwise
                    error('CDD not supported');          
        end
    end
else        % single Antenna Port (SISO)
    precode_y = layer_x;
    UE.CDD = 0; % no CDD for single antenna transmission
    W=1;
    D=1;
    U = 1;
end
%figure(2)
%plot(precode_y,'x');