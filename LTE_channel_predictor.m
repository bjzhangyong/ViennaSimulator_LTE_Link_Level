function H_est = LTE_channel_predictor(prediction_buffer,delay,fading,predict)
% channel predictor based on 1D / 2D extrapolation 
% Author: Stefan Schwarz
% stefan.schwarz@nt.tuwien.ac.at
if predict
    numb = numel(find(sum(prediction_buffer(:,1,:,1,1),1)~=0));
%     surf(abs(reshape(prediction_buffer(:,:,1:numb,1,1),size(prediction_buffer,1),size(prediction_buffer,2)*numb)))
    switch fading
        case 'BlockFading'
            if sum(prediction_buffer(:,1,3,1,1)) == 0 || delay == 0 % predict only if there is some past data available
                H_est = reshape(prediction_buffer(:,:,end,:,:),size(prediction_buffer,1),size(prediction_buffer,2),size(prediction_buffer,4),size(prediction_buffer,5));
                return
            end
            for i1 = 1:size(prediction_buffer,4)  % linear prediction on every subcarrier
                for i2 = 1:size(prediction_buffer,5)
                    for i3 = 1:size(prediction_buffer,1)
                        data = reshape(prediction_buffer(i3,1,end-numb+1:end,i1,i2),1,numb);
                        H_temp = interp1(1:numb,data,numb+1:numb+delay,'linear','extrap');
                        H_est(i3,:,i1,i2)=repmat(H_temp,[1,size(prediction_buffer,2)]);
                    end
                end
            end
    %         for i1 = 1:size(prediction_buffer,4) % 2D prediction
    %             for i2 = 1:size(prediction_buffer,5)
    %                     data = reshape(prediction_buffer(:,1,end-numb+1:end,i1,i2),size(prediction_buffer,1),numb);
    %                     H_temp = griddata(1:numb,1:size(data,1),data,numb+1:numb+delay,1:size(data,1),'v4');
    %                     H_est(:,:,i1,i2)=repmat(H_temp,[1,size(prediction_buffer,2)]);
    %             end
    %         end
        case 'FastFading'
            if delay == 0 % predict only if there is some past data available
               H_est = reshape(prediction_buffer(:,:,end,:,:),size(prediction_buffer,1),size(prediction_buffer,2),size(prediction_buffer,4),size(prediction_buffer,5));
               return
            end
            for i1 = 1:size(prediction_buffer,4) % number of receive antennas   
                for i2 = 1:size(prediction_buffer,5)    % number of transmit antennas
                    for i3 = 1:size(prediction_buffer,1)    % number of subcarriers
                        data = reshape(prediction_buffer(i3,:,end-numb+1:end,i1,i2),1,size(prediction_buffer,2)*numb);
                        H_temp = interp1(1:size(data,2),data,size(prediction_buffer,2)*numb+1:size(prediction_buffer,2)*(numb+delay),'linear','extrap');
                        H_est(i3,:,i1,i2)=H_temp;
                    end
                end
            end  
    end
else
    H_est = reshape(prediction_buffer(:,:,end,:,:),size(prediction_buffer,1),size(prediction_buffer,2),size(prediction_buffer,4),size(prediction_buffer,5));
end
% surf(abs(H_est(:,:,1,1)));

