function layer_x = LTE_common_layer_mapping(tx_mode,tx_user_symbols,nLayers,layer_x,nCodewords,cw,nAtPort)

switch tx_mode % transmission mode used for the UE
    case 1  % single antenna transmission, 3GPP TS 36.211-820, Section 6.3.3.1 page 36
        layer_x = tx_user_symbols;
    case 2  % transmit diversity, 3GPP TS 36.211-820, Section 6.3.3.3 page 37
        if (nLayers == 2)  % two layers
            layer_x(1,:)=tx_user_symbols(1:2:end);
            layer_x(2,:)=tx_user_symbols(2:2:end);
        else
            if (mod(length(tx_user_symbols),4)~=0)
                tx_user_symbols = [tx_user_symbols,0,0];
            end % see TS 36211-850 section 6.3.3.3
            layer_x(1,:)=tx_user_symbols(1:4:end);
            layer_x(2,:)=tx_user_symbols(2:4:end);
            layer_x(3,:)=tx_user_symbols(3:4:end);
            layer_x(4,:)=tx_user_symbols(4:4:end);
        end
    case {3,4,6} % spatial multiplexing, 3GPP TS 36.211-820, Section 6.3.3.2 page 36
        switch nLayers
                case 1
                    layer_x = tx_user_symbols;
                case 2
                    if (nCodewords == 2)
                        layer_x(cw,:)=tx_user_symbols;
                    else
                        if (nAtPort == 4)
                            layer_x(1,:)=tx_user_symbols(1:2:end);
                            layer_x(2,:)=tx_user_symbols(2:2:end);
                        else
                            error('number of layers not supported');
                        end
                    end
                case 3
                    if nCodewords == 1
                        layer_x(1,:) = tx_user_symbols(1:3:end);
                        layer_x(2,:) = tx_user_symbols(2:3:end);
                        layer_x(3,:) = tx_user_symbols(3:3:end);
                    else
                        if(cw == 1)
                            layer_x(1,:)=tx_user_symbols;
                        else
                            layer_x(2,:)=tx_user_symbols(1:2:end);
                            layer_x(3,:)=tx_user_symbols(2:2:end);
                        end
                    end
                case 4
                    if nCodewords == 1
                        layer_x(1,:) = tx_user_symbols(1:4:end);
                        layer_x(2,:) = tx_user_symbols(2:4:end);
                        layer_x(3,:) = tx_user_symbols(3:4:end);
                        layer_x(4,:) = tx_user_symbols(4:4:end);
                    else
                        layer_x(2*cw-1,:)=tx_user_symbols(1:2:end);
                        layer_x(2*cw,:)=tx_user_symbols(2:2:end);
                    end
                case 5
                    if cw == 1
                        layer_x(1,:)=tx_user_symbols(1:2:end);
                        layer_x(2,:)=tx_user_symbols(2:2:end);
                    else
                        layer_x(3,:) = tx_user_symbols(1:3:end);
                        layer_x(4,:) = tx_user_symbols(2:3:end);
                        layer_x(5,:) = tx_user_symbols(3:3:end);
                    end
                case 6
                    layer_x(3*cw-2,:) = tx_user_symbols(1:3:end);
                    layer_x(3*cw-1,:) = tx_user_symbols(2:3:end);
                    layer_x(3*cw,:) = tx_user_symbols(3:3:end);
                case 7
                    if cw == 1
                        layer_x(1,:) = tx_user_symbols(1:3:end);
                        layer_x(2,:) = tx_user_symbols(2:3:end);
                        layer_x(3,:) = tx_user_symbols(3:3:end);
                    else
                        layer_x(4,:) = tx_user_symbols(1:4:end);
                        layer_x(5,:) = tx_user_symbols(2:4:end);
                        layer_x(6,:) = tx_user_symbols(3:4:end);
                        layer_x(7,:) = tx_user_symbols(4:4:end);
                    end
                case 8
                    layer_x(4*cw-3,:) = tx_user_symbols(1:4:end);
                    layer_x(4*cw-2,:) = tx_user_symbols(2:4:end);
                    layer_x(4*cw-1,:) = tx_user_symbols(3:4:end);
                    layer_x(4*cw,:) = tx_user_symbols(4:4:end);              
                otherwise
                    error('Layer number not supported');
        end
end

end