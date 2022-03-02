function vortex_filt_pos = findVortexRegion(durationValid,WCenPos,pattStPos,Vx_flowmap_norm_phase,Vy_flowmap_norm_phase,curlz,paramsVor)
%%
tic
% Postive and negative porlarity is treated seperately: positve only (+ negative only)
expansion_threshold = paramsVor.expanPara;     %***** exapansion threshold: if data points count more in whole radias than filtered angle difference, break

% find angle difference between center vector (center of mass point to actual vector location) and actual vector

% load parcellation_template
% parcellation_template_01 = parcellation_template./parcellation_template;
        for irow = 1:size(curlz,1)
            for icol = 1:size(curlz,2)
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;
            end
        end
        
for ipatt = 1:size(WCenPos,2) %50; %**** select pattern ID
    % disp(['working on pattern ',num2str(ipatt)])
    for time_withinpatt = 1:durationValid(ipatt) %5;
        %temp1_time = absoluteTime{ipatt};
        %time = temp1_time(time_withinpatt);
        time = pattStPos(ipatt)+time_withinpatt-1 ;
        
        centerofmass_1patt_vxy_angle_dif = [];
        centerofmass_1patt_vxy_angle_dif_fil = [];
        
        
        
        temp1_vx = Vx_flowmap_norm_phase(:,:,time);
        temp1_vy = Vy_flowmap_norm_phase(:,:,time);
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy);
        
        temp1_centerofmass_1patt = WCenPos{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        
        vortex_filt = [];


        % x_tempalte = x_tempalte.*parcellation_template_01;
        % y_tempalte = y_tempalte.*parcellation_template_01;
        
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;
        centerofmass_1patt_vx_norm = centerofmass_1patt_vx./sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx + centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vy_norm = centerofmass_1patt_vy./sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx + centerofmass_1patt_vy.*centerofmass_1patt_vy);
        
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;
        
        for irow = 1:size(centerofmass_1patt_vxy_angle_dif,1)
            for icol = 1:size(centerofmass_1patt_vxy_angle_dif,2)
                temp1 = centerofmass_1patt_vxy_angle_dif(irow,icol);
                if temp1> pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) - 2*pi;
                elseif temp1< -pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) + 2*pi;
                    
                end
            end
        end
        
        % ******short cut, regardless polarity, assuming they won't go opposite directions, should've use polarity as indicator of angle diff
        centerofmass_1patt_vxy_angle_dif = abs(centerofmass_1patt_vxy_angle_dif);
        
        xlimit_low = round(centerofmass_1patt(1) - 10);
        xlimit_high = round(centerofmass_1patt(1) + 10);
        ylimit_low = round(centerofmass_1patt(2) - 10);
        ylimit_high = round(centerofmass_1patt(2) + 10);
        
        % angle difference filter
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2618; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2618; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;
        
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi;
            end
        end
        
        % radius filter
        % gradually increase radias until product sum equals to nan;
        for d = 2:0.5:100
            %         d
            centerofmass_1patt_vxy_abs_fil = centerofmass_1patt_vxy_abs; % radius filter
            centerofmass_1patt_vxy_abs_fil(centerofmass_1patt_vxy_abs_fil>d) = nan;
            temp1 = centerofmass_1patt_vxy_abs_fil./centerofmass_1patt_vxy_abs_fil;
            temp1_count_notnan = nansum(temp1(:));
            %     temp1_count_notnan
            
            temp2 = centerofmass_1patt_vxy_abs_fil .* centerofmass_1patt_vxy_angle_dif_fil; % combine 2 filters
            temp2_count_notnan = nansum(temp2(:)./temp2(:)); % count number of not nan data points
            
            %***** exapansion threshold: if data points count more in whole radias than filtered angle difference, break
            if temp1_count_notnan - temp2_count_notnan > expansion_threshold.*d ; % if data points count more in whole radias than filtered angle difference, break
                temp2_01 = temp2./temp2;
                temp2_01(isnan(temp2_01)) = 0;
                if centerofmass_1patt_round(1)>251 % remove the border centered vortex
                    break
                end
                if centerofmass_1patt_round(2)>176 % remove the border centered vortex
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) > 0
                    vortex_filt_pos{ipatt,time} = sparse(temp2_01.*2); % mark the polarity +1 if anticlock-wise
                elseif curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) < 0
                    vortex_filt_pos{ipatt,time} = sparse(-1.*temp2_01);    % mark the polarity -1 if clock-wise
                end
                break
            end
            
        end
    end
end
toc