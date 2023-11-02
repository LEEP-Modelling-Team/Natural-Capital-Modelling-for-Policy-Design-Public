function [c, ceq] = mycon_Biod(p, q, costs, elm_option, cnst_data, cnst_target)

    % Calculate uptake
    uptake = myfun_uptake(p, q, costs, elm_option);

    % Biodiversity target constraints: cnst_target - spgrp_chg <0
    c   = [];
    ceq = [];

    num_spgrp = length(cnst_target);
    
    spgrp_chg = zeros(num_spgrp,1);
    for k = 1:num_spgrp
        spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
    end
    
    c = [c; cnst_target - spgrp_chg];
            
end