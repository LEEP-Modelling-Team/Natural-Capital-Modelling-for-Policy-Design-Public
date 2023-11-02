% Function to stop ga algorithm once found biodiversity constrained 
% price solution that is within the budget

function [state, options, optchanged] = mystop_budget(options, state, flag, q, c, elm_options, budget, cnst_data, cnst_target)
    
    % Options not changed in procedure
    % ---------------------------------
    optchanged = false;

    % Objective function (spend) at best current solution
    % ---------------------------------------------------
    ibest = state.Best(end);
    ibest = find(state.Score == ibest, 1, 'last');
    bestx = state.Population(ibest,:);
    bestf = myfun_spend(bestx, q, c, elm_options);

    % Stopping Criterion
    % ------------------
    % Are we within budget?
    if bestf < budget
        uptake = myfun_uptake(bestx, q, c, elm_options); 
        num_spgrp = length(cnst_target);
        spgrp_chg = zeros(num_spgrp,1);
        for k = 1:num_spgrp
            spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
        end
        % Do all biodiversity constraints hold?
        if all((spgrp_chg-cnst_target) > 0) 
           state.StopFlag = 'y';
        end
    end
end