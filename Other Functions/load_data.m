function [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, vars_price, new2kid] = load_data(num_sample, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year)

    % (1) Set up
    %  ==========

    % Load ELM option results from .mat file
    % --------------------------------------
    % Generated in script2_run_elm_options.m
    % Depends on carbon price
    load(data_path);
    num_cells_all = length(cell_info.new2kid);
           
    % Remove cells 
    % ------------
    % Cells that are majority urban (excluding water area)
    cell_remove_ind = (cell_info.baseline_lcs.urban_ha./(cell_info.baseline_lcs.wood_ha + cell_info.baseline_lcs.farm_ha + cell_info.baseline_lcs.sngrass_ha + cell_info.baseline_lcs.urban_ha) ...
                            > urban_pct_limit); 
    
    % Cells where no farm land
    cell_remove_ind = or(cell_remove_ind, (cell_info.baseline_lcs.farm_ha < 1)); 
    
    % Cells where no land cost
    cell_remove_ind = or(cell_remove_ind, (costs.arable_reversion_sng_noaccess(:,1) + costs.destocking_sng_noaccess(:,1) == 0)); 
    
    % Remove Cells
    for k = 1:length(elm_options)
        elm_option_k = elm_options{k};
        benefits.(elm_option_k)             = benefits.(elm_option_k)(~cell_remove_ind,:);
        benefits_table.(elm_option_k)       = benefits_table.(elm_option_k)(~cell_remove_ind,:,:);
        benefit_cost_ratios.(elm_option_k)  = benefit_cost_ratios.(elm_option_k)(~cell_remove_ind,:);
        costs.(elm_option_k)                = costs.(elm_option_k)(~cell_remove_ind,:);
        costs_table.(elm_option_k)          = costs_table.(elm_option_k)(~cell_remove_ind,:,:);
        es_outs.(elm_option_k)              = es_outs.(elm_option_k)(~cell_remove_ind,:,:);
        env_outs.(elm_option_k)             = env_outs.(elm_option_k)(~cell_remove_ind,:,:);
        elm_ha.(elm_option_k)               = elm_ha.(elm_option_k)(~cell_remove_ind);
        biodiversity_constraints.(elm_option_k).data_20 = biodiversity_constraints.(elm_option_k).data_20(~cell_remove_ind, :);
        biodiversity_constraints.(elm_option_k).data_30 = biodiversity_constraints.(elm_option_k).data_30(~cell_remove_ind, :);
        biodiversity_constraints.(elm_option_k).data_40 = biodiversity_constraints.(elm_option_k).data_40(~cell_remove_ind, :);
        biodiversity_constraints.(elm_option_k).data_50 = biodiversity_constraints.(elm_option_k).data_50(~cell_remove_ind, :);
    end  
    new2kid = cell_info.new2kid(~cell_remove_ind);
    num_cells = length(new2kid);
    
    % Choose sample of cells to go into price search
    % ----------------------------------------------   
    if isnumeric(num_sample)
        cell_perm = randperm(num_cells);
        cell_sample_ind = (cell_perm <= num_sample)';
        budget = unscaled_budget * (num_sample/num_cells_all);
        num_cells = num_sample;
    elseif strcmp(num_sample, 'no')
        cell_sample_ind = true(num_cells,1);
        num_sample = num_cells;
        budget = unscaled_budget;
    else
        fprintf('''sample_num'' can only assume a numeric value or ''no'' if no sampling is required\n'); 
    end

    
    % Names of vars with quantity prices
    % ----------------------------------
    switch payment_mechanism
        case 'fr_es'
            vars_price = vars_es_outs;
            quantities = es_outs;
        case 'fr_env'
            vars_price = vars_env_outs;
            quantities = env_outs;
        case {'fr_act', 'fr_act_pctl', 'fr_act_pctl_rnd', 'oc_pay', 'up_auc'}
            vars_price = elm_options;
            quantities = elm_ha;  
    end
    if bio_as_prices 
       vars_price = [vars_price biodiversity_constraints.names_grp'];
    end  
    
    
    % Remove vars in 'drop_vars'
    % -------------------------- 
    if ~isempty(drop_vars)
        
        for var = drop_vars

            % Remove from Quantities
            % ----------------------
            if ~any(strcmp(payment_mechanism, {'fr_act', 'fr_act_pctl', 'fr_act_pctl_rnd', 'oc_pay', 'up_auc'}))
                [indvar, idxvar] = ismember(var, vars_price);            
                if indvar
                    % Remove from var list
                    vars_price(idxvar) = [];                                
                    % Remove from quantities                
                    for k = 1:length(elm_options)                 
                        quantities.(elm_options{k})(:,idxvar,:) = [];
                    end
                end
            end
            
            % Remove from Benefits
            % --------------------
            [indvar, idxvar] = ismember(var, vars_benefits);            
            if indvar   
                % Remove from benefits var list
                vars_benefits(idxvar) = [];                   
                % Remove from quantities
                for k = 1:length(elm_options) 
                    % Subtract away benefits for this var
                    benefits_drop = squeeze(benefits_table.(elm_options{k})(:, idxvar, :));                    
                    benefits.(elm_options{k}) = benefits.(elm_options{k}) - benefits_drop;
                    % Remove from benefits_table
                    benefits_table.(elm_options{k})(:,idxvar,:) = [];
                    % Recalculate benefits/costs ratio
                    benefit_cost_ratios.(elm_options{k}) = benefits.(elm_options{k}) ./ costs.(elm_options{k});                     
                end
            end                
        end        
    end
        
    % (2) Extract Sample Data
    % =======================
    
    % Select scheme year data
    % ------------------------
    % Note: All data has a 40 year time series for each cell for scheme  
    %       costs & benefits when implemented in each of those years.
    %       Those are expressed in npv in the base_year in £price_year.
    %       Model currently extracting just single year's data given by 
    %       data_year.
    benefits_year = nan(length(cell_sample_ind), num_elm_options);
    costs_year    = nan(length(cell_sample_ind), num_elm_options);
    for k = 1:num_elm_options
        benefits_year(:, k) = benefits.(elm_options{k})(:, data_year);
        costs_year(:, k)    = costs.(elm_options{k})(:, data_year);
        if ~any(strcmp(payment_mechanism, {'fr_act', 'fr_act_pctl', 'fr_act_pctl_rnd', 'oc_pay', 'up_auc'})) % ha quantities do not change across years of scheme
            quantities.(elm_options{k}) = quantities.(elm_options{k})(:, :, data_year);
        end
    end

    % Reduce data to cell sample
    % --------------------------
    % Extract relevant rows from above arrays/structures
    b = benefits_year(cell_sample_ind, :);    
    c = costs_year(cell_sample_ind, :);
    c = c .* markup; 
    for k = 1:num_elm_options
        hectares.(elm_options{k})       = elm_ha.(elm_options{k})(cell_sample_ind, :);
        quantities.(elm_options{k})     = quantities.(elm_options{k})(cell_sample_ind, :);
        bio_quantities.(elm_options{k}) = biodiversity_constraints.(elm_options{k}).data_20(cell_sample_ind, :);
    end 
    new2kid   = new2kid(cell_sample_ind);
    num_cells = length(new2kid);
    lu_data   = ones(num_cells, num_elm_options);
    
    % Construct data matrices by cell or by parcel
    % --------------------------------------------       
    if byparcel
        b = cell2pcl(b, elm_options);
        c = cell2pcl(c, elm_options); 
        for k = 1:num_elm_options
            if contains(elm_options{k}, 'arable')
                hectares.(elm_options{k})       = [hectares.(elm_options{k});       zeros(size(hectares.(elm_options{k})))];
                quantities.(elm_options{k})     = [quantities.(elm_options{k});     zeros(size(quantities.(elm_options{k})))];
                bio_quantities.(elm_options{k}) = [bio_quantities.(elm_options{k}); zeros(size(bio_quantities.(elm_options{k})))];
            else
                hectares.(elm_options{k})       = [zeros(size(hectares.(elm_options{k})));       hectares.(elm_options{k})];
                quantities.(elm_options{k})     = [zeros(size(quantities.(elm_options{k})));     quantities.(elm_options{k})];
                bio_quantities.(elm_options{k}) = [zeros(size(bio_quantities.(elm_options{k}))); bio_quantities.(elm_options{k})];
            end
        end        
        lu_data   = cell2pcl(lu_data, elm_options);
        new2kid   = repmat(new2kid,2,1);
        num_cells = length(new2kid);
    end
    
    % Area
    % ----
    hectares  = struct2table(hectares);

    % Quantities: Including add biodiversity as a quantity to be priced
    % -----------------------------------------------------------------       
    if any(strcmp(payment_mechanism, {'fr_act', 'fr_act_pctl', 'fr_act_pctl_rnd', 'oc_pay', 'up_auc'}))    
        if bio_as_prices 
            % Need to go to wide format for activity payments if adding
            % biodiversity as prices
            q = [];
            for k = 1:num_elm_options
                ha_option     = quantities.(elm_options{k});
                q_option      = zeros(num_cells, num_elm_options);
                q_option(:,k) = ha_option;
                q = cat(3, q, [q_option,bio_quantities.(elm_options{k})]);
            end
        else
            q = table2array(struct2table(quantities));    
        end    
    else        
        q = [];
        for k = 1:num_elm_options
            if bio_as_prices 
                quantities.(elm_options{k}) = [quantities.(elm_options{k}),bio_quantities.(elm_options{k})];
            end
            q = cat(3, q, quantities.(elm_options{k}));
        end                    
    end
    
    
    % Constraint Data
    % ---------------
    cnst_data = nan(num_cells, num_elm_options);
    for k = 1:num_elm_options
        cnst_data(:, :, k) = bio_quantities.(elm_options{k});
    end    
    cnst_data = permute(cnst_data, [2,3,1]);   % reorientate so rows: sp_grp, cols: option, depth: cells 
    % Target is for an absolute increase in national biodiversity defined 
    % as a percentage of current levels where that percent is in
    % bio_constraint. If bio_constraint = 0 then all cnst_target go to zero.
    % Constraints are scaled to the size of the sample.
    cnst_target = biodiversity_constraints.targets_20/num_cells_all * num_sample;
    cnst_target = cnst_target * bio_constraint;   

    
end

function pcldata = cell2pcl(celldata, elm_options)

    num_cells = length(celldata);
    num_options = length(elm_options);
    arable_options = find(contains(elm_options, 'arable'));
    grass_options  = find(contains(elm_options, 'destocking'));   
    
    pcldata = zeros(num_cells*2, num_options);
    pcldata(1:num_cells,arable_options)    = celldata(:, arable_options); 
    pcldata(num_cells+1:end,grass_options) = celldata(:, grass_options); 

end