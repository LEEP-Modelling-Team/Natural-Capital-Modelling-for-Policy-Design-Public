% write_summary_results.m
% =======================

% 1. Initialise
% -------------
clear
rng(23112010)

% Model runs to document
% ----------------------
policies = { ...
        'oc_pay_1bill', ...
        'fr_act_pctl_rnd_1bill', ...
        'fr_act_1bill', ...
        'fr_act_1bill_15pct', ...
        'fr_act_1bill_15pct_pbio', ...
        'fr_env_1bill', ...
        'fr_env_1bill_15pct', ...
        'fr_env_1bill_15pct_pbio', ...
        'fr_es_1bill', ...
        'fr_es_1bill_15pct_pbio'};

% Model
% -----
mechanisms = {'fr_env', 'fr_es', 'fr_act', 'fr_act_pctl_rnd', 'oc_pay'};
urban_pct_limit = 0.5;
carbon_price_string = 'non_trade_central';
drop_vars = {'habitat_non_use', 'biodiversity'};
scheme_year = 1;

% Markup
% ------
markup = 1.15;

% Paths to Data Dir
% -----------------
data_folder  = 'D:\myGitHub\defra-elms\Data\';
data_path = [data_folder, 'elm_data_', carbon_price_string, '.mat'];

% Path to Solutions
% -----------------
soln_folder  = 'D:\myGitHub\Natural-Capital-Modelling-for-Policy-Design\Data\';

% Result Spreadsheet
% ------------------
xls_folder  = 'D:\myGitHub\Natural-Capital-Modelling-for-Policy-Design\Data\';
filename = [xls_folder 'results_summary.xlsx'];


% 2. Initialise
% -------------

% Load ELM option results from .mat file
% --------------------------------------
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

% Spreadsheet
% -----------

% (a) Title
% ---------
sheet = 'Summary';
xlswrite(filename, {'Summary'}, sheet, 'A1');

% (b) Column Headers
% ------------------
xlswrite(filename, {'Scheme'}, sheet, 'A3');
xlswrite(filename, {'Budget'}, sheet, 'B3');
xlswrite(filename, {'Bio Constraint'}, sheet, 'C3');
xlswrite(filename, {'Bio Pricing'}, sheet, 'D3');
xlswrite(filename, {'Benefits'}, sheet, 'E3');
xlswrite(filename, {'Spend'}, sheet, 'F3');
xlswrite(filename, {'Costs'}, sheet, 'G3');
xlswrite(filename, {'Ben:Spend'}, sheet, 'H3');
xlswrite(filename, {'Ben:Cost'}, sheet, 'I3');
xlswrite(filename, {'Uptake'}, sheet, 'J3');
xlswrite(filename, {'Surplus per Farm'}, sheet, 'K3');
xlswrite(filename, {'Surplus per Ha'}, sheet, 'L3');
xlsrow = 4;

% Loop through policies
% ---------------------
for i = 1:numel(policies)
    
    policy = policies{i};
    tokens = split(policy, '_');
    
    fprintf('\nPolicy: %s \n', policy);    
    fprintf('------------------\n');    

    % Payment mechanism for policy
    % ----------------------------
    for i = 1:numel(mechanisms)
        if contains(policy, mechanisms{i});
            mechanism = mechanisms{i};
        end
    end

    % Payment mechanism for policy
    % ----------------------------
    budget = tokens{contains(tokens, 'bill')};
    budget = str2double(erase(budget,'bill'))*1e9;

    % Policy has biodiversity constraint
    % ----------------------------------
    bio_constraint  = 0;  % 0 if no biodiversity constraint
    if any(contains(tokens, 'pct'))
        bio_constraint = tokens{contains(tokens, 'pct')};
        bio_constraint = str2double(erase(bio_constraint,'pct'));
        if isnan(bio_constraint) 
            bio_constraint = 0;
        else
            bio_constraint = bio_constraint/100;
        end
    end

    % Policy has biodiversity prices
    % ------------------------------
    bio_as_prices = false;  % only set to true if have a biodiversity const
    if bio_constraint > 0 && any(contains(tokens, 'pbio'))
         bio_as_prices = true; 
    end

    % Names of vars with quantity prices
    % ----------------------------------
    switch mechanism
        case 'fr_es'
            vars_price = vars_es_outs;
            quantities = es_outs;
        case 'fr_env'
            vars_price = vars_env_outs;
            quantities = env_outs;
        case {'fr_act', 'fr_act_pctl_rnd', 'oc_pay'}
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
            if ~any(strcmp(mechanism, {'fr_act', 'fr_act_pctl_rnd', 'oc_pay'}))
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
    
    cell_sample_ind = true(num_cells,1);


    % Load Solution for this Payment Mechanism
    % ----------------------------------------
    clear solution
    load([soln_folder 'solution_' policy '.mat'], 'solution');   

    % Number of cells chosen in this year
    num_chosen = size(solution.new2kid, 1);
    fprintf('Num Farmers with Agreement: %.0f \n', num_chosen);    

    % Identify Areas
    % --------------



    % XLS: Write Summary
    % ------------------
    sheet = 'Summary';

    % Summary Stats
    % -------------
    xlswrite(filename, {mechanism}, sheet, ['A' num2str(xlsrow)]);
    xlswrite(filename, budget, sheet, ['B' num2str(xlsrow)]);
    xlswrite(filename, bio_constraint, sheet, ['C' num2str(xlsrow)]);
    xlswrite(filename, bio_as_prices, sheet, ['D' num2str(xlsrow)]);
    xlswrite(filename, full(solution.fval), sheet, ['E' num2str(xlsrow)]);
    xlswrite(filename, solution.spend, sheet, ['F' num2str(xlsrow)]);
    xlswrite(filename, full(sum(solution.farm_costs)), sheet, ['G' num2str(xlsrow)]);
    xlswrite(filename, full(solution.fval)/solution.spend, sheet, ['H' num2str(xlsrow)]);
    xlswrite(filename, full(solution.fval)/sum(solution.farm_costs), sheet, ['I' num2str(xlsrow)]);
    xlswrite(filename, sum(solution.uptake_ind), sheet, ['J' num2str(xlsrow)]);
    surplus  = full(sum(solution.farm_payment-solution.farm_costs));
    xlswrite(filename, surplus/sum(solution.uptake_ind), sheet, ['K' num2str(xlsrow)]);
    xlswrite(filename, surplus/sum(solution.hectares), sheet, ['L' num2str(xlsrow)]);
    
    xlsrow = xlsrow + 1;

end

