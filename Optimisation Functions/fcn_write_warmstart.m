% Write Warmstart File for Cplex Optimiser
% ----------------------------------------
%
% Purpose
% -------
%   Takes a (full or partial) solution to a MILP and writes this to a
%   'warmstart' file to be read subsequently by Cplex as a starting
%   solution.
%
% Inputs
% ------
%   sln      - numeric vector of solution values. Can be a solution for
%              every choice variable or only certain choice variables. If
%              the latter, Cplex searches for a full solution using the
%              partial solution as a starting point.
%   idx      - integer vector of same length as sln giving the index in the
%              set of choice variables of each value passed in sln. The
%              first choice variable is index '0', the second '1' etc.
%   filename - path and filename to location where warmstart file will be
%              written. Filename should preferably have extension .mst
%   probname - name of the cplex problem, but uncertain this is important 
%              or has to match the name given to the cplex object to be
%              optimised.

function fcn_write_warmstart(sln, idx, filename, probname)
    
    % nvars = size(solutions,1);
    nslns = size(sln,2);
    
    % Open file for writing
    fid = fopen(filename, 'w');
    
    fprintf(fid, '<?xml version = "1.0" encoding="UTF-8" standalone="yes"?>\n');
    fprintf(fid, '<CPLEXSolutions version="1.2">\n');
    for k = 1:nslns
        fprintf(fid, ' <CPLEXSolution version="1.2">\n');
        fprintf(fid, '  <header problemName="%s" solutionsutionName="m%u" solutionsutionIndex="%u" MIPStartEffortLevel="4" writeLevel="1"/>\n', probname, k, k);           
        fprintf(fid, '  <variables>\n');
        fprintf(fid, '   <variable index="%u" value="%f"/>\n', [idx full(sln(:,k))]');
        fprintf(fid, '  </variables>\n');
        fprintf(fid, ' </CPLEXSolution>\n');
    end   
    fprintf(fid, '</CPLEXSolutions>\n');

end
