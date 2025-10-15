function [flow] = buildFlow(v_limit,dist,rho)
%BUILD FLOW - Flow struct build
%   Function building the flow struct used as input to execute the flutter
%   analysis of the wing
%
%   SYNTAX:
%       [flow] = buildFlow(v_limit,dist,rho)
%
%   INPUT:
%       v_limit,  double: vector containing min and max vactor to simulate
%       dist(*),  struct: two fileds defining velocity vector
%                           > type: logspace ('log') or linspace ('lin')
%                           > elem: number of velocity discrete value
%       rho(*),   double: air density to be considered
%
%   OUTPUT:
%       flow,  struct: flow propriety organized in fields:
%                       > v_min and v_max
%                       > rho
%                       > dist_type, velocity distribution type
%                       > v_vec, velocity vector values
%
%   OPTIONAL INPUT:
%       dist: set by default as 100 elements linear distribution
%       rho:  set by default at 1.225 kg/m3
%
%
    
    % Optional input
    if nargin < 3
        if nargin < 2
            dist.type = 'lin';
            dist.elem = 100;
        end
        rho = 1.225;
    end

    % Flow struct initialization
    flow = struct;
        flow.v_min = v_limit(1);
        flow.v_max = v_limit(2);
        flow.rho = rho;

    % Velocity distribution definition
    if isequal(dist.type,'lin')
        flow.v_vec  = linspace(v_limit(1),v_limit(2),dist.elem);
        flow.dist_type = 'linear';
    elseif  isequal(dist.type,'log')
        flow.v_vec = logspace(v_limit(1),v_limit(2),dist.elem);
        flow.dist_type = 'logarithmic';
    else
        error('Invalid distribution mode');
    end

end

