function [delta] = addActuatorTerm(term,sweep,dihed,Phi,nodeLabel,actNode)
%ADD ACTUATOR TERM - Compute dK or dC contribute of the actuator
%   Function computing the contribute of a local added 'term' (stiffness or
%   damping) on the hinge actuator mode on the complete modal system of the
%   aircraft
%
%   SYNTAX:
%       [delta] = addActuatorTerm(term,sweep,dihed,Phi,nodeLabel,actNode)
%
%   INPUT:
%       term,       double: scalar value of physical K or C to add
%       sweep,      double: wing local reference system sweep angle in Rad
%       dihed,      double: wing local reference system sweep angle in Rad
%       Phi,        double: original modal shape of the aircraft provided by 'LabData_Session1.mat'
%       nodeLabel,  double: vector of all the nodes label of the system provided by 'LabData_Session1.mat'
%       actNode(*), double: vector of the two nodes associated to the hinge
%   
%   OUTPUT:
%       delta,  struct: struct containing the following fields
%                       > local: delta 'term' matrix in local reference
%                       > global: delta 'term' matrix in global reference
%                       > FEM: delta 'term' matrix in FEM reference
%                       > modal: delta 'term' matrix in original modal reference
%
%   OPTIONAL INPUT:
%       actNode: by default pick 1205 and 2205 as hinge nodes
%
%

    % Optional input
    if nargin < 6
        actNode = [1205, 2205];
    end

    % Initialize output strcut
    delta = struct;
        delta.local  = [];
        delta.global = [];
        delta.FEM    = [];
        delta.modal  = [];

    % Added 'term' local contribute
    delta.local = zeros(3);
    delta.local(2,2) = term;
    
    % Rotation tensor to the local wing reference
    R1 = [  cos(sweep), -sin(sweep),    0;
            sin(sweep),  cos(sweep),    0;
                     0,           0,    1];
    
    R2 = [  1,           0,            0;
            0,  cos(dihed),  -sin(dihed);
            0,  sin(dihed),   cos(dihed)];
    
    Rtot = R2' * R1';   % assembly total rotation tensor
    
    % Added 'term' global contribute
    delta.global = Rtot' * delta.local * Rtot;
    delta.FEM = [delta.global, -delta.global; -delta.global, delta.global];
    
    % Hinge node position
    hingeNode_Pos1 = find(actNode(1) == nodeLabel);
    hingeNode_Pos2 = find(actNode(2) == nodeLabel);
    
    % Hinge node modal base components
    hingeNode_Phi1 = squeeze(Phi(hingeNode_Pos1,4:6,:));
    hingeNode_Phi2 = squeeze(Phi(hingeNode_Pos2,4:6,:));
    
    % Added 'term' modal contribute
    delta.modal = [hingeNode_Phi1; hingeNode_Phi2]'* delta.FEM * [hingeNode_Phi1; hingeNode_Phi2];

end