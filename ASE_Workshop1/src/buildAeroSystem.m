function [sys] = buildAeroSystem(sys,aero,c,rho,v_inf)
%BUILD AERO SYSTEM - Builds coupled system
%   Function building the coupled state space system which has as variables
%   the vector X = {q,q',xa}, where q are the general coordinates for
%   displacement, q' are the general coordinates for displacement speed and
%   xa are the aerodynamic coordinates in state space form.
%
%   SYNTAX:
%       [sys] = buildAeroSystem(sys,aero,c,rho,v_inf)
%
%   INPUT
%       sys,   struct: contains as fields M,K and C matrix of the dynamical system
%       aero,  struct: contains as fields the aerodynamic matrices
%       c,     double: chord of the wing
%       rho,   double: air density
%       v_inf, double: free stream velocity
%
%   OUTPUT:
%       sys,  struct: updated struct of the system containing as field:
%                     > V and A, state space matricies
%                     > dA_dU and dV_dU, derivatives of V and A wrt to speed
%
%

    q_inf = 0.5*rho*v_inf^2;    % compute dynamic pressure

    sys.V = [ eye(length(sys.M)),                        zeros(length(sys.M)),                               zeros(length(sys.M),length(aero.Aaero));
              zeros(length(sys.M)),                      sys.M-aero.D2aero .*q_inf.*(c./(2.*v_inf)).^2       zeros(length(sys.M),length(aero.Aaero));
              zeros(length(aero.Aaero),length(sys.M))    zeros(length(aero.Aaero),length(sys.M))             eye(length(aero.Aaero))];

    sys.A = [ zeros(length(sys.M)),                      eye(length(sys.M)),                                 zeros(length(sys.M),length(aero.Aaero));
              -sys.K+aero.D0aero.*q_inf                  -sys.C+aero.D1aero.*c.*q_inf./(2.*v_inf)            aero.Caero.*q_inf;
              aero.Baero.*(2.*v_inf)./c                  zeros(length(aero.Aaero),length(sys.M))             aero.Aaero.*(2.*v_inf)./c];

    sys.dA_dU = [  zeros(length(sys.M)),                 zeros(length(sys.M)),                               zeros(length(sys.M),length(aero.Aaero));
                   aero.D0aero.*rho*v_inf                aero.D1aero.*c.*rho/4                               aero.Caero.*v_inf*rho;
                   aero.Baero*2/c                        zeros(length(aero.Aaero),length(sys.M))             aero.Aaero*2/c];

    sys.dV_dU = zeros((length(sys.M)*2)+length(aero.Aaero), (length(sys.M)*2)+length(aero.Aaero));

end

