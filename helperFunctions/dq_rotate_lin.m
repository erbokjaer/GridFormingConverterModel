function dx_c = dq_rotate_lin(xd, xq, ddelta, xd0, xq0, delta0)
%DQ_ROTATE_LIN Linearized dq rotation around operating point
%
% Inputs:
%   xd, xq     : actual signals
%   ddelta     : angle perturbation (delta - delta0)
%   xd0, xq0   : operating point
%   delta0     : operating angle
%
% Output:
%   dx_c = [d(xdc); d(xqc)]

    % Operating point trig
    c0 = cos(delta0);
    s0 = sin(delta0);

    % Form perturbations (CRITICAL FIX)
    dxd = xd - xd0;
    dxq = xq - xq0;

    % Linear state contribution
    d_xdc =  c0 .* dxd - s0 .* dxq;
    d_xqc =  s0 .* dxd + c0 .* dxq;

    % Angle contribution
    d_xdc = d_xdc + (-s0 .* xd0 - c0 .* xq0) .* ddelta;
    d_xqc = d_xqc + ( c0 .* xd0 - s0 .* xq0) .* ddelta;

    dx_c = [d_xdc + xd0; d_xqc + xq0];
end