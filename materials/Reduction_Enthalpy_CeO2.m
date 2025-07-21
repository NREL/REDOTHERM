function dh = Reduction_Enthalpy_CeO2(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: h in [J/mol]
% Bulfin et al. (2015)
% dh = polyval([-64929 23368 1790 -1158 478],delta).*1e3;
% Bulfin et al. (2016)
% dh = (395-31.4*log10(delta)).*1e3;
% Ackermann et al. (2017)
% dh = 0.5.*(969.408715407529-503.738744939872.*delta.^0.5).*1e3;
% Janna's Panlener's data fitting
% if delta > 0
%     log10delta = log10(delta);
% else
%     log10delta = 1e-10;
% end
% 
% if log10delta <= -3.6
%     dh = 428.5;
% elseif log10delta <= -1.03
%     dh = 302.618 + 27.8543*log10delta + 253.263*(log10delta^2) + ...
%          184.585*(log10delta^3) + 54.0575*(log10delta^4) + 5.83517*(log10delta^5);
% elseif log10delta <= -0.56
%     dh = -643.250 - 5.62764e3*log10delta - 1.08066e4*(log10delta^2) - ...
%          8.77426e3*(log10delta^3) - 2.55864e3*(log10delta^4);
% else
%     dh = 408.5;
% end
% dh = dh*1e3;
% Vectorized - Option 1
% log10delta = log10(max(delta, 1e-300));  % Avoid log10(0) or negative input
% 
% dh = zeros(size(log10delta));
% 
% % Region 1: log10delta <= -3.6
% idx1 = log10delta <= -3.6;
% dh(idx1) = 428.5;
% 
% % Region 2: -3.6 < log10delta <= -1.03
% idx2 = log10delta > -3.6 & log10delta <= -1.03;
% x2 = log10delta(idx2);
% dh(idx2) = 302.618 + 27.8543*x2 + 253.263*(x2.^2) + ...
%            184.585*(x2.^3) + 54.0575*(x2.^4) + 5.83517*(x2.^5);
% 
% % Region 3: -1.03 < log10delta <= -0.56
% idx3 = log10delta > -1.03 & log10delta <= -0.56;
% x3 = log10delta(idx3);
% dh(idx3) = -643.250 - 5.62764e3*x3 - 1.08066e4*(x3.^2) - ...
%            8.77426e3*(x3.^3) - 2.55864e3*(x3.^4);
% 
% % Region 4: log10delta > -0.56
% idx4 = log10delta > -0.56;
% dh(idx4) = 408.5;
% Vectorized - Option 2
log10delta = log10(max(delta, 1e-300));
dh = zeros(size(log10delta));

% Precompute powers (optional)
x = log10delta;

idx1 = x <= -3.6;
idx2 = x > -3.6 & x <= -1.03;
idx3 = x > -1.03 & x <= -0.56;
idx4 = x > -0.56;

dh(idx1) = 428.5;

if any(idx2)
    x2 = x(idx2);
    x2_2 = x2 .* x2;
    x2_3 = x2_2 .* x2;
    x2_4 = x2_3 .* x2;
    x2_5 = x2_4 .* x2;
    dh(idx2) = 302.618 + 27.8543*x2 + 253.263*x2_2 + ...
               184.585*x2_3 + 54.0575*x2_4 + 5.83517*x2_5;
end

if any(idx3)
    x3 = x(idx3);
    x3_2 = x3 .* x3;
    x3_3 = x3_2 .* x3;
    x3_4 = x3_3 .* x3;
    dh(idx3) = -643.250 - 5.62764e3*x3 - 1.08066e4*x3_2 - ...
               8.77426e3*x3_3 - 2.55864e3*x3_4;
end

dh(idx4) = 408.5;

dh = dh.*1e3;

end