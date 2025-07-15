function ds = Reduction_Entropy_CeO2(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent from
% Bulfin et al. (2016)
% Input: delta in [-]
% Output: s in [J/mol-K]
R = 8.31446261815324;   % Universal gas constant [J/mol-K]
% --- Choose only one option ---
% Option 1 from Bulfin et al. (2016) - Table 1
% delta_m = 1/2.31;
% delta_s_th = 165;       % Thermal entropy change [J/mol-K]
% a = 1;
% Option 2 from Bulfin et al. (2016) - Table 2
% delta_m = 1/2.9;
% delta_s_th = 160;       % Thermal entropy change [J/mol-K]
% a = 1.83;
% Calculate entropy of reduction
% ds = delta_s_th+(a/delta_m)*R*(log(delta_m-delta)-log(delta));
% Janna's fitting of Panlener's data
% log10delta = log10(delta);
% if log10delta <= -3.6
%     ds = 291.1;
% elseif log10delta <= -1.03
%     ds = 40.1831 - 60.8925*log10delta + 150.565*(log10delta^2) + ...
%          126.019*(log10delta^3) + 39.6964*(log10delta^4) + 4.47333*(log10delta^5);
% elseif log10delta <= -0.56
%     ds = -73.0801 - 726.348*log10delta - 859.650*(log10delta^2) - ...
%          368.734*(log10delta^3);
% else
%     ds = 129.6;
% end
% Vectorized - Option 1
% log10delta = log10(max(delta, 1e-300));  % Avoid log10(0) or negative input
% 
% ds = zeros(size(log10delta));
% 
% % Region 1: log10delta <= -3.6
% idx1 = log10delta <= -3.6;
% ds(idx1) = 291.1;
% 
% % Region 2: -3.6 < log10delta <= -1.03
% idx2 = log10delta > -3.6 & log10delta <= -1.03;
% x2 = log10delta(idx2);
% ds(idx2) = 40.1831 - 60.8925*x2 + 150.565*(x2.^2) + ...
%            126.019*(x2.^3) + 39.6964*(x2.^4) + 4.47333*(x2.^5);
% 
% % Region 3: -1.03 < log10delta <= -0.56
% idx3 = log10delta > -1.03 & log10delta <= -0.56;
% x3 = log10delta(idx3);
% ds(idx3) = -73.0801 - 726.348*x3 - 859.650*(x3.^2) - 368.734*(x3.^3);
% 
% % Region 4: log10delta > -0.56
% idx4 = log10delta > -0.56;
% ds(idx4) = 129.6;
% Vectorized - option 2
log10delta = log10(max(delta, 1e-300));
ds = zeros(size(log10delta));

x = log10delta;

idx1 = x <= -3.6;
idx2 = x > -3.6 & x <= -1.03;
idx3 = x > -1.03 & x <= -0.56;
idx4 = x > -0.56;

ds(idx1) = 291.1;

if any(idx2)
    x2 = x(idx2);
    x2_2 = x2 .* x2;
    x2_3 = x2_2 .* x2;
    x2_4 = x2_3 .* x2;
    x2_5 = x2_4 .* x2;
    ds(idx2) = 40.1831 - 60.8925*x2 + 150.565*x2_2 + ...
               126.019*x2_3 + 39.6964*x2_4 + 4.47333*x2_5;
end

if any(idx3)
    x3 = x(idx3);
    x3_2 = x3 .* x3;
    x3_3 = x3_2 .* x3;
    ds(idx3) = -73.0801 - 726.348*x3 - 859.650*x3_2 - 368.734*x3_3;
end

ds(idx4) = 129.6;
end