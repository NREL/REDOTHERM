function DeltaH_WS = Reaction_Enthalpy_WS(T)
% This function calculates the water splitting reaction enthalpy as a
% function of temperature
% Input:
% T - Temperature [K] (vector [1xN])
% Output:
% DeltaH_WS - water splitting reaction enthalpy [J/mol]
% Index:
% 1 - H2
% 2 - H2O
% 3 - O2
%%
T_mod = T./1000;    % Temperature divided by 1000 as needed in Shomate equations

flag = zeros(3,length(T));  % Flag to see if temperature is in low or high value range
A_LT = zeros(3,8);          % Initialize low temperature Shomate coefficients array
A_HT = zeros(3,8);          % Initialize high temperature Shomate coefficients array

flag(1,:) = ceil(fix(T./1000)./100);	% A flag to see if T>1000K
flag(2,:) = ceil(fix(T./1700)./100);	% A flag to see if T>1700K
flag(3,:) = ceil(fix(T./700)./100);     % A flag to see if T>700K

h_0 = [0 ; -241.826 ; 0];       % Standard enthalpy of formation [kJ/mol]

% A_H2 = [33.066178 -11.363417 11.432816 -2.772874 -0.158558 -9.980797 172.707974 0 ; ...
%     18.563083 12.257357 -2.859786 0.268238 1.977990 -1.147438 156.288133 0];
% A_H2O = [30.09200 6.832514 6.793435 -2.534480 0.082139 -250.8810 223.3967 -241.8264 ; ...
%     41.96426 8.622053 -1.499780 0.098119 -11.15764 -272.1797 219.7809 -241.8264];
% A_O2 = [31.32234 -20.23531 57.86644 -36.50624 -0.007374 -8.903471 246.7945 0 ; ...
%     30.03235 8.772972 -3.988133 0.788313 -0.741599 -11.32468 236.1663 0];

% Shomate equation coefficients (for 2 temperature ranges for each species
% LT - low temperature, HT - high temperature
A_LT(1,:) = [33.066178 -11.363417 11.432816 -2.772874 -0.158558 -9.980797 172.707974 0 ];
A_HT(1,:) = [18.563083 12.257357 -2.859786 0.268238 1.977990 -1.147438 156.288133 0];
A_LT(2,:) = [30.09200 6.832514 6.793435 -2.534480 0.082139 -250.8810 223.3967 -241.8264 ];
A_HT(2,:) = [41.96426 8.622053 -1.499780 0.098119 -11.15764 -272.1797 219.7809 -241.8264];
A_LT(3,:) = [31.32234 -20.23531 57.86644 -36.50624 -0.007374 -8.903471 246.7945 0 ];
A_HT(3,:) = [30.03235 8.772972 -3.988133 0.788313 -0.741599 -11.32468 236.1663 0];

%%

delta_H = (A_LT(:,1).*T_mod+(A_LT(:,2)/2).*T_mod.^2+(A_LT(:,3)/3).*T_mod.^3+(A_LT(:,4)/4).*T_mod.^4-A_LT(:,5)./T_mod+A_LT(:,6)-A_LT(:,8)).*(1-flag)+...
    (A_HT(:,1).*T_mod+(A_HT(:,2)/2).*T_mod.^2+(A_HT(:,3)/3).*T_mod.^3+(A_HT(:,4)/4).*T_mod.^4-A_HT(:,5)./T_mod+A_HT(:,6)-A_HT(:,8)).*flag;
delta_H = (delta_H+h_0).*1000;    % Convert from [kJ/mol] to [J/mol]

DeltaH_WS = delta_H(1,:)+0.5*delta_H(3,:)-delta_H(2,:);

end