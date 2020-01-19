% ============ Strain driven cyclic loading ================= %
%  by Mingzhao Zhuo Jan 2020, @TU Delft
%  copyright reserved
% =========================================================== %

clear;

% -------------------------------------------------------------
% material parameters
% -------------------------------------------------------------
Ea = 25e9;   % young's modulus
Em = 19.4e9; % unit Pa

k = 6.8e6; 
SNL = 0.043; % transformation strain
ds0 = k * SNL; % unit Pa/K

% transformation temperature 
Mf = 261; Ms = 262.5; As = 281; Af = 283;
du0 = (Ms + Mf + As + Af) / 4 * ds0;
w = (Ms - Mf + Af - As) / 4 * ds0;
c = 3.2e6; % heat capacity (Pa/K)
alpha = 11e-6; % expansion coefficient (/K)

% -------------------------------------------------------------
% manipulating parameter \lambda
% ratio of two time scale
% -------------------------------------------------------------
mu = 1.0e3; 

% -------------------------------------------------------------
% preallocate time vector
% -------------------------------------------------------------
% number of cycles
N = 1; 

dt = 0.0001;
t = 0 : dt : 2 * N;

% -------------------------------------------------------------
% prescribe the applied strain 
% -------------------------------------------------------------
% maximum strain
SN0 = 0.055; 
SN = SN0 / 2 * (1 - cos(pi * t));

% -------------------------------------------------------------
% preallocate other variables 
% -------------------------------------------------------------
len = length(SN);

% xi is the martensite phase fraction
xi = zeros(1, len); 
ximin = xi(1); ximax = xi(1);

% T0 is the initial temperature of the sample
T0 = 298; 
T = ones(1, len) * T0;

% SS is stress
SS = zeros(1, len); 

% beta is the controllable driving force
beta = ones(1, 1) * (SS(1) * SNL - T(1) * ds0); 

% -------------------------------------------------------------
% Cycling
% -------------------------------------------------------------
for i = 1:length(SN) - 1

    % Assume elastic deformation

    xi(i + 1) = xi(i);

    T(i + 1) = (-alpha / c * T(i) / (xi(i + 1) / Em + (1 - xi(i + 1)) / Ea) * (SN(i + 1) - SN(i)) - mu * (T(i) - T0) * dt) ...
        / (1 - alpha^2 / c * T(i) / (xi(i + 1) / Em + (1 - xi(i + 1)) / Ea)) + T(i);

    SS(i + 1) = (SN(i + 1) - xi(i + 1) * SNL - alpha * (T(i + 1) - T0)) / (xi(i + 1) / Em + (1 - xi(i + 1)) / Ea);

    beta(i + 1) = SS(i + 1) * SNL - T(i + 1) * ds0;

    %%  Judge if phase transition
    if beta(i + 1) > beta(i)

        if beta(i + 1) >= ximin * (Ms - Mf) * ds0 - Ms * ds0 && beta(i) <- Mf * ds0

            while 1 
                xirecord = xi(i + 1); 
                SSrecord = SS(i + 1);

                % Note that xi(i+1) must be placed before SS(i+1) otherwise divergency happens
                xi(i + 1) = (SN(i + 1) - SS(i + 1) * (xi(i + 1) / Em + (1 - xi(i + 1)) / Ea) - alpha * (T(i + 1) - T0)) / SNL;

                SS(i + 1) = (xi(i + 1) * (Ms - Mf) * ds0 - (Ms - T(i + 1)) * ds0) / SNL;

                T(i + 1) = -alpha / c * T(i) * (SS(i + 1) - SS(i)) + (SS(i + 1) * SNL + du0 - w * (2 * xi(i + 1) - 1)) / c * (xi(i + 1) - xi(i)) - mu * (T(i) - T0) * dt + T(i);

                if abs(xirecord - xi(i + 1)) < 1e-10 && abs((SSrecord - SS(i + 1)) / SS(i + 1)) < 1e-10
                    break
                end

            end

            beta(i + 1) = SS(i + 1) * SNL - T(i + 1) * ds0;
            ximax = xi(i + 1);
        end

    else

        if beta(i) >- Af * ds0 && beta(i + 1) <= ximax * (Af - As) * ds0 - Af * ds0

            while 1 
                xirecord = xi(i + 1); SSrecord = SS(i + 1);

                xi(i + 1) = (SN(i + 1) - SS(i + 1) * (xi(i + 1) / Em + (1 - xi(i + 1)) / Ea) - alpha * (T(i + 1) - T0)) / SNL;

                SS(i + 1) = (xi(i + 1) * (Af - As) * ds0 - (Af - T(i + 1)) * ds0) / SNL;

                T(i + 1) = -alpha / c * T(i) * (SS(i + 1) - SS(i)) + (SS(i + 1) * SNL + du0 - w * (2 * xi(i + 1) - 1)) / c * (xi(i + 1) - xi(i)) - mu * (T(i) - T0) * dt + T(i);

                if abs(xirecord - xi(i + 1)) < 1e-10 && abs((SSrecord - SS(i + 1)) / SS(i + 1)) < 1e-10
                    break
                end

            end

            beta(i + 1) = SS(i + 1) * SNL - T(i + 1) * ds0;
            ximin = xi(i + 1);
        end

    end

end

% -------------------------------------------------------------
% plot
% -------------------------------------------------------------
%%
figure
hold on
plot(SN, SS, '-', 'LineWidth', 2)