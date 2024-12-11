%{
    UPDATE: CODE BLOCK WORK, HOWEVER,  PREDICTION ARE WRONG. INITIAL LAUNCH
    ANGLE CANNOT BE VARIED INDIVIDUALLY AS IT IS A DEPENDENT VARIABLE.

    DATE: 12/10/2024
    TIME STAMP: 10:30 PM MT
%}


%% Initialize test/verification data
clc; clear; close all;

%% Get Constants and Attach Changing Variables 

% attach variables 
Angle0 = linspace(0, 89, 100); % Initial Angle

% get constants
const = getConst();


%% Initial State Vector
statevector_0 = [const.x0, 0, const.z0, 0, const.mr0, const.Vair, const.Mair];

%% Time Span
tspan = [0 5];

%% Solve Combined Phases with ode45

% Create a matrix that accounts for the change in distance and height travelled with respect to the Initial Angle
dataAngle0_MaxDist = ones(length(Angle0), 1);
dataAngle0_MaxHeight= ones(length(Angle0), 1);

% Get data readings
for i = 1:length(Angle0)
    const.theta = Angle0(i);
    [~, state] = ode45(@(t, statevector) rocket_phases(t, statevector, const), tspan, statevector_0);
    dataAngle0_MaxDist(i) = max(state(:, 1));
    dataAngle0_MaxHeight(i) = max(state(:, 3));
end

%% Plot dataAngle0 vs VAngle0
figure;
plot(dataAngle0_MaxDist, Angle0, 'Color', '[0.4940, 0.1840, 0.5560]', 'LineWidth', 2, LineStyle='-')
hold on;

% Define x value at which the vertical line is needed
requiredDistance = 92;
% Interpolate the y-value (Angle) for x = 92
y_interp = interp1(dataAngle0_MaxDist, Angle0, requiredDistance, 'linear', 'extrap');

% Plot vertical line from x-axis to the interpolated point
plot([requiredDistance, requiredDistance], [0, y_interp], 'Color', '[0.4660, 0.6740, 0.1880]', 'LineWidth', 1, LineStyle='--');

% Plot horizontal line from y-axis to the interpolated point
plot([0, requiredDistance], [y_interp, y_interp], 'Color', '[0.4660, 0.6740, 0.1880]', 'LineWidth', 1, LineStyle='--');

% Highlight the intersection point
plot(requiredDistance, y_interp, 'Color', '[0.4660, 0.6740, 0.1880]', 'MarkerSize', 7, 'LineWidth', 1, 'MarkerFaceColor', '[0.6350, 0.0780, 0.1840]', Marker='o');

% Annotate the intersection point
text(requiredDistance + 1, y_interp, sprintf('(%.2f, %.2f)', requiredDistance, y_interp), 'Color', 'k');

% Display values that meet requirement
correspondingHeight = interp1(dataAngle0_MaxDist, dataAngle0_MaxHeight, requiredDistance, 'linear', 'extrap');
correspondingPressure0= y_interp;

% Display Results
fprintf("Required Max Distance: %2.1f ", requiredDistance);
fprintf("\n\nPredicted Corresponding Initial Angle: %2.3f ", correspondingPressure0);
fprintf("\nPredicted Corresponding Max Height: %2.1f ", correspondingHeight);

% plots
xlabel('Distance (m)');
ylabel('Initial Volume (degrees)');
title('Change in Distance Travelled vs Initial Angle');
xlim([-1 130])
ylim([-1 95])
xticks([0, 50 75 92 110])
yticks([0, 5, y_interp, 60, 90])
yline(0, '--', 'Color', 'blue', Label={'Minimum Angle'});
yline(y_interp, '--', 'Color', [0.4660, 0.6740, 0.1880], Label={'Extrapolated Value'});
yline(89, '--', 'Color', 'blue', Label={'Maximum Angle'});
grid on;
hold off;




%% Plot Trajectory and Thrust Profile With The Predicted Value of Initial Angle 

% Run ODE with predicted initial Angle 
const.theta = y_interp; % Angle prediction 
[t, state] = ode45(@(t, statevector) rocket_phases(t, statevector, const), tspan, statevector_0); 

Thrust = zeros(length(t), 1);
for i = 1:length(t)
    const.theta = y_interp; % change the initial Angle to the predicted value of Angle0
    [~, Thrust(i), ~] = rocket_phases(t(i), state(i, :), const);
end

Phase = zeros(length(t), 1);
for i = 1:length(t)
    [~, ~, Phase(i)] = rocket_phases(t(i), state(i, :), const);
end

% Get Phase times
t2_idx = find(Phase(:,1) == 2, 1);
t3_idx = find(Phase(:,1) == 3, 1);
distanceState = state(:, 1);
Phase2Time = t(t2_idx);
Phase3Time = t(t3_idx);
Phase2Dist = distanceState(t2_idx);
Phase3Dist = distanceState(t3_idx);

% Save Predictions
maxHeight = max(state(:, 3));
maxDistance = max(state(:, 1));

% Display Results
fprintf("\n\nAt prediction: ")
fprintf("\nActual Max Distance: %2.1f ", maxDistance);
fprintf("\nActual Max Height: %2.1f ", maxHeight);


% Trajectory Plot
figure;
plot(state(:, 1), state(:, 3), 'Color', '[0.4940, 0.1840, 0.5560]', 'LineWidth', 2, LineStyle='-');
xlabel('Distance (m)');
ylabel('Height (m)');
title('Rocket Trajectory');
xline(Phase2Dist, '--', 'Color', 'blue', Label={'Phase 2'})
xline(Phase3Dist, '--', 'Color', 'blue', Label={'Phase 3'})
xlim([0 95]);
ylim([0 30]);
grid on;
grid minor;

% Thrust vs Time Plot
figure;
plot(t, Thrust, 'Color', '[0.4940, 0.1840, 0.5560]', 'LineWidth', 2, LineStyle='-'); 
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Thrust vs Time');
xticks(0:0.05:0.2);
xline(Phase2Time, '--', 'Color', 'blue', Label={'Phase 2'})
xline(Phase3Time, '--', 'Color', 'blue', Label={'Phase 3'})
xlim([0 0.2]);
ylim([0 max(Thrust)]);
grid on;
grid minor;


%% Rocket Phases Function
function [d_statevector_dt, F_thrust, Phase] = rocket_phases(t, statevector, const)
 
    % Extract state variables
    x = statevector(1); % x-position
    vx = statevector(2); % x-velocity
    z = statevector(3); % z-position
    vz = statevector(4); % z-velocity
    mr = statevector(5); % Mass of rocket
    Vair = statevector(6); % Volume of air
    mAir = statevector(7); % Mass of air
    
    vel = [vx; vz];
    velMag = norm(vel);

    % Determine the heading vector
    if norm([x, z - const.z0]) > const.ls
        h_hat = vel/velMag;
    else
        h_hat = [cosd(const.theta); sind(const.theta)];
    end

    % Compute Drag Force
    D = 0.5 * const.rho_air * (velMag^2) * const.CD * const.A_body;
    
    % Set Pressure for Phase change
    if Vair > const.Vb
        pEnd = const.p0 * ((const.Vair / const.Vb)^const.gamma);
        p = pEnd * ((mAir / const.Mair)^const.gamma);
    end


    % Phase determination
    if Vair < const.Vb % Phase 1: Water exhaustion
        p = const.p0 * ((const.Vair / Vair)^const.gamma);
        F_thrust = 2*const.cdis*const.A_throat*(p - const.pa);
        mDot_r = -1*const.cdis*const.A_throat*sqrt(2*const.rho_w*(p - const.pa)); % change in rocket mass
        vDot_Air = const.cdis * const.A_throat *sqrt((2/const.rho_w)*(const.p0*((const.Vair/Vair)^const.gamma)-const.pa)); % change in air volume

        % Net forces
        Fx = (F_thrust * h_hat(1)) - (D * h_hat(1)); 
        Fz = (F_thrust * h_hat(2)) - (D * h_hat(2)) - (mr * const.g);

        % State derivatives
        ax = Fx / mr;
        az = Fz / mr;
        
        Phase = 1; % send phase instance 

        % Return derivatives
        d_statevector_dt = [vx; ax; vz; az; mDot_r; vDot_Air; 0];

    elseif p > const.pa % Phase 2: Air exhaustion

        rhoODE = mAir / const.Vb;
        Temp = p / (rhoODE * const.Rair);
        pCritical = p*(2/(const.gamma+1))^(const.gamma/(const.gamma-1)); % Critical Pressure 


        if pCritical > const.pa

            pExit = pCritical;

            TempExit = (2 / (const.gamma + 1)) * Temp;

            vExit = sqrt(const.gamma * const.Rair * TempExit);

            rhoExit = pExit / (const.Rair * TempExit);

        elseif pCritical < const.pa
            pExit = const.pa;

            MachExit = sqrt((2/(const.gamma-1))*(((p/const.pa)^((const.gamma-1)/const.gamma))-1));

            TempExit = Temp/(1+((const.gamma-1)/2)*MachExit^2);
            
            rhoExit = pExit / (const.Rair * TempExit);

            vExit = MachExit * sqrt(const.gamma * const.Rair * TempExit);
        end

        mDot_air = -1*const.cdis * rhoExit * const.A_throat * vExit; % change in air mass

        F_thrust = -1*mDot_air * vExit + (pExit - const.pa) * const.A_throat;

        % Net forces
        Fx = (F_thrust * h_hat(1)) - (D * h_hat(1));
        Fz = (F_thrust * h_hat(2)) - (D * h_hat(2)) - (mr * const.g);

        % Accelerations
        ax = Fx / mr;
        az = Fz / mr;

        Phase = 2; % send phase instance 


        % Return derivatives
        d_statevector_dt = [vx; ax; vz; az; mDot_air; 0; mDot_air];


    else % Phase 3: Ballistic flight

        % Thrust is zero
        F_thrust = 0;

        % Net forces
        Fx = -D * h_hat(1);
        Fz = -D * h_hat(2) - (mr * const.g);

        % Accelerations
        ax = Fx / mr;
        az = Fz / mr;

        if z < 0 && x > 0 
            vx = 0;
            vz = 0;
        end
        
        Phase = 3; % send phase instance 

        % Return derivatives
        d_statevector_dt = [vx; ax; vz; az; 0; 0; 0];

    end

end

%% Constants Function
function const = getConst()
    const.g = 9.81; % Gravity (m/s^2)
    const.cdis = 0.78; % Discharge coefficient
    const.rho_air = 0.961; % Air density (kg/m^3)
    const.Vb = 0.002; % Bottle volume (m^3)
    const.pa = 83426.563088; % Atmospheric pressure (Pa)
    const.gamma = 1.4; % Specific heat ratio
    const.rho_w = 1000; % Water density (kg/m^3)
    const.de = 0.021; % Nozzle diameter (m)
    const.dB = 0.105; % Bottle diameter (m)
    const.Rair = 287; % Specific gas constant for air (J/(kg·K))
    const.mB = 0.15; % Mass of empty bottle (kg)
    const.T0 = 310; % Initial air temperature (K)
    const.v0 = 0.0; % Initial velocity of rocket (m/s)
    const.theta = 40; % Launch angle (degrees)
    const.x0 = 0.0; % Initial x-position (m)
    const.z0 = 0.25; % Initial z-position (m)
    const.ls = 0.5; % Length of test stand (m)

    const.Vi_water = 0.0005; % Initial water volume (m^3)
    const.CD = 0.425; % Drag coefficient
    const.p0 = 330948.34944 + const.pa; % Initial air pressure (Pa)

    const.A_throat = pi * (const.de / 2)^2; % Nozzle throat area (m^2)
    const.A_body = pi * (const.dB / 2)^2; % Bottle cross-sectional area (m^2)
    const.Vair = const.Vb - const.Vi_water; % Initial volume of air
    const.Mair = (const.p0 * const.Vair) / (const.Rair * const.T0); % Initial mass of air
    const.mWater = const.rho_w * (const.Vb - const.Vair); % Initial mass of water
    const.mr0 = const.mB + (const.rho_w*(const.Vb - const.Vair)) + const.Vair*(const.p0/(const.Rair * const.T0)); % Initial mass of rocket
end