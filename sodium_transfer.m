%% Some simple model of sodium transfer 
% no real PDE solver, just a simple X(N+1) = X(N) + dX(N)
% updating values based on functions that depend on current values
% fps controls how many times per second it evaluates

% The model has only outer sphere, it is split in two parts: the part that
% interfaces nitrogen and the part that is in contact with sodium. It also
% has the lid and the bottom part of the sphere (~3m2 each) these parts
% considered perfect isolators and do not participate in heat exchange.

% Pressure inside 3M is hardcoded to 5 psig, that leads to ~2 L/s flow
% Heat exchange is happening by adding 75K nitrogen (calculated in the
% assumption of the adiabatic expansion of 2500psi 300K nitrogen into 20psi
% with gamma = 7/5, or in other words - we can not get any worse)

% Heat exchange N2 and the wall of 3M, heat exchange between N2 and Na, 
% and heat exchange between Na and the 3M walls. 
% Taken convective heat exchange model:
% q = h_i*(T2-T1)
% q - energy flux density [W/m2], h_i - convective heat exchange coeff [W/m2/K]
% T2 and T1 are temperatures of the medias [K]
% h_n2 ~ 2, h_na ~50 

% the model DOES NOT include thermal expansion of the nitrogen, it just
% assumes the density ~1.6 kg/m3 for the set up pressure 20 psi, 
% or in other words - we can not get any worse, the thermal expancion
% should increase pressure inside so we wouldn't need to add so much
% nitrogen hence less "cooling" will be in the real life.

% Perevalov A. Dec 2019

%%
% evaluations per second
fps = 5;

% physical parameters
r=1.46;                 % outer sphere redius (hardcoded in other functions)
V0 = 4/3*pi*r^3;        % outer sphere volume

% initial volume of the N2
V_n2_0 = 0.070;         % [m3]      initial N2 volume in 3M      
Tn2_0 = 125+273;        % [K]       initial temperature of N2

T_cool_n2 = 75;         % from gas_temp [K] 
Cp_n2 = 1039;           %(J/(kg K))  heat capacity of nitrogen
Cp_na = 1380;           %(J/(kg K))  heat capacity of sodium

% 3M walls 
T_wall = 125+273;       % temperature of the walls

% initial sodium
rho_na = 930;           % kg/m3 density of sodium
V_na = V0-V_n2_0;       % initial volume of sodium
Tna0 = 125+273;         % initial temperature of sodium


%% here we go
step = 0;

t_vector = [0];                         % time vector to save
V_vector = [V_na];                      % sodium volume array
h_vector = [sodium_height(V_na)];       % height on the axis
s_vector = [n2_surface(h_vector(1))];   % interface of N2 and the wall
T_N2_vector=[Tn2_0];                    % temperature of N2 in the 3M
m_N2_vector=[V_n2_0*1.5];               % mass of the N2
T_Na_vector = [Tna0];                   % temperature of sodium
power_vector = [0, 0, 0, 0];            % powers vector
dT_vector = [0];                        % temperature drop in the pipe 

while V_na > 0.1 && step < 10^5
    
    % checking if sodium is still ok
    if T_Na_vector(end) < 105
        warning('sodium is too cool, stopping now')
        break
    end
    
    step = step+1;                      % adding another step
    time = t_vector(end)+1/fps;         % writing time
    t_vector = [t_vector;time];         
    
    % pressure inside as a control parameter
    % ____________________________________________________________
    
    pn2 = 20;   % lets set it so the flux is ~2 L/s
    % ____________________________________________________________
    % density of the nitrogen inside
    rho_n2 = 1.2*pn2/15; % kg/m3
    
    %
    flux = -flux_pressure(pn2-15);      % evaluating the flux L/s
    flux = flux/1000;                   % switching to SI
    
    dV = flux/fps;                      % volume change per step
    V_na = V_na + dV;                   % evaluating the new volume
    V_vector = [V_vector; V_na];        % writing the new volume
    
    h = sodium_height(V_na);            % the height from the sodium level to the top of the sphere
    h_vector = [h_vector; h];
    
    s_3m_n2 = n2_surface(h);            % evaluating surface of 3M to N2
    s_vector = [s_vector; s_3m_n2];
    
    s_na_n2 = interface_n2_na(h);
    
    % 
    dm = -dV*rho_n2;                    % N2 mass change per step
    m = m_N2_vector(end);               % total mass of N2 inside
    m_N2_vector = [m_N2_vector; m + dm];
        
    T_n2=T_N2_vector(end);              % getting the last written temperature of N2 inside
    T_na=T_Na_vector(end);              % getting the last written T_Na
    
    %% now energy exchange
    % Convective Heat Transfer Coefficients [W/m2/K]
    h_n2 = 2;                           % 2 seems like an average number for N2
    h_na = 50;                          % the lowest I found for a liquid to metal is 50
    
    %% N2 heat exchange
    
    interface_N2_wall = max(0,s_3m_n2-n2_surface(0.3));
    dE1 = h_n2*interface_N2_wall*(T_wall-T_n2)/fps;   % exchange between the wall and n2
    dE2 = h_n2*s_na_n2*(T_na-T_n2)/fps;     % exchange between the Na and N2
    
    dEn2 = dE1+dE2;                         % adding dE's
   
    T_n2_next = (Cp_n2*(m*T_n2+dm*T_cool_n2)+dEn2)/(m+dm)/Cp_n2;  % heat exch eq
    T_N2_vector = [T_N2_vector; T_n2_next];
    
    %% now to sodium heat exchange
    dE3 = - dE2;                            % cooling due heating N2
    % evaluating efffective wall to Na heating area
    interface_na_wall_eff = max(0,(4*pi*r^2 - s_3m_n2) - n2_surface(0.3));
    
    dE4 = h_na* interface_na_wall_eff* (T_wall - T_na)/fps;  % heating sodium due the wall heat
    dEna = dE3 + dE4;
    
    T_na_next = T_na + dEna/Cp_na/V_na/rho_na;  % evaluating the next Na temperature
    T_Na_vector = [T_Na_vector; T_na_next];
    
    
    %% and about the pipe
    pipe_surf = pi*1.5*25.4/1000*h;
    power_loss_tube = pipe_surf*h_n2*(T_na_next-T_n2);
    
    power_vector = [power_vector; [dE1, dE2, dE4,power_loss_tube/fps]*fps];
    
    s_pipe = pi/4*(1.5*25.4/1000)^2;
    v_pipe = (-flux)/s_pipe;
    pipe_dT = 2*h*h_n2*(T_na-T_n2)/(Cp_na*rho_na*1.5*25.4/1000*v_pipe);
    
    dT_vector = [dT_vector; pipe_dT];
end


%% Plotting
figure(1)
plot(t_vector,V_vector/V0*100,'r',t_vector,h_vector/2/r*100,'b','LineWidth',2)
legend('Volume of sodium','sodium level','Location','south')
xlabel('Time, s')
ylabel('Normalized by maximum values, %')
title('The volume and the level of the sodium')
set(gca,'FontSize',15)


figure(2)
plot(t_vector,T_N2_vector-273,'r',t_vector,T_Na_vector-273,'b','LineWidth',2)
xlabel('Time, s')
ylabel('Temperature, C')
title('Temperatures of sodium and nitrogen in the sphere')
set(gca,'FontSize',15)
legend('T_{N2}','T_{Na}' ,'Location','southeast')

figure(3)
plot(t_vector,power_vector,'LineWidth',2)
legend('wall-N2','Na-N2','wall-Na','diptube','Location','north')
xlabel('Time, s')
ylabel('Power, W');
title('Heat exchange power')
set(gca,'FontSize',15)

figure(4)
plot(t_vector,dT_vector,'r','LineWidth',2)
xlabel('Time, s')
ylabel('Temperature drop in the diptube, K');
title('Temperature drop')
set(gca,'FontSize',15)