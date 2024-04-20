%% Parameters to change
rn=0.004; % Nozzle radius [m]

Vw0=0.250e-3; % Initial water mass [m³]

p_i=8e5; % initial pressure [Pa]

me=0.1; % Empty weight [kg]

Vtot=0.5e-3; % Total bottle volume [m³]

t=linspace(0,3,4000); % Simulation time

drag=0.000;


%% Constants

g=9.8 % Gravitational acceleration [m/s²]
gamma = 7/5; % Adiabatic exponent of air.
penv=1e5; % Environmental pressure [Pa]
rho=1e3; % Water density [kg/m^3]

%% Calculations

An=pi*rn^2 % Nozzle area [m²]

m0=Vw0*rho; % Initial water mass [kg]

V0=Vtot-Vw0; % Initial air Volume [m³]

p0=p_i+penv; % Absolute initial pressure [Pa]

V=@(p) V0*(p0/p)^(1/gamma); % Air volume [m³]

Vw=@(p) Vtot-V(p); % Water volume [m³]

m=@(p) Vw(p)*rho; % Water mass [kg]

qm=@(p) (m(p)>0)*0.6*An*sqrt(2*rho*(p-penv)); % Water mass flow rate [kg/s]

pdot=@(p) -p0*V0^gamma * gamma * V(p)^(-gamma-1) * qm(p)/rho; 

%% 


F=@(p) qm(p)^2/rho/An;

a=@(p) F(p)/(me+m(p))-g;


pvy=lsode(@(pvy,t) [pdot(pvy(1)),a(pvy(1))-drag*sign(pvy(2))*pvy(2)^2/(me+m(pvy(1))),pvy(2)],[p0,0,0],t);

hold off;
plot(t,arrayfun(a,pvy(:,1))/10)
hold on;
plot(t,pvy(:,2:3));
plot(t,pvy(:,1)/1e5-1);

title(sprintf("Maximum Height: %.2f",max(pvy(:,3))))
legend("Acceleration [dam/s²]","Speed [m/s]","Height [m]","Pressure [bar]")
xlabel("Time [s]")
