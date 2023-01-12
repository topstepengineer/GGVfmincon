function [gLatMax, xopt] = pureGLat(setup, Vx, x0)

%set the options for FMINCON
options = optimoptions("fmincon",...
    'Display','off',...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "MaxFunctionEvaluations",500,...
    "SubproblemAlgorithm","cg");

%set environmental conditions
g =             9.807;
rho_air =       1.22;

%set the number of min/max iterations to establish convergence
iterations = [2 20]; %[min max]

%get the static vertical loads
Fz0_fl = g*setup.mCar*setup.rWeightbal/2;
Fz0_rl = g*setup.mCar*(1-setup.rWeightbal)/2;
Fz0_fr = Fz0_fl; %left and right are equal
Fz0_rr = Fz0_rl;
Fz0 = [Fz0_fl, Fz0_fr, Fz0_rl, Fz0_rr];

%initialize the lateral accelation
accy = 0;

%set the lower and upper bounds for the slip angles
lb = [-0.2,-0.2];
ub = [0,0];

%preallocate some variables
gLatMax = zeros(1,length(Vx));

%initialize some values
gLat = [];
n = 1;
setup.u = Vx/3.6;

%now the aero loads at this speed
FzAero_f = 0.5*rho_air*setup.Cz*setup.rAerobal*setup.u^2; %Front downforce
FzAero_r = 0.5*rho_air*setup.Cz*(1-setup.rAerobal)*setup.u^2; %Rear downforce
FxAero = 0.5*rho_air*setup.Cx*setup.u^2; %Drag force
FAero = [FzAero_f, FzAero_r, FxAero]; %combine aero forces into one vector

while n < iterations(2) %loop until max iterations are reached

    %get the vertical loads with the latest lateral/longitudinal accelerations
    [Fz_fl, Fz_fr, Fz_rl, Fz_rr] = GetFz(setup, Fz0, FAero, accy, 0);
    Fz = [Fz_fl,Fz_fr,Fz_rl,Fz_rr];

    %define some anonymous functions so we can pass extra variables
    fun = @(x)objective(x,setup,Fz);
    con = @(x)constraint(x,setup,Fz);

    %solve for the optimum slip angles
    [xopt,fval,exitflag,output] = fmincon(fun, x0, [],[],[],[],lb,ub,con,options);

    %calculate the yaw rates and steering angles
    r = accy/setup.u; %yaw rate
    beta = xopt(2) + setup.b*r/setup.u; %chassis slip angle
    rwa = beta + setup.a*r/setup.u - xopt(1); %road wheel angle
    [Fy_f, Fy_r] = calcFy(xopt,setup,Fz); %lateral forces at each axle
    Fy_f = Fy_f*cos(rwa); 
    accy = (Fy_f+Fy_r)/setup.mCar; %lateral acceleration result

    %save result
    gLat(n) = accy;

    if n > iterations(1)
        if abs(1-accy/gLat(n-1)) <= 0.001
            msg = ['gLat converged in ',num2str(n),' iterations to ', num2str(Vx), ',',num2str(accy)];
            disp(msg)
            disp(xopt)
            gLatMax = gLat(end);
            break
        end
    end
    n = n + 1;

    if n >= iterations(2)
        disp('max iterations reached')
    end

end
end

%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, ceq] = constraint(x,setup,Fz)
c = x(1)-x(2); %ensure and US balance
ceq =yawEquilibrium(x,setup,Fz); %ensure zero yaw moment
end

function obj = objective(x,setup,Fz)
obj = -calcFyTires(x,setup,Fz);
end

function yawEq = yawEquilibrium(x,setup,Fz)
[~, Fy_fl] = GetFxy(setup.tirParams_f,Fz(1),x(1),0,setup.camber_f,setup.u);
[~, Fy_fr] = GetFxy(setup.tirParams_f,Fz(2),-1*x(1),0,setup.camber_f,setup.u);
Fy_fr = Fy_fr*-1;
[~, Fy_rl] = GetFxy(setup.tirParams_r,Fz(3),x(2),0,setup.camber_r,setup.u);
[~, Fy_rr] = GetFxy(setup.tirParams_r,Fz(4),-1*x(2),0,setup.camber_r,setup.u);
Fy_rr = Fy_rr*-1;
yawEq = (Fy_fl+Fy_fr)*setup.a - (Fy_rl+Fy_rr)*setup.b;
end

function FyTires = calcFyTires(x,setup,Fz)
[~, Fy_fl] = GetFxy(setup.tirParams_f,Fz(1),x(1),0,setup.camber_f,setup.u);
[~, Fy_fr] = GetFxy(setup.tirParams_f,Fz(2),-1*x(1),0,setup.camber_f,setup.u);
Fy_fr = Fy_fr*-1;
[~, Fy_rl] = GetFxy(setup.tirParams_r,Fz(3),x(2),0,setup.camber_r,setup.u);
[~, Fy_rr] = GetFxy(setup.tirParams_r,Fz(4),-1*x(2),0,setup.camber_r,setup.u);
Fy_rr = Fy_rr*-1;
FyTires = (Fy_fl+Fy_fr+Fy_rl+Fy_rr);
end

function [Fy_f, Fy_r] = calcFy(x,setup,Fz)
[~, Fy_fl] = GetFxy(setup.tirParams_f,Fz(1),x(1),0,setup.camber_f,setup.u);
[~, Fy_fr] = GetFxy(setup.tirParams_f,Fz(2),-1*x(1),0,setup.camber_f,setup.u);
Fy_fr = Fy_fr*-1;
[~, Fy_rl] = GetFxy(setup.tirParams_r,Fz(3),x(2),0,setup.camber_r,setup.u);
[~, Fy_rr] = GetFxy(setup.tirParams_r,Fz(4),-1*x(2),0,setup.camber_r,setup.u);
Fy_rr = Fy_rr*-1;
Fy_f = (Fy_fl+Fy_fr);
Fy_r = (Fy_rl+Fy_rr);
end