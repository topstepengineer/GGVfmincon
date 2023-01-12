function [gLongMax, xopt] = combinedTraction(setup, accy, accx, Vx, x0)

%set the options for FMINCON
options = optimoptions("fmincon",...
    'Display','none',...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "ConstraintTolerance",1,...
    "MaxFunctionEvaluations",500,...
    "SubproblemAlgorithm","cg");

%set environmental conditions
g =             9.807;
rho_air =       1.22;

%set the number of min/max iterations to establish convergence
iterations = [2 20]; %[min max]

%get the static vertical loads
Fz0_fl = g*setup.mCar*setup.rWeightbal/2;
Fz0_fr = Fz0_fl;
Fz0_rl = g*setup.mCar*(1-setup.rWeightbal)/2;
Fz0_rr = Fz0_rl;
Fz0 = [Fz0_fl, Fz0_fr, Fz0_rl, Fz0_rr];

%set the lower and upper bounds for the slip angles
lb = [-0.2,-0.2, 0, 0];
ub = [0,0, 0.2, 0.2];

%preallocate some variables
gLongMax = zeros(1,length(Vx));

%initialize some values
n = 1;
setup.u = Vx/3.6;

%now the aero loads at this speed
FzAero_f = 0.5*rho_air*setup.Cz*setup.rAerobal*setup.u^2; %Front downforce
FzAero_r = 0.5*rho_air*setup.Cz*(1-setup.rAerobal)*setup.u^2; %Rear downforce
FxAero = 0.5*rho_air*setup.Cx*setup.u^2; %Drag force
FAero = [FzAero_f, FzAero_r, FxAero]; %combine aero forces into one vector

while n < iterations(2) %loop until max iterations are reached

    %get the vertical loads with the latest lateral/longitudinal accelerations
    [Fz_fl, Fz_fr, Fz_rl, Fz_rr] = GetFz(setup, Fz0, FAero, accy, accx);
    Fz = [Fz_fl,Fz_fr,Fz_rl,Fz_rr];

    %define some anonymous functions so we can pass extra variables
    fun = @(x)objective(x,setup,accy,Fz);
    con = @(x)constraint(x,setup,accy,Fz);

    %solve for the optimum slip angles and slip ratios
    [xopt,fval,exitflag,output] = fmincon(fun, x0, [],[],[],[],lb,ub,con,options);

    accx_old = accx;
    accx_new = (-fval-FAero(3))/setup.mCar;
    accx = (1/3)*accx_old + (2/3)*accx_new; % smooth the convergence...

    %save result
    gLong(n) = accx;

    if n > iterations(1)
        if abs(1-accx/gLong(n-1)) <= 0.001
            msg = ['gLong converged in ',num2str(n),' iterations to ', num2str(Vx), ',',num2str(accx)];
            disp(msg)
            disp(xopt)
            gLongMax = gLong(end);
%             [c, ceq] = constraint(xopt,setup,accy,Fz)
            break
        end
    end

    n = n + 1;
    if n >= iterations(2)
        disp('max iterations reached')
        plot(gLong)
    end

end

end


function obj = objective(x,setup,accy,Fz)
obj = -calcFxTires(x,setup,accy,Fz);
end

function [c, ceq] = constraint(x,setup,accy,Fz)
r = accy/setup.u; %yaw rate
beta = x(2) + setup.b*r/setup.u; %chassis slip angle
rwa = beta + setup.a*r/setup.u - x(1); %road wheel angle
Fpreload = 500; % differential preload
coeff_lock = 1.2; % differential locking coeff
[~, Fy_fl] = GetFxy(setup.tirParams_f,Fz(1),x(1),0,setup.camber_f,setup.u);
[~, Fy_fr] = GetFxy(setup.tirParams_f,Fz(2),-1*x(1),0,setup.camber_f,setup.u);
Fy_fr = Fy_fr*-1;
[Fx_rl, Fy_rl] = GetFxy(setup.tirParams_r,Fz(3),x(2),x(3),setup.camber_r,setup.u);
[Fx_rr, Fy_rr] = GetFxy(setup.tirParams_r,Fz(4),-1*x(2),x(4),setup.camber_r,setup.u);
Fy_rr = Fy_rr*-1;
c(1) = x(1)-x(2); %ensure and US balance
c(2) = (Fx_rr - Fx_rl)/2 - max(Fpreload/2,coeff_lock*(Fx_rl+Fx_rr)/2); %limited slip differential
ceq(1) = (Fy_fl+Fy_fr)*cos(rwa)*setup.a - (Fy_rl+Fy_rr)*setup.b + Fx_rl*setup.track/2 - Fx_rr*setup.track/2; %ensure zero yaw moment
ceq(2) = Fy_fl+Fy_fr+Fy_rl+Fy_rr-setup.mCar*accy; %Fy = mAy
end


function FxTires = calcFxTires(x,setup,accy,Fz)
r = accy/setup.u; %yaw rate
beta = x(2) + setup.b*r/setup.u; %chassis slip angle
rwa = beta + setup.a*r/setup.u - x(1); %road wheel angle
[Fx_fl, Fy_fl] = GetFxy(setup.tirParams_f,Fz(1),x(1),0,setup.camber_f,setup.u);
[Fx_fr, Fy_fr] = GetFxy(setup.tirParams_f,Fz(2),-1*x(1),0,setup.camber_f,setup.u);
Fy_fr = Fy_fr*-1;
[Fx_rl, ~] = GetFxy(setup.tirParams_r,Fz(3),x(2),x(3),setup.camber_r,setup.u);
[Fx_rr, ~] = GetFxy(setup.tirParams_r,Fz(4),-1*x(2),x(4),setup.camber_r,setup.u);
Fsd = (Fy_fl+Fy_fr)*cos(pi/2-rwa); %steering drag force
FxTires = Fx_fl+Fx_fr+Fx_rl+Fx_rr-Fsd;
end
