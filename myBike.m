clear

w = warning ('off','all');

% Vx = linspace(70,250,5); %speed in kph
Vx = 250;

%Car geometries and weights
setup.hCoG =          0.5; %center of gravity height from ground
setup.track=          1.6; %track width
setup.wheelbase =     3.105; %wheelbase
setup.mUnsprung_f =    41; %unsprung mass of each front corner
setup.mUnsprung_r =    42; %unsprung mass of each rear corner
setup.Iw =            1000; %wheel rotational inertia
setup.Re =            0.3479; %baseline rolling radius, will be updated by Pacejka model after 1 iteration
setup.steering_ratio = 11;

%Car setup
frontTIRfile = 'MagicFormula52_Paramerters.tir';
rearTIRfile = 'MagicFormula52_Paramerters.tir';
setup.tirParams_f = mfeval.readTIR(frontTIRfile);
setup.tirParams_r = mfeval.readTIR(frontTIRfile);
setup.rMechbal =      0.55; %LLTD mechanical balance
setup.rWeightbal =    0.47; %weight balance
setup.rAerobal =      0.43; %aero balance
setup.Cz =            4.5; %downforce coefficient
setup.Cx =            1.5; %drag coefficient
setup.mCar =          1170; %complete car mass (uncluding unsprung elements)
setup.camber_f =      2.8/57.3; %front camber angle
setup.camber_r =      2/57.3; %rear camber angle
setup.brbal =         0.64; %brake balance forward percent
setup.splitfactor =   1; % 1=RWD, 0=FWD
setup.a = setup.wheelbase*(1-setup.rWeightbal);
setup.b = setup.wheelbase*setup.rWeightbal;

accy_sweep = linspace(0.1,0.9,10);
accx = 0;

for i = 1:length(Vx)
    %find the max accy for each Vx
    alpha_f_guess = -0.1;
    alpha_r_guess = -0.05;
    x0 = [alpha_f_guess alpha_r_guess];
    [accyMax(i), xopt] = pureGLat(setup, Vx(i), x0);

    %now we have the max accy
    alpha_f_guess = xopt(1);
    alpha_r_guess = xopt(2);
    kappa_rl_guess = 0.006;
    kappa_rr_guess = 0.001;
    x0 = [alpha_f_guess alpha_r_guess kappa_rl_guess kappa_rr_guess]; %combine the slip angles and ratios in one vector
    for j = 1:length(accy_sweep)
        [accxMax(i,j), xopt] = combinedTraction(setup, accyMax(i)*accy_sweep(j), accx, Vx(i), x0);
        accx = accxMax(i,j);
        x0 = [xopt(1) xopt(2) xopt(3) kappa_rr_guess];
    end

figure(2)
hold on
plot([accy_sweep.*accyMax(i) accyMax(i)],[accxMax(i,:) 0],'-o')
plot([-1.*accy_sweep.*accyMax(i) -1*accyMax(i)],[accxMax(i,:) 0],'-o')
end

