function [Fx, Fy] = GetFxy(tirParams,Fz,alpha,kappa,gamma,Vx)

    nPoints = length(kappa);

    Fz = Fz.*ones(nPoints,1);
    alpha = alpha.*ones(nPoints,1);
%     kappa =     kappa.*ones(nPoints,1);
    gamma =     gamma.*ones(nPoints,1);
    Vx =        Vx.*ones(nPoints,1);
    phit =      0.*ones(nPoints,1);               % turnslip            	(1/m)

    %force the Fy values to zero at alpha = 0
    tirParams.LHX = 0;
    tirParams.LVX = 0;
    tirParams.LHY = 0;
    tirParams.LVY = 0;
    
    % Select a Use Mode
    useMode = 111;
    
    % Wrap all inputs in one matrix
    inputs = [Fz kappa alpha gamma phit Vx];

    % Store the output from mfeval in a 2D Matrix
    output = mfeval(tirParams, inputs, useMode);
        
    % Extract variables from output MFeval. For more info type "help mfeval"
    Fx = output(:,1);
    Fy = output(:,2);
end



