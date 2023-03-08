function runCS3opt(ofun, F_VTR)     % F_VTR		"Value To Reach" (stop when ofunc < F_VTR)

switch ofun
    case 'STHX2D2'                  % example: shell-and-tube heat exchanger optimization from Saari et al(2019)
        varfile='STHX2varD2.txt';
        ofun='objfunSTHX2D2';
    case 'myFunction'               % add your own function here
        varfile='myVar.txt';
        ofun='objfun_myFunction';
end
 
% This bit here reads problem-related data;
%1) number of variables, 2) min.bounds, 3) max.bounds 
fid = fopen(varfile,'r');
I_D = fscanf(fid, '%g', 1);             % I_D: number of decision vars
FVr_minbound = fscanf(fid, '%g', I_D);  % vector of I_D values, min boundaries of each decision variable  
FVr_maxbound = fscanf(fid, '%g', I_D);  % vector of I_D values, max boundaries of each decision variable 
IVr_discr = fscanf(fid, '%g', I_D).';   % vector of I_D values, 1 for discrete, 0 for continuous decision variable
fclose(fid);

% Set up the struct for the optimizer with problem data. 
S_struct.F_VTR=F_VTR;
S_struct.I_D=I_D;
S_struct.FVr_minbound=FVr_minbound.';
S_struct.FVr_maxbound=FVr_maxbound.';
S_struct.IVr_discr=IVr_discr.';
        
% Set algorithm tuning and control parameters & outputs for the optimizer 
S_struct.I_NP = 30;         % I_NP > 5      population size; 2*I_D...6*I_D depending on difficulty seem work well
S_struct.F_CR = 0.1;        % 0<F_CR<1      crossover probability; F_CR ~ 0.1 for non-separable, high for separable.
S_struct.F_alpha = 0.1;     % F_alpha << 1  scaling factor. 0.01 ... 0.1, high values emphasizes global exploration 
S_struct.F_beta = 1.5;      % F_beta ~ 1.5  Lévy exponent, currently hard-coded to 1.5 in CS3.m --> value has no effect.
S_struct.I_itermax = 10000; % 10^2...6      maximum number of iterations (generations)
S_struct.I_bnd_constr = 1;  % 0/1           one if constraints are bound, 0 if 
S_struct.I_refresh = 2;     % >=0           output interval in generations
S_struct.I_plotting = 0;    % [0,1]         0: no plot, 1: print to file for post-processing 
S_struct.I_strategy = 1;    % I_strategy    differential mutation strategy, 1: rand/1/bin 2: local-to-best
S_struct.I_NFEmax = S_struct.I_itermax*S_struct.I_NP*2; % estimate, and print out, max function evaluations based...
fprintf(1,'\n Max number of function evaluations: %d \n',S_struct.I_NFEmax); % ...on pop.size & iterations allowed.


%********************************************************************
% Start of optimization
%********************************************************************

     
if (strcmp(ofun,'STHX2D2')==1)     
    % varfile='STHX2varD2.txt';
	ofun='objfunSTHX2D2';
    [FVr_x,S_y,I_NFE] = CS3opt(ofun,S_struct); 

    
elseif (strcmp(ofun,'myFunction')==1)
    % run your own function from here
    % could have separate varfile that's read here, overwriting lines 14-26.
    [FVr_x,S_y,I_NFE] = CS3opt(ofun,S_struct);
    
else                    
    [FVr_x,S_y,I_NFE] = CS3opt(ofun,S_struct);
    save unk_optres.mat FVr_x S_y I_nf  
end

%print out the result:
fprintf(1,'\n Optimization finished. Best found value in %d evaluations %5.0f\n ',I_NFE, S_y.FVr_oa(1));
fprintf(1,'Decision variable vector achieving best result:\n  ');
save(convfile);
for n=1:S_struct.I_D
    fprintf(1,'  %d:%6.4f ',n,FVr_x(n));
end


end



