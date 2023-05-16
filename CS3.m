% -----------------------------------------------------------------
% Modified Cuckoo Search (CS3) algorithm 
% by Jussi Saari, C.Mendoza, J. Kaikko, E.Sermyagina, A.Mankonen. 
% based on original Cuckoo Search by Xin-She Yang and Suash Deb  
% Programmed by Saari et al. at LUT University, 2017-2023        
% Last revised: 2023                                        
% ----------------------------------------------------------------
% First published in:
% Saari J, Martinez CM, Kaikko J, Sermyagina E, Mankonen A, Vakkilainen E.
% Techno-economic optimization of a district heat condenser in a small 
% cogeneration plant with a novel greedy cuckoo search 
% Energy 2022, 239, 122622, ISSN 0360-5442
% doi.org/10.1016/j.energy.2021.122622.
%
% Saari J, Suikkanen H, Martinez CM, Hyvärinen J.
% Optimization of natural circulation district heating reactor 
% primary heat exchangers
% Energies 2023.
% -------------------------------------------------------------- %
%

function [FVr_bestnest,S_bestval,I_NFE] = CS3(ofun,S_struct)

% Input variables that MUST be available from 
I_NP            = S_struct.I_NP;            % population size
I_D             = S_struct.I_D;             % number of decision variables
FVr_minbound    = S_struct.FVr_minbound;    % decision var. minimum bounds
FVr_maxbound    = S_struct.FVr_maxbound;    % decision var. maximum bounds
F_VTR           = S_struct.F_VTR;           % termination crit, Value To Reach 
I_itermax       = S_struct.I_itermax;       % term.crit, max number of iterations

% Input variables checked from input struct 
if (isfield(S_struct,'I_plotmax'))          % maximum index for convergence output
    I_plotmax       = S_struct.I_plotmax;   % taken from input struct if available
else
    I_plotmax       = 100;                  % otherwise default to 100
end
if (isfield(S_struct,'I_refresh'))          % output refresh cycle
    I_refresh       = S_struct.I_refresh;   
else
    I_refresh       = floor(I_itermax/I_plotmax);
    I_refresh       = max(I_refresh,1);     % can't be < 1.
end
if (isfield(S_struct,'I_plotting'))         % are we plotting or not?
    I_plotting         	= S_struct.I_plotting;
else
    I_plotting          = 0;                % default no.
end
if (isfield(S_struct,'I_PopPrint'))         % are we printing out whole populations?
    I_PopPrint         	= S_struct.I_PopPrint;
else
    I_PopPrint          = 0;                % default to no.
end
if (isfield(S_struct,'I_NFEmax'))           % max.number of function evaluations
    I_NFEmax            = S_struct.I_NFEmax;
else
    I_NFEmax            = floor(I_itermax*I_NP*2);
end
if (isfield(S_struct,'F_pa1'))              % probability to "abandon nest" 
    F_pa1           = S_struct.F_pa1;       % i.e., to apply Lévy flights
else                                        % in FIRST generation
    F_pa1           = 0.2;                  % default 20%
end
if (isfield(S_struct,'F_paMax'))            % pa can be set to change linearly
    F_paMax         = S_struct.F_paMax;     % to F_paMax at the max iteration
else                                        % criteria. Default: constant pa,
    F_paMax      	= F_pa1;                % i.e. pa at end = pa at start.
end
if (isfield(S_struct,'F_CR'))
    F_CR             = S_struct.F_CR;
else
    S_struct.F_CR    = 0.1;
end

% alpha and mutation strategy are checked also, 
% but they gets sent around as part of S_struct, so...
if (isfield(S_struct,'F_alpha'))            
    %F_alpha         = S_struct.F_alpha;    % ..not copied into variable..
else
    S_struct.F_alpha = 0.05;                % ..but default saved to struct
end
if (isfield(S_struct,'I_strategy'))            
    %F_alpha         = S_struct.F_alpha;    
else
    S_struct.I_strategy = 1;                % default to rand/1/bin
end

fprintf(1,'\nStarting optimization.\nNFE(max) = %d \n\n',I_NFEmax);

prkl            = size(FVr_minbound);
FM_nestpop      = zeros(I_NP,prkl(2));
for i=1:I_NP                                                % Random initial population
    FM_nestpop(i,:) = FVr_minbound + (FVr_maxbound-FVr_minbound).*rand(size(FVr_minbound));
end

fitness         = (2^64)*ones(I_NP,1);                      % initialize kukkuupopulation's obj.f.values to smth very poor 
S_tempval       = feval(ofun,FM_nestpop(1,:),S_struct);     % evaluate first one ...
S_bestval       = S_tempval;                                % ...just to have one to be called the best
[fmin,FVr_bestnest,FM_nestpop,fitness,S_bestval,I_NFE] = get_best_nest(FM_nestpop,FM_nestpop,fitness,ofun,S_struct,S_bestval);  % Get the current best

I_RefreshCycle  = 0;
I_iter          = 0;
I_plotindex     = 0;

FM_convplot        = (10^10)*ones(2,I_plotmax);
FM_convplot(1,1:I_plotmax)=I_NFEmax+1;
FM_convplot(1,1)   = 1;
FM_convplot(2,1)   = S_bestval.FVr_oa(1);

fp_convhist = fopen('CS3Convhist.txt','w');        
fprintf(fp_convhist,'F_pa1: %f,  I_NP: %d\n\n',F_pa1,I_NP);
fclose(fp_convhist);
fp_pophist = fopen('CS3population_hist.txt','w');
fseek(fp_pophist, 0, 'eof');
fprintf(fp_pophist,'Population history.\nF_pa1: %f,  I_NP: %d\n\n',F_pa1,I_NP);
fprintf(fp_pophist,'Iteration %d\n',I_iter);  
for i=1:I_NP
    fprintf(fp_pophist,'%3.0f : %8.0f ',i, fitness(i));
    for j=1:I_D
        fprintf(fp_pophist,'%7.3f ',FM_nestpop(i,j));
    end
    fprintf(fp_pophist,' \n');
end

fprintf(fp_pophist,'Best so far: f = %g\n',S_bestval.FVr_oa(1));            
for nn=1:I_D
    fprintf(fp_pophist,'best(%d) = %g\n',nn,FVr_bestnest(nn));
end

fprintf(fp_pophist,' \n');
fclose(fp_pophist);





done = 0;
while (done==0)   
    I_iter = I_iter+1; % iteration counter incremented
    
    % First, do the main bit here: first differential mutation to all of
    % the population, then Lévy flights for worst F_pa 
    
    % 1) Generate new trial solutions by differential mutation...
    FM_new_nestpop=diff_mutation(FVr_bestnest,FM_nestpop,FVr_minbound,FVr_maxbound,F_CR,S_struct); 
    fprintf(1,'\nSuccess rate in differential mutation -- ');    
    [~,~,FM_nestpop,fitness,S_bestval,NFEincr] = get_best_nest(FM_nestpop,FM_new_nestpop,fitness,ofun,S_struct,S_bestval); 
    I_NFE=I_NFE+NFEincr;                               % Update the counter again        

    % 2) Then the Lévy flights:
    % 2.1 update fraction of "abandoned nests" (Lévy flights performed):
    S_struct.F_pa = ( ((I_NFEmax-I_NFE)/I_NFEmax)*F_pa1) + ((I_NFE/I_NFEmax)*F_paMax);   
    % 2.2 Generate new trial solutions by Lévy flights:
    FM_new_nestpop=get_cuckoos(FM_nestpop,FVr_minbound,FVr_maxbound,S_struct);   
    % ...and evaluate the generated set of candidate solutions:
    fprintf(1,'\nSuccess rate in Lévy flights (NP %d; pa %4.2f) -- ', I_NP, S_struct.F_pa);      
    [fnew,FVr_bestmem,FM_nestpop,fitness,S_bestval,NFEincr] = get_best_nest(FM_nestpop,FM_new_nestpop,fitness,ofun,S_struct,S_bestval);     
   	I_NFE=I_NFE+NFEincr;                               % Update the counter   

    if fnew<fmin                                        % Find the best objective so far  
        fmin=fnew;
        FVr_bestnest=FVr_bestmem;
    end    
    
    
    % PROGRESS OUTPUT
    % ---------------
    % best member objective function value and parameter vector are saved in convhist every I_refresh iterations
    % obj. function value and parameter vectors of entire cuckoo population  saved in pophist every 2*I_refresh iters 
    % Screen printing: best member objective function value and parameter vector every I_refresh iterations
    % ----------------------------------------------------------------------------------------------------------------
    
    I_RefreshCycle = I_RefreshCycle+1;    
    
    % in first iterations always output, then only once per I_refresh:
    if ((I_RefreshCycle>I_refresh)|| (I_iter<I_refresh))          
        I_RefreshCycle=0;  % reset refresh counter
        
        %print out the best found objective fuction value first:
        fprintf(1,'\n iter %d, Best:%5.0f,   ',I_iter,S_bestval.FVr_oa(1));    
        format long e;        
        
        % print out best so far decision variable vector
        for n=1:I_D
            fprintf(1,'  %d:%6.4f ',n,FVr_bestmem(n));
        end
            
        % update the matrix for creating a convergence plot:
        I_plotindex = I_plotindex+1;        % data point index
        FM_convplot(1,I_plotindex)=I_NFE;      % N of function evaluations so far
        FM_convplot(2,I_plotindex:I_plotmax)=S_bestval.FVr_oa(1); % best found value so far           
        if (I_plotting==1)                      % if plotting while running ...
            plot(FM_convplot(2,:),FM_convplot(1,:));  % ... then do plot.
        end
        
        % check if also populations are printed:
        if (I_PopPrint>0)            
            fp_pophist = fopen(fp_pophist,'r+');  % open and...
            fseek(fp_pophist, 0, 'eof');       % find end of file to start adding.            
            fprintf(fp_pophist,'Iteration %d\n',I_iter); % iteration first.
            
            % one by one, print solution index, obj.function value, and all
            % decision variable values.
            for i=1:I_NP
                fprintf(fp_pophist,'%3.0f : %8.0f ',i, fitness(i));
                for j=1:I_D
                    fprintf(fp_pophist,'%7.3f ',FM_new_nestpop(i,j));
                end
                fprintf(fp_pophist,' \n');
            end
            
            % at end write also best value & best cand.solution vector.
            fprintf(fp_pophist,'Best so far: f = %g\n',S_bestval.FVr_oa(1));            
            for nn=1:I_D
                fprintf(fp_pophist,'best(%d) = %g\n',nn,FVr_bestmem(nn));
            end            
            fprintf(fp_pophist,' \n');
            fclose(fp_pophist); % ...aand we're done for this generation.
        end
    end
    
    % check if any of the termination criteria has fired yet:
    if (I_iter >= I_itermax) 
        fprintf(1,'Maximum generations reached, finishing.');
        done =1;
    end
    if (I_NFE >= I_NFEmax) 
        fprintf(1,'Maximum NFEs reached, finishing.');
        done =1;
    end
    if (fmin <= F_VTR)
        fprintf(1,'VTR reached, finishing. ');
        done =1;
    end
    %print out some data out at end of each iteration:
    save('LastFinishedCS3iter.mat','FVr_bestnest','S_bestval','I_NFE','FM_convplot');

end

% now that we're done, a few last minor bits:
% add convplot to the best found value struct
S_bestval.convplot=FM_convplot; 

% and print what was done:
disp(strcat('\n\nTotal number of function evaluations = ',num2str(I_NFE))); 
disp(strcat('\nBest found obj.f.value = ',num2str(fmin)));

end


% Generate new cuckoos via Lévy-distribute ramdom walk (Lévy flight)
function FM_nestpop=get_cuckoos(FM_nestpop,FVr_minbound,FVr_maxbound,S_struct)
% perform Levy flights

F_alpha = S_struct.F_alpha;
F_pa = S_struct.F_pa;
I_NP=size(FM_nestpop,1);

% Levy exponent and coefficient; for details see eq. (2.21), pg.16 (Ch 2) of "Nature-Inspired Metaheuristic Algorithms", 2nd Ed. (2010).
% or Saari et al.2022, https://doi.org/10.1016/j.energy.2021.122622 eq(23) 
F_beta=3/2;
F_sigma=(gamma(1+F_beta)*sin(pi*F_beta/2)/(gamma((1+F_beta)/2)*F_beta*2^((F_beta-1)/2)))^(1/F_beta);

%how many eggs are deemed unworthy and discarded?
I_Nfound = max(1,round(F_pa*I_NP));
for j=1:I_Nfound
    FVr_temp=FM_nestpop(j,:);
    % This is a simple way of implementing Levy flights. For standard random walks, use step=1; Levy flights by Mantegna's algorithm
    FVr_q=randn(size(FVr_temp))*F_sigma;
    FVr_p=randn(size(FVr_temp));
    FVr_step=FVr_q./abs(FVr_p).^(1/F_beta);
    FVr_step=FVr_step.*(FVr_maxbound-FVr_minbound);
    
    I_rndindex = ceil((I_NP-1)*rand(1));
    FVr_r1 = FM_nestpop(I_rndindex,:);
  
    % F_alpha is recommended 0.01 with classic CS. In the CS3 the Lévy
    % flights are intended to emphasize global exploration, and several
    % times larger values tend to work better. Not the the vector
    % difference is also to a random vector, not the population best.
    FVr_stepsize=F_alpha*FVr_step.*(FVr_temp-FVr_r1);       
    FVr_temp=FVr_temp+FVr_stepsize.*randn(size(FVr_temp));
    
    % enforce the bounderies: big alpha => big risk to fly off the map:
    if (S_struct.I_bnd_constr==1)
        FM_nestpop(j,:)=simplebounds(FVr_temp,FVr_minbound,FVr_maxbound);
    end
end

end

function [bestF,FVr_bestmem,FM_nestpop,fitness,S_bestval,NFEincr] = get_best_nest(FM_nestpop,FM_newnestpop,fitness,ofun,S_struct,S_bestval) % Find the current best nest
% find the best value so far and sort the population.
% function largely follows the cuckoo search .m code presented in
% Yang, X-S.: Nature-Inspired Optimization Algorithms, 2014.

NFEincr=0;
tries=0;
succ=0;
for j=1:size(FM_nestpop,1)                                      % Evaluating all new solutions
    diff = sum(FM_newnestpop(j,:)-FM_nestpop(j,:));
    if (diff~=0)||(fitness(j)==(2^64))                          % if fitness = 2^64 --> first round, initialize the f's of the initial population  
        S_tempval=feval(ofun,FM_newnestpop(j,:),S_struct);      % fobj(FM_newnestpop(j,:));
        tries=tries+1;                
        NFEincr=NFEincr+1;
        if (S_tempval.FVr_oa(1) <= fitness(j))
            succ = succ+1;
            fitness(j)  =S_tempval.FVr_oa(1);
          
            FM_nestpop(j,:) = FM_newnestpop(j,:);                % if newnestpop has a better member than the nestpop --> update the nestpop with better solution.
            if (S_tempval.FVr_oa(1) <= S_bestval.FVr_oa(1))
                S_bestval = S_tempval;
            end
        end
    end
end
fprintf(1, ' %d/%d \n',succ,tries);                             % print out the success rate
FM_sortpop = horzcat(fitness,FM_nestpop);
FM_sortpop = sortrows(FM_sortpop,-1);                           % FM_sortpop = sortrows(FM_sortpop,1,'descend');
szdim = size(FM_nestpop);                                       % size of the nestpop, szdim(1) = NP, szdim(2) = decision variables
FM_nestpop=FM_sortpop(:,2:(szdim(2)+1));                        % dimension + 1
fitness = FM_sortpop(:,1);

[bestF,~]=min(fitness);                                         % Find the current best
FVr_bestmem=FM_nestpop(szdim(1),:);

end

% perform differential mutation - main search method in CS3.
function FM_new_nestpop=diff_mutation(FVr_bestnest,FM_nestpop,FVr_minbound,FVr_maxbound,F_CR,S_struct) 
% differential mutation almost as original (Storn); minor tweaks such as...
% ... randomize mutation weight factor individually for each each. 
% minimum and maximum bounds can be tweaked for problem at hand.  
    I_strategy = S_struct.I_strategy;
    F_weight = max(0.5,(0.4+0.5*randn(size(FM_nestpop))));
    I_NP=size(FM_nestpop,1);      %FM_bestpop = repmat(FVr_bestnest,I_NP,1);
    FM_bestpop = repmat(FVr_bestnest,I_NP,1);

    FVr_rot  = (0:1:I_NP-1);                    % rotating index array (size I_NP)
    FVr_ind = randperm(4);                      % index pointer array

    FVr_a1  = randperm(I_NP);                   % shuffle locations of vectors
    FVr_rt  = rem(FVr_rot+FVr_ind(1),I_NP);   	% rotate indices by ind(1) positions
    FVr_a2  = FVr_a1(FVr_rt+1);                	% rotate vector locations

    FM_pm1 = FM_nestpop(FVr_a1,:);              % shuffled population 1
    FM_pm2 = FM_nestpop(FVr_a2,:);           	% shuffled population 2

    % then include crossover: K indicates if crossover performed or not for
    % each decision variable value of each candidate solution.
    K=rand(size(FM_nestpop))>F_CR; 
    I_strategy = min(2,max(I_strategy,1));
    if (I_strategy == 1)                        % rand/1 strategy
        stepsize=K.*(F_weight.*(FM_pm1-FM_pm2)); 
    elseif (I_strategy == 2)                    % local-to-best
        stepsize=K.*(F_weight.*(FM_bestpop-FM_nestpop)+F_weight.*(FM_pm1-FM_pm2));
    end
    FM_new_nestpop = FM_nestpop+stepsize;

    if (S_struct.I_bnd_constr==1)
        for j=1:size(FM_new_nestpop,1)
            FVr_temp=FM_new_nestpop(j,:);
            FM_new_nestpop(j,:)=simplebounds(FVr_temp,FVr_minbound,FVr_maxbound); 
        end   
    end
end

% Application of simple constraints
% mirror all violations back into the box
% note that with very long Lévy flights or very large F, sometimes
% mirroring can throw value off the other end. If obj.function cannot
% handle this, it should be checked against and prevented; here very low
% probability of "wasting" a function evaluation is permitted under the
% assumption the objective function evaluation can handle such violations.
function FVr_temp=simplebounds(FVr_temp,FVr_minbound,FVr_maxbound)  
	ns_tmp=FVr_temp;
    I=ns_tmp<FVr_minbound;                              % Apply the lower bound  %ns_tmp(I)=FVr_minbound(I);
    ns_tmp = FVr_temp + 2*I.*(FVr_minbound-FVr_temp);    %ns_tmp(I) = FVr_minbound(I) + (rand(1)*(FVr_maxbound(I)-FVr_minbound(I)));    
    
    J=ns_tmp>FVr_maxbound;                              % Apply the upper bounds     %ns_tmp(J)=FVr_maxbound(J);     %ns_tmp(J) = FVr_minbound(J) + (rand(1)*(FVr_maxbound(J)-FVr_minbound(J)));    
    ns_tmp = FVr_temp - 2*J.*(FVr_temp-FVr_maxbound);    
    FVr_temp=ns_tmp;                                           % Update this new move 
end
