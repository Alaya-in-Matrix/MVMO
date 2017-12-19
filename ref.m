function mvmo(fhd,iii,jjj,args)
    global proc
    global ps %disc
    global parameter table
    global changed_best
    global stat_shape
    persistent amax
         
    
    proc.finish           = 0;
    proc.i_eval           = 0;
	proc.last_improvement = 1;

    % only the following parameters are required from the user    
    parameter.n_par           = 1;
    parameter.n_tosave        = 25;
    parameter.fs_factor_start = .5;
    parameter.fs_factor_end   = 20;
    parameter.delta_Shape_dyn = 0.25;
    parameter.local_prob      = 0;
    min_eval_LS               = 0.25*proc.n_eval;
    max_eval_LS               = 0.25*proc.n_eval+10;
    parameter.ratio_gute_max  = 0.5;
    parameter.ratio_gute_min  = 0.5;
    parameter.n_random_ini    = round(ps.D/6);
    parameter.n_random_last   = round(ps.D/2);
% end of parameter input

    parameter.scaling=(ps.x_max-ps.x_min);
    proc.opt(2)=1e-10;
      
      indpendent_runs=parameter.n_par*2  ;
 

    %%----------------- Create initial random population ------------------   
    xx=zeros(parameter.n_par,ps.D);
    x_norm=xx;
    for iijj=1:parameter.n_par %Start initialization: x_normalized
        for jjkk=1:ps.D
            xx(iijj,jjkk)=ps.x_min (jjkk) + rand*(ps.x_max(jjkk)-ps.x_min(jjkk));
        end
        x_norm(iijj,:)=(xx(iijj,:)-ps.x_min)./parameter.scaling;
    end % End initialization
    
    x_normalized=x_norm;
    x_normalized_save=x_norm;
    
    %% ------------------ Initialize control parameters --------------------
    n_eval = proc.n_eval;
    n_par = parameter.n_par; 
    D=ps.D ; 
    n_to_save = parameter.n_tosave; 
    fs_factor_start = parameter.fs_factor_start;
    fs_factor_end = parameter.fs_factor_end;
    Shape_dyn = ones(n_par,D);%
    for i=1:n_par
        for k=1:D
          Shape_dyn(i,k)=Shape_dyn(i,k)*40.*rand ;
        end
    end
    delta_Shape_dyn = parameter.delta_Shape_dyn; 
    local_search0_percentage=parameter.local_prob;
    ratio_gute_min=parameter.ratio_gute_min;
    ratio_gute_max=parameter.ratio_gute_max;
    n_randomly_ini = parameter.n_random_ini;
    n_randomly_last = parameter.n_random_last;
    
    %% --------------------- Data structure for the table ---------------------
    table.bests       = NaN*zeros(parameter.n_tosave,ps.D,n_par);
    table.fitness     = Inf*ones(parameter.n_tosave,1,n_par);
    table.objective   = Inf*ones(parameter.n_tosave,1,n_par);
    table.feasibility = zeros(parameter.n_tosave,1,n_par);
    
    %best_history=NaN*zeros(proc.n_eval,ps.D);

    %% ----------------------------- Mapping ----------------------------------
    shape = zeros(n_par,D);
    local_search=zeros(n_par);
    meann = x_normalized;
    meann_app = x_normalized;

   %% ------------------------ Variable selection ----------------------------
    izm = zeros(1,n_par);
    izz = zeros(1,n_par);
    considered = true(n_par,ps.D);
    variance =ones(n_par,ps.D);
    probab=ones(n_par,ps.D);
    values_noneq = zeros(n_to_save,1);
    
   %% ---------------------- Checking correctness ------------------------   
   
    if (n_randomly_last<1)
        n_randomly_last=1;
    end
    
    if (n_randomly_last>ps.D)
        n_randomly_last=ps.D;
    end
    
    if (n_randomly_ini>ps.D)
        n_randomly_ini=ps.D;
    end  
    if (n_eval<=0)
        n_eval=100000d0;
    end

    if (fs_factor_start<=0)
        fs_factor_start=1d0;
    end
    
    if (fs_factor_end<=0)
        fs_factor_end=1d0;
    end

    fs_factor0=fs_factor_start;
    
    yes_n_randomly=true;
    if (n_randomly_ini==n_randomly_last)
        n_randomly=n_randomly_last;
        yes_n_randomly=false;
    end
    
    if (n_to_save<=1)
        n_to_save=2d0;
    end
    
    if (delta_Shape_dyn<=0)
        delta_Shape_dyn=0.2d0;
    end
    
    delta_Shape_dyn0=delta_Shape_dyn;
    delta_Shape_dyn1=1.d0+delta_Shape_dyn0;
    
    yes_fs_factor=true;
    if (fs_factor_start==fs_factor_end)
        yes_fs_factor=false;
    end

    %% ----------------------------- Counters ---------------------------------
    proc.i_eval=0;

    isfc1=0;
    isfc2=0;

    no_in = zeros(1,n_par);
    no_inin = zeros(1,n_par);
    n_par_last=n_par ;
    
    local_search0=local_search0_percentage/100; % Probability of local search (percentage / number of optimization variables)
    goodbad=zeros(n_par,1);
    make_sense=n_eval;
    firsttime=true;
    

%     %% --------------------------------- parameter-SH ---------------------------------
%     if parameter.sn
%     fprintf('===============================================================================\n');
%     fprintf('                   Mean-Variance Mapping Optimization Results                  \n');
%     fprintf('===============================================================================\n');
%     fprintf('Nomenclature:                                                                  \n');
%     fprintf('FE-No. => Objective function evaluation number                                 \n');
%     fprintf('NP     => Number of active particles                                           \n');
%     fprintf('AGBP   => Actual global best particle                                          \n');
%     fprintf('OF     => Global best objective function                                       \n');
%     fprintf('Fit    => Global best fitness                                                  \n');
%     fprintf('Feas   => Feasibilty: 0=No (Constraint violation),  1=Yes                      \n');
% % %     fprintf('No_fmincon  => Calls of local search                                           \n');
% % %     fprintf('Ev_fmincon  => Number of function evaluations performed by local search        \n');
%     fprintf('  FE-No.      NAP        AGBP            OF           Fitness          Feas?    \n');
%     fprintf('  ------     ------     ------        ---------      ---------     ------------ \n');
%     fprintf('\n');
% % %     fprintf('  FE-No.      NAP        AGBP            OF           Fitness       No_fmincon       Ev_fmincon  \n');
% % %     fprintf('  ------     ------     ------        ---------      ---------     ------------      ----------  \n');
%     end
   delta_nrandomly=n_randomly_ini-n_randomly_last;
   A=zeros(n_par,1); 
   local_search_called=0;
   
   while 1       
        %Evaluating the particles.....
        
        for ipp=1:parameter.n_par   
            ff=real(proc.i_eval/n_eval);
            ff2=ff*ff;
           one_minus_ff2=1.d0-ff2*0.99;
           streu=one_minus_ff2*0.5;
                        
            %Determining the relative number of particles belonging to the group of
            %good particles
            border_gute0=ratio_gute_max-ff*(ratio_gute_max-ratio_gute_min);
            border_gute=round(n_par_last*border_gute0);
            if border_gute < 3 || n_par-border_gute < 1
                border_gute=n_par;
            end
            
            
            
            %Selecting the subset of variables to be mutated
            if yes_n_randomly
                n_randomly_X=round(n_randomly_ini-ff*delta_nrandomly);
                n_randomly=round(n_randomly_last+rand*(n_randomly_X-n_randomly_last));
            end
            %Quadratic variation of fs_factor0
            if yes_fs_factor
                fs_factor0=fs_factor_start+ff*(fs_factor_end-fs_factor_start); 
            end
            
            %Denormalize to the original [min, max] range 
            x_normalized(ipp,:) = ps.x_min+parameter.scaling.* x_normalized(ipp,:);
            
            if  local_search(ipp)
                [msgstr, msgid] = lastwarn ;
                TFrcond = strcmp('MATLAB:nearlySingularMatrix',msgid); % Only informative from 'fmincon' function 
                if TFrcond~=0
                    rcond_value0=str2num(msgstr(81:end-1));
                end
                [ffx,oox,ggx,x_normalized(ipp,:),FEVALS] = LocalSearchMVMOSH(fhd,iii,jjj,args,x_normalized(ipp,:),proc.n_eval-proc.i_eval); %Local search
                local_search(ipp)=0;
                local_search_called=local_search_called+1;
                xt.fitness=ffx ; % Constraint handling outsourced so far
                xt.objective=oox;
                
                if xt.fitness==xt.objective
                    xt.feasibility= true; 
                else
                    xt.feasibility= false;
                end

                [msgstr1, msgid1] = lastwarn;
                TFrcond1 = strcmp('MATLAB:nearlySingularMatrix',msgid1);
                if TFrcond1~=0
                    rcond_value=str2num(msgstr1(81:end-1));
                    if (exist('rcond_value0','var')==1) && (rcond_value0 ~= rcond_value) 
                        local_search0=0;
                    end
                end
            else
                [ffx,oox,ggx,x_normalized(ipp,:)]=feval(fhd,iii,jjj,args,x_normalized(ipp,:)); %Problem evaluation
                xt.fitness=ffx; % Constraint handling outsourced so far   
                xt.objective=oox;
                if xt.fitness==xt.objective
                    xt.feasibility= true; %
                else
                    xt.feasibility= false;
                end
            end
            
            if proc.finish
                return;
            end
            
            x_normalized(ipp,:) = (x_normalized(ipp,:)-ps.x_min)./parameter.scaling;
                
            % Store the n-best solution corresponding to the corresponding particle's archive
            Fill_solution_archive();
            meann_app(ipp,:)= meann(ipp,:);
            
            
            if  proc.i_eval  > indpendent_runs  &&  border_gute < n_par; %   n_par  > 1   
%Determining the proportion of good particles
                 if changed_best || firsttime
                    A(1:n_par)=table.fitness(1,:,1:n_par);
                    if firsttime   
                      amax=max(A);
                    end
                    firsttime=false  ;
                    for ia=1:n_par
                      if ~table.feasibility(1,:,ia)
                        A(ia)=A(ia)+amax;
                      end
                    end 
                    [aaa,IX] = sort(A);
                  end
                                        
                    goodbad(IX(border_gute+1:n_par),1)=0;
                    goodbad(IX(1:border_gute),1)=1;
                    iec=randi(border_gute-2,1,1);                                      

                    bestp=IX(1); %First (global best) particle of the group of good particles
                    onep=IX(iec+1); %Randomly selected intermediate particle of the group of good particles
                    worstp=IX(border_gute); %Last particle of the group of good particles
                   
                    %Multi-parent strategy for bad particles
           if ~goodbad(ipp)    
                           shift=streu;
                            beta1 = (2.5)*(rand - shift) ;   
                        for jl=1:D
                               x_normalized(ipp,jl) =table.bests(1,jl,onep)+ beta1*(table.bests(1,jl,bestp)-table.bests(1,jl,worstp));
                               while x_normalized(ipp,jl) > 1.0d0 ||  x_normalized(ipp,jl)<0.0d0    
                                  beta2 = (2.5)*(rand - shift) ;  
                                  x_normalized(ipp,jl) =table.bests(1,jl,onep)+beta2*(table.bests(1,jl,bestp)-table.bests(1,jl,worstp));
                               end
                        end   
                                  meann_app(ipp,1:D) =x_normalized(ipp,1:D); 
                   else
                        x_normalized(ipp,1:D)= table.bests(1,1:D,ipp);   % x_normalized_best(ipp,:); %Local best-based parent assignment for good particles
                    end
            else
                     x_normalized(ipp,:)=table.bests(1,1:D,ipp); %Local best-based parent assignment during independent evaluations
            end
            considered(ipp,1:ps.D) = false;
           
            if rand < local_search0 &&  proc.i_eval < make_sense &&  proc.i_eval > min_eval_LS &&  proc.i_eval < max_eval_LS
              local_search(ipp)=1;
            else
              local_search(ipp)=0;
            end    
            if local_search(ipp)< 1
              VariableSelect1(); % Call random variable selection strategy
            end
            
            %Generation of random input for the mapping function
            for ivar=1:D
                if considered(ipp,ivar) 
                    x_normalized(ipp,ivar) = rand();
                    if shape(ipp,ivar) > 0.d0 
                      sss1 = shape(ipp,ivar);
                      sss2 = sss1;
                      delta_ddd_x=delta_Shape_dyn0*(rand()-0.5d0)*2.0d0+delta_Shape_dyn1; 
                      if (shape(ipp,ivar) > Shape_dyn(ipp,ivar))
                        Shape_dyn(ipp,ivar) = Shape_dyn(ipp,ivar)*delta_ddd_x;
                      else
                        Shape_dyn(ipp,ivar) = Shape_dyn(ipp,ivar)/delta_ddd_x;
                      end
                        if rand < 0.5
                       sss1=Shape_dyn(ipp,ivar) ;
                        else
                       sss2=Shape_dyn(ipp,ivar);  
                        end
                     
             
                     if ipp==1 & ivar==1
                        isfc1=isfc1+1;
                        stat_shape.iter(:,isfc1)=proc.i_eval;
                        stat_shape.sfaco1(:,isfc1)=sss1;
                        stat_shape.sfaco2(:,isfc1)=sss2;
                        stat_shape.shapetr(:,isfc1)=shape(ipp,ivar);
                        stat_shape.Shape_dyntr(:,isfc1)=Shape_dyn(ipp,ivar);
                     end
                     
                       if rand >0.5   
                         fs_factor=fs_factor0*(1.d0+rand);
                       else
                         fs_factor=1.d0+fs_factor0*(1.d0-rand)*0.25; 
                       end
                      sss1=sss1*fs_factor;
                      sss2=sss2*fs_factor;        
                      x_normalized(ipp,ivar)=...
                            h_function(meann_app(ipp,ivar),sss1,sss2,x_normalized(ipp,ivar)); %Mapping function
                      if x_normalized(ipp,ivar)>0.98 && rand < 0.2
                          x_normalized(ipp,ivar)=1.0;
                      elseif x_normalized(ipp,ivar)<0.02 && rand < 0.2
                          x_normalized(ipp,ivar)=0.0;
                      end
                      
                        
                    if ipp==1 & ivar==1
                        isfc2=isfc2+1;
                        stat_shape.sfac1(:,isfc2)=sss1;
                        stat_shape.sfac2(:,isfc2)=sss2;
                        stat_shape.xrmean(:,isfc2)=meann_app(ipp,ivar);
                        stat_shape.xrnorm(:,isfc2)=x_normalized(ipp,ivar);
                    end
                        
                        
                    end
                end
            end
        end %End n_par loop
    end %End while loop
    
        %% ----------------------- Complementary functions ------------------------

     function [ffx,oox,ggx,xn_out,FEVALS] = LocalSearchMVMOSH(testcase,iii,jjj,args,xx_yy,FEsAllowed)
        global PPL GGL
        if FEsAllowed <= 0, return, end
        lb=ps.x_min;
        ub=ps.x_max;
        Aeq=[]; %AeqX=Beq
        Beq=[];
        AA=[]; %AX<=B
        BB=[];
        
        
        options=optimset('Display','off','algorithm','interior-point','UseParallel','never','MaxFunEvals',FEsAllowed) ; %,'MaxFunEvals',500);  %,'TolFun',1e-15
        [Xsqp, FUN , ~ , output]=...
            fmincon(@(xx_yy)LSearch(xx_yy,testcase,iii,jjj,args),xx_yy,AA,BB,Aeq,Beq,lb,ub,[],options);
        
        FEVALS=output.funcCount ;
        for nvar=1:size(xx_yy,2)
            if isnan(Xsqp(1,nvar))
                Xsqp=xx_yy;
                break;
            end
        end
        
        xn_out = Xsqp;
        ffx=FUN;
        oox=PPL; 
        ggx=GGL;
    end

    function J=LSearch(xx_yy2,testcase,iii,jjj,args)
        global PPL GGL 
        [J,PPL,GGL,~] = feval(testcase,iii,jjj,args,xx_yy2);
        
    end
        
        
        
        function Fill_solution_archive()
        
  
        no_in(ipp) = no_in(ipp)+1;
        changed = false;
        changed_best=false;

        if no_in(ipp) ==1 % the solution coming to the table for the first time large
            table.fitness(1:n_to_save,:,ipp) = 1.e200;
            table.feasibility(1:n_to_save,:,ipp) = 0;
            table.bests(1,:,ipp)   = x_normalized(ipp,:) ;
            table.fitness(1,:,ipp) = xt.fitness ;
            table.objective(1,:,ipp) = xt.objective;
            table.feasibility(1,:,ipp) = xt.feasibility    ;  %repmat(,n_to_save,1);

            no_inin(ipp)=no_inin(ipp)+1;
            changed_best=true;
            
        else % not for the first time and check for the update
           i_position =0;
               for ij=1:n_to_save 
                   if (xt.fitness < table.fitness(ij,:,ipp) && xt.feasibility == table.feasibility(ij)) || (table.feasibility(ij,:,ipp) <  xt.feasibility)                                 
                       i_position = ij;
                       changed =true;
                       if (ij<n_to_save)
                           no_inin(ipp) = no_inin(ipp)+1; % how many times good solutions were found   
                       end
                       break;
                   end
               end
        end

        if changed   % if the new individual is better than any archived individual.
                     % Move the individuals and corresponding fitness values
                     % downward so that the individuals are sorted based on the
                     % fitness value in a descending order             
            nnnnnn = n_to_save;
            if i_position==1
              changed_best=true;
            end     
            if (no_inin(ipp) < n_to_save); nnnnnn = no_inin(ipp); end
            isdx=nnnnnn:-1:i_position+1;
            table.bests(isdx,:,ipp) = table.bests(isdx-1,:,ipp);
            table.fitness(isdx,:,ipp)= table.fitness(isdx-1,:,ipp);
            table.objective(isdx,:,ipp)= table.objective(isdx-1,:,ipp);
            table.feasibility(isdx,:,ipp)= table.feasibility(isdx-1,:,ipp);
 
            % save the new best
            table.bests(i_position,:,ipp) = x_normalized(ipp,:);
            table.fitness(i_position,:,ipp) = xt.fitness;
            table.objective(i_position,:,ipp) = xt.objective;
            table.feasibility(i_position,:,ipp) = xt.feasibility;

            % calculation of mean and variance
            if ((no_inin(ipp)>=2))
                for ivvar = 1:D
                    [meann(ipp,ivvar),variance(ipp,ivvar)] = mv_noneq(nnnnnn,table.bests(1:nnnnnn,ivvar,ipp));
                end
                id_nonzero = (variance(ipp,:) > 1.1d-100);
                shape(ipp,id_nonzero)=-log(variance(ipp,id_nonzero));
            end
        end
        end
        
    

    
    
        function VariableSelect1()
          mode=4 ;
          if rand > 1
              mode=4;
          end
          n_var=ps.D;
          switch mode
            case 1
                for ii=1:n_randomly
                    isrepeat = false;
                    while ~isrepeat
                        inn=round(rand*(n_var-1))+1;
                        if (~considered(ipp,inn))
                            isrepeat = true;
                        end
                    end
                    considered(ipp,inn)=true;
                end
            case 2
                in_randomly=0;
                isrepeat = false;
                izz(ipp)=round(rand*(n_var-1))+1; %NEWWWWW
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                        izz(ipp)=n_var;
                    end
                    considered(ipp,izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly)) 
                        isrepeat = true;
                    end
                end
            case 3
                in_randomly=0;
                izm(ipp)=izm(ipp)-1;
                
                             izm(ipp)=round(rand*(n_var-1))+1;
               
                
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                izz(ipp)=izm(ipp);
                isrepeat = false;
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                         izz(ipp)=n_var;
                    end
                    considered(ipp,izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly)) 
                        isrepeat = true;
                    end
                end   
            case 4
                izm(ipp)=izm(ipp)-1;
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                considered(ipp,izm(ipp))=true;
                if (n_randomly>1)  
                    for ii=1:n_randomly-1
                        isrepeat = false;
                        while ~isrepeat
                            inn=round(rand*(n_var-1))+1;
                            if (~considered(ipp,inn))
                                isrepeat = true;
                            end
                        end
                        considered(ipp,inn)=true;
                    end
                end
            case 5
                  summep=sum(probab(ipp,:));
                  wahr=probab(ipp,:)/summep;
                  SS0=0.d0;
                  SS=zeros(1,n_var);
                  for imm=1:(n_var-1)
                      SS0=SS0+wahr(imm);
                      SS(imm+1)=SS0;
                  end
                  for ijlr=2:n_var
                      wahr(ijlr)=wahr(ijlr)+SS(ijlr);
                  end 
                  for ltt=1:n_randomly
                       isrepeat = false;
                       while  ~isrepeat
                        rnn=rand;
                        for irw=1:n_var
                          if considered(ipp,irw)
                           continue
                          end
                          if (irw==1)
                              unten=0.d0;
                          else
                              unten=wahr(irw-1);
                          end
                         if (rnn>=unten)&&(rnn<wahr(irw))
                              isrepeat = true;
                             considered(ipp,irw)=true;
                             break;
                          end
                        end
                       end
                  end
          end
            
        end
    
    
    
    
    
    
         function [vmean,vvariance] = mv_noneq(n_to_save,values)
            iz =1; 
            values_noneq(iz)=values(1);
            for ii_jj=2:n_to_save
                izz = iz;
                gleich = false;
                for kk_ii=1:izz
                    if  abs(values_noneq(kk_ii) - values(ii_jj)) <  1.d-70 ; 
                        gleich = true;
                        break;
                    end
                end
                if ~gleich;
                    iz = iz+1;
                    values_noneq(iz)=values(ii_jj);
                end
            end
             vmean = values_noneq(1);
            if (iz>1)
               for kk_ii=2:iz
                    vmean = vmean+values_noneq(kk_ii);
               end
               
               vmean = vmean/iz;
           end
           vvariance = 0.d0;
           if (iz>1)
                for kk_ii=1:iz
                    vvariance =vvariance+(values_noneq(kk_ii)-vmean)*(values_noneq(kk_ii)-vmean);
                end
                  vvariance = vvariance/iz;
           else
                  vvariance=1.0d-100 ;
            end
        end 
    end

    %% Evacuated h-function
    function x = h_function(x_bar,s1,s2,x_p)
        H = x_bar .* (1.d0 - exp(-x_p.*s1)) + ...
            (1.0d0 - x_bar) .* exp(-(1.d0-x_p).*s2);              
        H0 = (1.d0 - x_bar) .* exp(-s2);
        H1 = x_bar .* exp(-s1);
        x = H + H1 .* x_p + H0 .* (x_p - 1.d0);
    end
    
    

