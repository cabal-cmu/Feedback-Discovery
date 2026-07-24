
nscans=round(nsecs/tr);  %round given that TR is not integer     
res=1/200;           % 1/20    raw data resolution in seconds ????????
resfactor=tr/res;
ntpts=nsecs/res;

ts=[];               %declare the variables
poop=[];             %declare the variables   




for subject_number = 1:Nsub       %the cycle runs through the number of subjects 

% cd ~/WOOLIE; 
% homedir=pwd; 
% datadir=sprintf('%s/nma_rsn_data_%d',homedir,testrun);
%datadir=sprintf('thesis_simulations/nma_rsn_data_markov');
datadir=strcat('thesis_simulations/nma_',outputfolder);
[urgh1,urgh2,urgh3]=rmdir(datadir,'s'); 
mkdir(datadir);
subjectname=sprintf('subject%d',subject_number);
mkdir(strcat(datadir,'/',subjectname));







% setup confound evs
confounds=mwdct(round(nsecs/tr),1);  % 14 for 'normal real data', 1 for simulations, using little highpass filtering
                                %this function adjuste the timepoints into
                                %a sum of cosine functions. 1 is only a
                                %straighline
                                
 tmp=devar(confounds(:,2:end),1);   %do not know what is devar (it is not in the file or define in matlab)
 tmp2=[confounds(:,1), tmp];       %adjunct the confound with tmp
 confounds=tmp2;                    %DEFINE CONFOUNDS AS TEMP2
 save(strcat(datadir,'/',subjectname,'/confound_evs.txt'),'confounds','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% irrelevant stuff  JUST CREATE CELLS
N=length(a);  %THE LENGHT OF A IS THE NUMBER OF ROIS. 
P=N; 
roinames=cell(1,N);  %CREATE A CELL WITH THE NAME OF THE ROIS, 1 ROW AND N COLUMNS, EACH CONTAING A ROI NAME: 1 THROUGH N
inputevnames=cell(1,N); %CREATE A CELL FOR THE INPUT EV NAMES. THIS WILL CONTAIN THE INPUTS LATER.
for i=1:N   
    roinames{i}=sprintf('%d',i);    %THE FIRST ELEMENT OF THE CELL IS "1", THEN "2", ETC...STRINGS. 
    inputevnames{i}=sprintf('%d',i); %THE SAME FOR INPUTEVENAMES
end
model=1;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a(1,2) is the input of node 1 into node 2

%HERE EXTRACT THE A MATRIX

a=squeeze(aa(subject_number,:,:));   %EXTRACT FOR EACH SUBJECT MATRIX A. SQUEEZE IS NEEDED DUE TO THE WAY THE MATRIX AA WAS SAVED IN TESTING.M

% b(1,2,3) is stimulation 3 modulating the input of node 1 into node 2
% b=0;
 
 %b=zeros(N,N,1);
b=zeros(N,N,P);  %IN TESTING.M B IS NOT EXPLICITLY DEFINE. HERE WE DEFINE IT FOR THE FORWARD MODEL. IT IS RELATED TO DOMODEXT AS DEFINED IN TESTING.
                    %IN THIS CASE P IS EQUAL TO N. IT IS RELATED TO THE
                    %NUMBER OF INPUTS.


%HERE EXTRACT THE B MATRIX                    
if exist('DOMODext')    %SEARCH IN TESTING.M IF DOMODEXT EXIST...IT IS ONLY DEFINED IN THE CASES OF MODULATION. CHECK THOSE IN THE SWITCH STRUCTURE.
  b=DOMODext;
end


%HERE EXTRACT THE C MATRIX

% c(1,2) is the input of stimulation 2 into node 1
c=squeeze(cc(subject_number,:,:)); %AS WITH AA IT IS NECESSARY TO SQUEEZE IT BECAUSE OF THE WAY IT IS SAVED IN TESTING.M THIS WILL GIVE US THE NXN MATRIX WE WANT



%D MATRIX ???????
% d(1,2,3)
%d=0;
%d=zeros(N,N,1);
d=zeros(N,N,N);    %WHEN MODULATING IS 5 OR LARGER IN TESTING.M
if exist('DOMODint')
  d=DOMODint;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save models


%MAKE A CELL CALL MODELS THAT CONTAINS FOR EACH MODEL AN STRUCTURE AND THE STRUCTURE CONTAINS THE CONNECTIVITY MATRICES A,B,C,D

models{model}.a=a;  
models{model}.b=b;
models{model}.c=c;
models{model}.d=d;

for m=1:length(models)  %WHAT IS THE LENGHT OF MODELS...IF MODELS IS A CELL?????
model_name=strcat('model_art_',num2str(m));  %NAME THE MODEL
mkdir(strcat(datadir,'/',model_name));   %MAKE A DIRECTORY FOR THAT MODEL

save(strcat(datadir,'/',model_name,'/','sigmaa_prior_mean.txt'),'sigma_prior_mean','-ascii'); %THIS IS DEFINED IN TESTING.M
tmp=100;  % originally was 100
save(strcat(datadir,'/',model_name,'/','sigmaa_prior_var.txt'),'tmp','-ascii');  %THIS IS DEFINED IN THE ABOVE LINE.


%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%THE CYLCE RUNS FROM 1 TO N...FOR EACH ROI IT GENERATE DIFFERENT
%BALLOON_PRIORS

for i=1:length(roinames)
  balloon_priors{i}=balloon_defaults();  %FIRST UNLOAD THE DEFAULT VALUES TO USE THEM TO GENERATE NEW VALUES ACCORDING TO THE DISTRIBUTION OF THE PARAMETERS
  %MAKE A CELL WITH STRUCTURES AND SAVE IN EACH STRUCTURE THE PRIORS OF THE BALLON MODEL.
  
  %HERE WE GENERATE NUMBERS FROM A MULTIVARIATE NORMAL DISTRIBUTION
  %MVNRND(MU,SIGMA) WHERE SIGMA IS THE COVARIANCE MATRIX.
  %WE ARE GENERATING DIFFERENT PARAMETERS FOR THE HRF MODEL ACCORDING TO
  %THE DISTRIBUTION OF THE PARAMETERS.
  
  %KAPPA GAMMA
  balloon_priors{i}.balloon_cbf_mean=mvnrnd(balloon_priors{i}.balloon_cbf_mean,balloon_priors{i}.balloon_cbf_cov)';
  %CHOOSE THE MAX BETWEEN THE RESULT OF THE OPERATION ABOVE OR 0.1
  balloon_priors{i}.balloon_cbf_mean=max(balloon_priors{i}.balloon_cbf_mean,0.1);   %SET THE MINIMUM AT 0.1
  
  %HERE AGAIN USE MVNRND USING THE SELECTED CBF_MEAN FROM ABOVE AND THE
  %BALLON_COV MULTIPLIED BY HRF_COVAR
  %TAU ALPHA E0
  balloon_priors{i}.balloon_mean=mvnrnd(balloon_priors{i}.balloon_mean,HRF_covar*balloon_priors{i}.balloon_cov)';  % factor of 4 gives an HRF delay of +- 0.5s
  balloon_priors{i}.balloon_mean=max(balloon_priors{i}.balloon_mean,0.1);
  
  
  %LOGEPSILON
  balloon_priors{i}.balloon2_mean=mvnrnd(balloon_priors{i}.balloon2_mean,balloon_priors{i}.balloon2_cov)';
  %MAX BETWEEN ABOVE COMPUTATION OR 0.1
  balloon_priors{i}.balloon2_mean=max(balloon_priors{i}.balloon2_mean,0.1);
end


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%



%ALL THIS IS TO SAVE THE INFORMATION IN FILES.
for i=1:length(roinames)
tmp=balloon_priors{i}.balloon_cbf_mean;
save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_cbf_mean.txt'),'tmp','-ascii');
tmp=balloon_priors{i}.balloon_cbf_cov;
save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_cbf_cov.txt'),'tmp','-ascii');

tmp=balloon_priors{i}.balloon_mean;
save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_mean.txt'),'tmp','-ascii');
tmp=balloon_priors{i}.balloon_cov;
save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_cov.txt'),'tmp','-ascii');

tmp=balloon_priors{i}.balloon2_mean;
save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon2_mean.txt'),'tmp','-ascii');
tmp=balloon_priors{i}.balloon2_cov;
save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon2_cov.txt'),'tmp','-ascii');
end



%HERE IT SAVE THE MATRICES A,B,C,D IN FILES
tmp=models{m}.a'; %TRASPOSE THE A MATRIX. IT WAS INPUT AS A(1,2)= 1->2, BUT IT IS SAVED, SO 1->2=A(2,1). THIS IS THE CORRECT WAY TO INPUT IT IN THE MODEL.
save(strcat(datadir,'/model_art_',num2str(m),'/matA.txt'),'tmp','-ascii','-double');
tmp=models{m}.c;
save(strcat(datadir,'/model_art_',num2str(m),'/matC.txt'),'tmp','-ascii','-double');



for i=1:N
   %for i=1:size(models{m}.b,1),
    tmp=zeros(N,N);
    %tmp=squeeze(models{m}.b(:,i,:));
    %save(strcat(datadir,'/model_art_',num2str(m),'/matB_into_',roinames{i},'.txt'),'tmp');%,'-ascii','-double');
    dlmwrite(strcat(datadir,'/model_art_',num2str(m),'/matB_into_',roinames{i},'.txt'),tmp,'delimiter','\t');
    %fopen(strcat(datadir,'/model_art_',num2str(m),'/matB_into_',roinames{i},'.txt'),'w');   
end
  for i=1:N
   %for   i=1:size(models{m}.d,1),
     tmp=zeros(N,N);
      %tmp=squeeze(models{m}.d(:,i,:));
    %save(strcat(datadir,'/model_art_',num2str(m),'/matD_into_',roinames{i},'.txt'),'tmp');%,'-ascii','-double');
    dlmwrite(strcat(datadir,'/model_art_',num2str(m),'/matD_into_',roinames{i},'.txt'),tmp,'delimiter','\t');
    %fopen(strcat(datadir,'/model_art_',num2str(m),'/matD_into_',roinames{i},'.txt'),'w');   
  end
end


%%%%%%%%%%%%%%%%%%%%%END SAVING THE FILES WITH THE MODELS: IE. A,B,C,D
%%%%%%%%%%%%%%%%%%%%%MATRICES AND BALLON_PRIORS FOR EACH ROI.


%%%%%%%%THIS IS WHERE THE CREATION OF THE INPUT U, X=AX+CU, STARTS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
% sim stimuli from GMM 

ddd=1/(res*sigma_prior_mean);    % adjust neural stim(?????) strength according to neural time-constant
in_strength=80/sum(exp(-[0:10*ddd]/ddd));  %IN_STRENGHT SET THE MEAN VALUE OF THE AMPLITUDE OF THE STIMULUS SIGNAL.SEE MIX.CENTRES BELOW

%if (N>50)
%  in_strength=in_strength*100;
%end

for e=1:N   %THE CYCLE RUNS FROM 1 TO N ROIS...
    
    %CREATE AN ARRAY CALLED MIX.___

  mix.centres=[0 in_strength];   %MAKE A MATRIX WITH 0 AND IN_STRENGTH. ZERO IS THE BASELINE AND IN_STRENGHT THE MAX VALUE
  
  
  mix.stds=[in_strength/20 in_strength/20];  
  
  trans_prob= res * [0.1 0.4] * inputspeedup;  %INPUTSPEEDUP IS DEFINED IN TESTING.M. IT IS SET TO 1 RIGHT NOW.
  
   if e > Nreal
     trans_prob = res * [0.05 0.03];   %WHEN E IS GREATER THAN NREAL...N ARE NODES MADE FOR THE MODULATION MATRICES.
   end
  
  %for p = trans_prob     %%% check what the associated time constants are for the two bistable states
  %  counter=zeros(10000,1); for i=1:10000, while rand>p, counter(i)=counter(i)+1; end; end;
  %  [pm,pv]=poisstat(1/p);  [ p res * [ mean(counter) pm sqrt(pv) ] ]   % transition probability / empirical mean duration (secs) / theoretical mean & std
  %end
  
 
  
  states=zeros(ntpts,1);   %MAKE A HUGE 1 COLUMN MATRIX WITH ROWS=NTPTS=NSECS*200 FOR STATES.
  
  stim=zeros(ntpts,1);  %ANOTHER EQUAL 1 COLUMN MATRIX FOR STIM. STIM IS STATES PLUS NOISE
  
  
  %THE FIRST VALUE OF STIM IS A RANDOM VALUE OF A NORMAL DISTRIBUTION WITH
  %MEAN (0) AND STANDARD DEVIATION IN_STRENGTH/20
  %THIS IS THE FIRST POINT OF THE TIME SERIES....
  stim(1)=normrnd(mix.centres(1),mix.stds(1));  
  
 
  %THE CYCLE RUNS FROM THE SECOND POINT ALL THROUGH NTPTS...WHICH IS ALSO
  %THE ROW SIZE OF STATES AND STIM MATRICES
  
  
   for t=2:ntpts
       
     if(states(t-1)==0)  %TO ASSIGN VALUE TO EACH POINT CHECK THE VALUE OF THE PREVIOUS POINT.
       
         if(rand<trans_prob(1))
             states(t)=1;  %DEFINE AS 1 THE FIRST STATE IF THE SECOND IS EQUAL TO 1
         else
             states(t)=0;  %ELSE DEFINE AS ZERO
         end
     
     else
         
         if(rand<trans_prob(2))
             states(t)=0;
         else
             states(t)=1;
         end
     end
     
     if(states(t)==0)  %IF THE VALUE OF THE POINT THAT WE ARE EVALUATING IS ZERO
       stim(t)=normrnd(mix.centres(1),mix.stds(1));   %STIM(T) IS A RANDOM VALUE OF THE DISTRIBUTION WITH MEAN(0) AND STDDV IN_STRENGTH/20
     else
       stim(t)=normrnd(mix.centres(2),mix.stds(2));  %ELSE IF THE POINT=1 THEN STIM(T) IS A RANDOM VALUE OF THE DITRIBUTION WITH MEAN(IN_STRENGHT) AND STDDV IN_STRENGTH/20
     end
   end

  
  
%stim= sin( [1:ntpts]/16)';

%    stim=zeros(ntpts,1);
%    stim(2000)=in_strength;
  
%    stim=zeros(ntpts,1);
%    for t=2000:2000:ntpts-2000,
%      for tt=t:t+30
%        stim(tt)=t*in_strength/2000;
%      end
%      for tt=t+1000:t+1500
%        stim(tt)=t*in_strength/2000;
%      end
%    end

%%%%%%%%%%%%%JOE, THIS IS WHERE THE DEFINITION OF THE INPUT U ENDS%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  if inputspeedup<0   %DEFINE STIM MATRIX FOR THE CASES WHEN INPUTSPEED IS LESS THAN ZERO
    stim=3*in_strength*randn(ntpts,1);
  end

  if exist('DONEstim')  % external stimuli already setup in matrix DONEstim   THIS IS DEFINED OUTSIDE. NOT USEFUL IN OUR CASE.
    stim=DONEstim(:,e);
  end
  
  %CHECK THE GRAPHS AT PRUEBA2 TO SEE HOW STATE AND STIM TIMESERIES LOOK
  %ACROSS NTPTS POINTS.
  
  %STIM IS THE U VECTOR AND EACH ROI HAS ITS OWN U (INPUT/STIMULUS) VECTOR
  %IN SMITHS SIMULATIONS.
  
  dlmwrite(strcat(datadir,'/',subjectname,'/',inputevnames{e},'.txt'),stim,'delimiter','\t','precision','%.3f');
  %save(strcat(datadir,'/',subjectname,'/',inputevnames{e},'.txt'),'stim','-ascii','-double');  %CONTAINS THE INPUTS(STIMULUS) U TO EACH ROI.RECALL THAT e GOES FROM 1 TO N ROIS
   %save(strcat(datadir,'/',subjectname,'/',inputevnames{e},'s.txt'),'states','-ascii','-double');
end


%%%%%%%
% sim data
truemodel=1;

%ASSIGN TO A,B,C,D THE VALUES OF A,B,C,D THAT WE PRODUCE IN TESTING.M
a=models{truemodel}.a;
b=models{truemodel}.b;
c=models{truemodel}.c;
d=models{truemodel}.d;
stim_list='1';   %THIS IS JUST A LABEL
roi_list='1';    %THIS IS JUST A LABEL

for i=2:N   %CYCLE RUNS FROM 2 TO N ROIS...WHERE IS THE FIRST ROI????
    stim_list=sprintf('%s,%d',stim_list,i);  %1,2...1,3...1,4...ETC
    roi_list=sprintf('%s,%d',roi_list,i); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%JOE, THIS IS THE CALL TO THE C FUNCTIONS WHERE THE MODEL IS SIMULATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nma_status,nma_result] = system(sprintf('nma_2010_Feb23cut/nma_generate --ld=%s/tmp_model --dd=%s --md=%s/model_art_%i  --tr=%d --resfactor=%d --sub=%s --stim=%s --nn=%s --dm=0 --hm=balloon --ns=%i --sascf --to',datadir,datadir,datadir,truemodel,tr,resfactor,subjectname,stim_list,roi_list,nscans));

nma_status

if nma_status == 0

      for r=1:length(roinames)
        ytrue(:,r)=load(sprintf('%s/tmp_model/signal_%s_bold.txt',datadir,roinames{r})); 
        ztrue_nma_gen(:,r)=load(sprintf('%s/tmp_model/signal_%s_z.txt',datadir,roinames{r}));
      end
      
   
       save(strcat(static,'/BOLDtrue_',num2str(subject_number),'.txt'),'ytrue','-ascii');
       save(strcat(static,'/Ztrue_',num2str(subject_number),'.txt'),'ztrue_nma_gen','-ascii');
      
      
%       save(strcat(datadir,'/',subjectname,'/BOLDtrue.txt'),'ytrue','-ascii');
%       save(strcat(datadir,'/',subjectname,'/Ztrue.txt'),'ztrue_nma_gen','-ascii');

      
      
      % debug HRF variability
      poop=[poop ytrue];

      %%%%%%%% add noise
      for(i=1:size(ytrue,2))   %THE SIZE OF THE SECOND DIMENSION IS THE NUMBER OF ROIS.
        y=ytrue(:,i);  %EXTRACT ALL THE VALUES OF THE i ROI
        y2=y+randn(size(y))*std_BOLD_noise; %ADD NOISE
        y2s(:,i)=y2;  %CREATE Y2S THAT CONTAINS NOISE
        %dlmwrite(strcat(datadir,'/',subjectname,'/',roinames{i},'_data.txt'),y2,'delimiter','\t','precision','%.3f');
        %save(strcat(datadir,'/',subjectname,'/',roinames{i},'_data.txt'),'y2','-ascii');  %MAKE THE FILE WITH Y FOR ROI i
      end
      
      dlmwrite(strcat(static,'/BOLDnoise_',num2str(subject_number),'.txt'),y2s,'delimiter','\t','precision','%.4f');
      %save(strcat(static,'/50_BOLDnoise_',num2str(subject_number),'.txt'),'y2s','-ascii');


      ts=[ts;y2s];  %APPEND YS2...WHICH CONTAINS THE BOLD SIGNALS OF ALL THE ROIS.

else
     nma_result;   %THIS IS AN OUTPUT OF THE ABOVE FUNCTION

end

save(outname);

end

