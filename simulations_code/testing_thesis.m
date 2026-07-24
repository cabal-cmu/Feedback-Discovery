function []=testing_thesis(tr,Nsub,Amatrix,outputfolder)

%ADD PATHS TO VARIOUS EFFECTIVE CONNECTIVITY ANALYSES

% % % addpath('~/matlab/wtc-r16'); addpath('~/matlab/funcconn'); addpath('~/matlab/lingam-1.4.2/code')
% % % addpath('~/matlab/FastICA_25'); addpath('~/matlab/icasso122'); addpath('~/matlab/L1precision');
% % % addpath('~/matlab/GCCA_toolbox_sep21'); addpath('~/matlab/GCCA_toolbox_sep21/utilities');
% % % addpath('~/WOOLIE'); cd ~/matlab/biosig4octmat-2/biosig ; install ; cd ~/WOOLIE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std_BOLD_noise=1;   % added thermal (measurement) noise
%tr=1.20;               % temporal sampling of BOLD timeseries
nsecs=10*60;          % session length (per "subject") in seconds
%Nsub=1;            % number of "subjects"; must be even: needs to be at least as large as max(N*(N-1)/2,50) to get enough samples in null
HRF_covar=4;        % variability in HRF delay; 4 gives ~0.5s variability
Nsubnets=1;         % how many 5-node subnets do we setup?
sigma_prior_mean=20;    % 1/neural lag;  e.g. 20 or 10, originally was 1 in DCM/NMA - but that gives a 1s lag neurally!!
                        %THIS IS ALSO USED IN THE BOLD SIGNAL GENERATION
                        %FUNCTION

%INITIAL DEFINITION OF PARAMETERS FOR SIMULATIONS. SOME ARE CHANGED IN THE
%CASES OF TESTRUN
SharedInputs=0; %INPUTS MAY NOT BE EXCLUSIVE TO ITS CORRESPONDING ROI
GroupTwoDiff=0; 
ConfoundGlobalMean=0; 
BackwardsConnections=0; 
OneBackwardsConnections=0;
TwoTogetherBackwardsConnections=0;
TwoSeparateBackwardsConnections=0;
OneBackwardsConnectionsWithinCycle=0;
TwoTogetherBackwardsConnectionsWithinCycle=0;
TwoSeparateBackwardsConnectionsWithinCycle=0;
tennodes15edges=0;
tennodes19edges=0;
tennodes15edges2=0;
tennodes19edges2=0;
tennodes15edges3=0;
tennodes19edges3=0;
tennodes15edges4=0;
tennodes19edges4=0;

Cyclic=0; 
OneExternal=0; 
StrongLinks=0; 
BadROIs=0; 
MoreLinks=0; 
NegativeConnections=0; 
Modulation=0; 
ThreeNodes=0; 
NormICOV=0; 
inputspeedup=1; 
external_A=0;

simdesc='';  %THESE ARE ONLY LABELS OF THE SIMULATIONS

% 
% 
%switch testrun
 
%    case {1}
%     Nsubnets=2;  %NUMBER OF SUBNETS OF 5 NODES EACH
%     amatrix=1;    %all the A coeff positive
%     
%     case{2}
%         Nsubnets=2;
%         amatrix=2; %all the A coeff negative
%         
%     case{3}
%         Nsubnets=2;
%         amatrix=3; %all the A coeff random
%     
%   case {1}
%   case {2}
%     Nsubnets=2;
%   case {3}
%    Nsubnets=3;
%   case {4}
%     Nsubnets=10;
%   case {5}
%     nsecs=3600;  %NUMBER OF SECONDS
%   case {6}
%     nsecs=3600; Nsubnets=2;
%   case {7}
%     nsecs=15000;
%   case {8}
%     SharedInputs=1; simdesc='shared inputs';
%   case {9}
%     nsecs=15000; SharedInputs=1; simdesc='shared inputs';
%   case {10}
%     ConfoundGlobalMean=0.5; simdesc='global mean confound';
%   
%   case {11}
%     BadROIs=1; Nsubnets=2; simdesc='bad ROIs (timeseries mixed with each other)';
%   case {12} % BadROIs - not mix with each other (like above), but mix with other random shit
%     BadROIs=2; Nsubnets=2; simdesc='bad ROIs (new random timeseries mixed in)';
%   
%     
%  case {13}
%     BackwardsConnections=1; simdesc='backwards connections';
%     
%     %%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%TWO CYCLES
%     
%  case {131}      
%     OneBackwardsConnections=1; %implement just one backwards connection
%  case {132}    
%     TwoTogetherBackwardsConnections=1; %implement two backward connections together
%  case {133}
%     TwoSeparateBackwardsConnections=1; %implement two backwards connections separated
%  
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%
%  case{134}
%     tennodes15edges=1;
%     Nsubnets=2;
%     
%  case{135}   
%     tennodes19edges=1;
%     Nsubnets=2;
%     
%  case{136}
%     tennodes15edges2=1;
%     Nsubnets=2;
%  
%  case{137}
%     tennodes19edges2=1;
%     Nsubnets=2;
% 
%     %%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%TWO CYCLES WITHIN CYCLES
%  case {138}      
%     OneBackwardsConnectionsWithinCycle=1; %implement just one backwards connection
%  case {139}    
%     TwoTogetherBackwardsConnectionsWithinCycle=1; %implement two backward connections together
%  case {140}
%     TwoSeparateBackwardsConnectionsWithinCycle=1; %implement two backwards connections separated
%     
%  case{141}
%     tennodes15edges3=1;
%     Nsubnets=2;
%     
%  case{142}   
%     tennodes19edges3=1;
%     Nsubnets=2;
%     
%  case{143}
%     tennodes15edges4=1;
%     Nsubnets=2;
%  
%  case{144}
%     tennodes19edges4=1;
%     Nsubnets=2;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%55   
%   case {14}
%     Cyclic=1; simdesc='cyclic connections';
%   case {15}
%     std_BOLD_noise=0.1; StrongLinks=1; simdesc='stronger connections';
%   case {16}
%     MoreLinks=1; simdesc='more connections';
%   case {17}
%     Nsubnets=2;  std_BOLD_noise=0.1;
%   case {18}
%     HRF_covar=0;
%   case {19}
%     tr=0.25; std_BOLD_noise=0.1; sigma_prior_mean=10; simdesc='neural lag=100ms';
%   case {20}
%     tr=0.25; std_BOLD_noise=0.1; sigma_prior_mean=10; HRF_covar=0; simdesc='neural lag=100ms';
% %   case {21}
% %     GroupTwoDiff=1; simdesc='2-group test';
%   case {22}
%     std_BOLD_noise=0.1; Modulation=1; simdesc='nonstationary connection strengths';
%   case {23}
%     std_BOLD_noise=0.1; Modulation=2; simdesc='stationary connection strengths';
%   case {24}
%     std_BOLD_noise=0.1; StrongLinks=1; OneExternal=1; simdesc='only one strong external input';
%   case {25}
%     nsecs=300;
%   case {26}
%     nsecs=150;
%   case {27}
%     nsecs=150; std_BOLD_noise=0.1;
%   case {28}
%     nsecs=300; std_BOLD_noise=0.1;
% % 
% %   case {34}
% %     inputspeedup=3;
% %   case {35}
% %     Nsubnets=2; inputspeedup=3;
% %   case {36}
% %     tr=0.25; std_BOLD_noise=0.1; sigma_prior_mean=10; simdesc='neural lag=100ms'; inputspeedup=3;
% %   case {40}
% %     Nsubnets=2; NegativeConnections=1; simdesc='one negative connection';
% %   case {50}   % for Reza's triplets work
% %     Nsubnets=2; std_BOLD_noise=0.1; Modulation=10; simdesc='modulatory connections'; Nsub=100;
% % 
% %   case {60}
% %     Nsubnets=10; nsecs=3600; tr=0.8; Nsub=50; simdesc='low-TR HCP';
% %   case {61}
% %     Nsubnets=1;  nsecs=3600; tr=0.8; Nsub=50; Modulation=1; simdesc='low-TR HCP with modulations A';
% %   case {62}
% %     Nsubnets=1;  nsecs=3600; tr=0.8; Nsub=50; Modulation=2; simdesc='low-TR HCP without modulations A';
% %   case {63}
% %     Nsubnets=10; nsecs=3600; tr=0.8; Nsub=50; Modulation=3; simdesc='low-TR HCP with modulations B';
% %   case {64}
% %     Nsubnets=10; nsecs=3600; tr=0.8; Nsub=50; Modulation=4; simdesc='low-TR HCP without modulations B';
% % 
% %   case {70}
% %     Nsubnets=40; nsecs=3600; tr=0.72; std_BOLD_noise=0.1;  simdesc='low-TR HCP';
% %      average ROI is ~1000 2mm voxels; stddev factor = /33
% %      typical HCP stddev = 2%;  2/33=0.06
% % 
% %   case {101}
% %     ThreeNodes=1;
% %   case {102}
% %     ThreeNodes=2;
% %   case {103}
% %     ThreeNodes=3;
% %   case {104}
% %     ThreeNodes=4;
% %   case {105}
% %     ThreeNodes=5;
% %   case {106}
% %     ThreeNodes=6;
% %   case {107}
% %     ThreeNodes=7;
% %   case {108}
% %     ThreeNodes=8;
% %   case {109}
% %     ThreeNodes=9;
% %   case {110}
% %     ThreeNodes=10;
% %   case {111}
% %     ThreeNodes=11;
% %   case {112}
% %     ThreeNodes=12;
% % 
% %   case {300}
% %     external_A=1;  %define an A matrix specified outside the code
% % 

%%%%% CASES FOR EXPERIMENTS: NEW PAPER SEPTEMBER
%%%%% 2013%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
  %  case{1}  % 5 NODES AND 2-CYCLE IN EACH CONNECTION. 10 SUBJECTS.
  %      Nsubnets=1;
  %      Nsub=10;


 %end
% 
 


outname=sprintf('netsim-thesis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now setup the underlying network topology, etc.


for n = 1:Nsub 
%    for n = 1:Nsub %this n is again refered at the end of the cycle. Nsub/2 is the half of the subjects (or the half of the simulations. number of subjects are interpreted as number of runs of the cycle.)

  clear a c   %CLEAR THE MATRIX A AND C

% if external_A>0   % pre-specified A from file
%   grot=load(sprintf('A_%d',testrun));  grot=grot.A>0;
%   N=size(grot,1); Nreal=N; Nsubnets=N/5;
%   a=max(min(randn(N)/20+0.5,0.6),0.4) .* grot;   % tighter range than originally used
%   a=a-eye(N);
%   c=eye(N);
% elseif Nsubnets>39    % huge Netsims for HCP
%   N=Nsubnets*5; Nreal=N;
%   a=-eye(N);
%   column_means=.3+rand(1,N)*.6; % each target has between 30%-90% of potential sources connected (but 30% will be removed....correct for that?)
%   for i=1:N      % sources
%     for j=i+1:N  % targets
%       if rand<column_means(j)
%         strength= 0.01+min(0.7,.1*exp(.5*randn));   % weights on entries: range 0.2:0.9 (with decreasing probability)
%         a(i,j)=strength;  a(j,i)=strength;
%         grot=rand;
%         if grot<0.3        % 10-30% of connections should be unidirectional
%           if rand<0.5 a(i,j)=0; else a(j,i)=0; end
%         elseif grot<0.7    % 40% of connections should be asymmetric in strength
%           if rand<0.5
%             a(i,j)=a(i,j)* 0.01*10.^(rand);   % multiply by a fraction from 0.01 to 0.1
%           else
%             a(j,i)=a(j,i)* 0.01*10.^(rand);
%           end 
%         end
%         if rand<0.9    % make 20% connections negative
%           if rand<0.5 a(i,j)=-a(i,j); else a(j,i)=-a(j,i); end
%         end
%       end
%     end
%   end
%   c=diag(rand(1,N));
% else
%define c matrix by making the diagonal with random numbers

 %this else is in case Nsubnets is less tha 39
    
    %Here we set the connection matrix A for Subnets of 5 nodes, where
    %1-2,2-3,3-4,4-5...no cycle and each node connects into itself with a
    %decay, ie. negative parameter.
    
    
%   if Nsubnets==1000;
%     
%   for I = 1:5:5*Nsubnets   %cycle goes from 1 to N nodes in steps from 5, eg. 1 6 11 ... This is to separate the N into subnets of 5.  
%     for i=I:I+4, a(i,i)=-1; end;    
% 
%     for i=I:I+3                     
%       a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%     end
%       a(I,I+4)=max(min(randn/10+0.4,0.6),0.2);
%       
%   end
    
% 
% %Build private-sets
% p=1:5:Nsubnets*5;
% for j=1:5:Nsubnets
%     for k=j:j+3
%         a(p(j),p(k+1))=max(min(randn/10+0.4,0.6),0.2);
%     end
% end
%     
% 
% %Build sargeant-sets
% s=1:25:Nsubnets*5;
% for j=1:5:Nsubnets/5;
%     for k=j:j+3
%         a(s(j),s(k+1))=max(min(randn/10+0.4,0.6),0.2);
%     end
% end;
%         
% 
% %Build lieutenant-sets
% l=1:125:Nsubnets*5;
% for j=1:5:Nsubnets/25;
%      for k=j:j+3
%         a(l(j),l(k+1))=max(min(randn/10+0.4,0.6),0.2);
%      end
% end
% 
% %Build colonel-ring 
% co=1:625:Nsubnets*5;
% for j=1:(Nsubnets/125)-1;
%      a(co(j),co(j+1))=max(min(randn/10+0.4,0.6),0.2);
% end
%      a(co(1),co(end))=max(min(randn/10+0.4,0.6),0.2);
%      
     
% %%%%% all 2-cylce

%   if Nsubnets==1
%         
%        for I = 1:5:5*Nsubnets
% 
%     for i=I:I+4, a(i,i)=-1; end;    
% 
%     for i=I:I+3                     
%       a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%       a(i+1,i)=max(min(randn/10+0.4,0.6),0.2);
%     end
%       a(I,I+4)=max(min(randn/10+0.4,0.6),0.2);
%       a(I+4,I)=max(min(randn/10+0.4,0.6),0.2);
%     
%        end







%      if Nsubnets==1
%         
%        for I = 1:5:5*Nsubnets
% 
%     for i=I:I+4, a(i,i)=-1; end;    
% 
%     for i=I:I+3                     
%       a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%     end
%       a(I,I+4)=max(min(randn/10+0.4,0.6),0.2);
%     
%        end
       
%
% %Build private-sets
% p=1:5:Nsubnets*5;
% for j=1:5:Nsubnets
%     for k=j:j+3
%         a(p(j),p(k+1))=max(min(randn/10+0.4,0.6),0.2);
%     end
% end
%     
% 
% %Build sargeant-sets
% s=1:25:Nsubnets*5;
% for j=1:5:Nsubnets/5;
%     for k=j:j+3
%         a(s(j),s(k+1))=max(min(randn/10+0.4,0.6),0.2);
%     end
% end;
%         
% 
% 
% %Build colonel-ring 
% co=1:125:Nsubnets*5;
% for j=1:(Nsubnets/25)-1;
%      a(co(j),co(j+1))=max(min(randn/10+0.4,0.6),0.2);
% end
%      a(co(1),co(end))=max(min(randn/10+0.4,0.6),0.2);


    
        
        
        
        
        
        
      
      
      
      
      
      
      
      
     %for i=I:I+4 
      %  a(i,i)=-1; 
     %end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.
     
    %  for i=I:I+3                     % 1st forward connection
     % a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
    %end

    
     
%      if amatrix==1; %ALL THE A MATRIX COEFFICIENTS (EXCEPT THE DIAGONAL) ARE POSITIVE
%          
%         for i=I:I+3                     % 1st forward connection.  the for-cycle goes from ROI1 to ROI4, from ROI6 to ROI9. Each subnet at the time.
%       a(i,i+1)=max(min(randn/10+0.4,0.6),0.2); %a matrix fro 1-2;2-3;3-4;4-5//6-7;7-8;8-9;9-10;etc....HERE we are setting the connection in each subnet.
%         end
%         a(3,8)=max(min(randn/10+0.4,0.6),0.2);
%         a(I,I+4)=max(min(randn/10+0.4,0.6),0.2);
%         
%      elseif amatrix==2;
%          for i=I:I+3
%        a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%        a(i,i+1)=-a(i,i+1);
%          end
%        a(3,8)=-max(min(randn/10+0.4,0.6),0.2);
%        a(I,I+4)=-max(min(randn/10+0.4,0.6),0.2);
%      
%          
%      else
%            for i=I:I+3
%                r=binornd(1,0.5);
%                     if r==1;
%                         a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%                     else
%                         a(i,i+1)=-max(min(randn/10+0.4,0.6),0.2);   
%                     end
%            end
%            r=binornd(1,0.5);
%            if r==1;
%               a(3,8)=max(min(randn/10+0.4,0.6),0.2);
%               a(I,I+4)=max(min(randn/10+0.4,0.6),0.2);
%            else
%               a(3,8)=-max(min(randn/10+0.4,0.6),0.2);
%               a(I,I+4)=-max(min(randn/10+0.4,0.6),0.2);
%            end
%             
%          
%            
%      end
%                
           
             
    
    
%     
%         if MoreLinks == 1  %look for case(16) 
%         a(I+1,I+3)=max(min(randn/10+0.4,0.6),0.2);   %this set an edges between ROI2-ROI4; ROI7-ROI9; etc.
%         a(I+2,I+4)=max(min(randn/10+0.4,0.6),0.2);      %ROI3-ROI5; ROI8-ROI10, etc.
%         end
% 
% %        for i=I:I+4, for j=I:I+4        % other forward connections, 20% chance of being negative
% %            if j~=i && a(j,i)==0 && a(i,j)==0 && rand>0.75
% %              a(i,j)= sign(rand-0.2) * max(min(randn/10+0.4,0.6),0.2);
% %            end
% %        end; end
% 
%    
%     if BackwardsConnections == 1        % backward connections
%       for i=I:I+4     %each subnet at a time 
%           for j=I:I+4
%             if a(j,i)>0 && rand>=0.5   %rand is to include a random quality to the simulations.
%                a(i,j)= -max(min(randn/10+0.4,0.6),0.2);  %make negative the opposite of a(i,j), ie. a(j,i) will be negative.
%             end
%            end; 
%       end
%     end
%     
%     
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CYCLES SECTION
%                                                   
%      if OneBackwardsConnections==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.  
%          for i=I:I+3                     % 1st forward connection
%            a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%          end
%            a(2,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1   
%      end
%                                                      
%      if TwoTogetherBackwardsConnections==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.
%          for i=I:I+3                     % 1st forward connection
%            a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%          end         
%           a(2,1)= max(min(randn/10+0.4,0.6),0.2); %make an edge going from ROI2->ROI1
%           a(3,2)= max(min(randn/10+0.4,0.6),0.2); %make an edge going from ROI3->ROI2        
%      end
%      
%      
%       if TwoSeparateBackwardsConnections==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.    
%          for i=I:I+3                     % 1st forward connection
%            a(i,i+1)=max(min(randn/10+0.4,0.6),0.2);
%          end
%            
%           a(2,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%           a(4,3)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI4->ROI3        
%      end
% 
%     if OneBackwardsConnectionsWithinCycle==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.  
%          
%          a(1,2)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(2,3)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(3,4)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(4,5)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(5,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          
%          a(2,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%       end
%           
%       if TwoSeparateBackwardsConnectionsWithinCycle==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.    
%          a(1,2)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(2,3)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(3,4)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(4,5)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(5,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%            
%          a(2,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(4,3)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI4->ROI3        
%      end
% 
%                                                      
%      if TwoTogetherBackwardsConnectionsWithinCycle==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;    % no self-connections(?)  set the diagonal to -1. This negative number models a decay.
%          a(1,2)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(2,3)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(3,4)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(4,5)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
%          a(5,1)= max(min(randn/10+0.4,0.6),0.2);  %make an edge going from ROI2->ROI1
% 
%          a(2,1)= max(min(randn/10+0.4,0.6),0.2); %make an edge going from ROI2->ROI1
%          a(3,2)= max(min(randn/10+0.4,0.6),0.2); %make an edge going from ROI3->ROI2        
%      end
% 
%      
%      if  tennodes15edges==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,4)=max(min(randn/10+0.4,0.6),0.2); 
%          a(1,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(1,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%       if  tennodes19edges==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,4)=max(min(randn/10+0.4,0.6),0.2); 
%          a(1,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(1,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(2    ,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          %two cycles
%          a(3,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(10,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,9)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%       if  tennodes15edges2==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,3)=max(min(randn/10+0.4,0.6),0.2); 
%          a(10,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%       if  tennodes19edges2==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,3)=max(min(randn/10+0.4,0.6),0.2); 
%          a(10,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          
%          %two cycles
%          a(8,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(10,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,7)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%          
%      if  tennodes15edges3==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,4)=max(min(randn/10+0.4,0.6),0.2); 
%          a(1,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%       if  tennodes19edges3==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,4)=max(min(randn/10+0.4,0.6),0.2); 
%          a(1,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          %two cycles
%          a(2,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(10,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,9)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%       if  tennodes15edges4==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,3)=max(min(randn/10+0.4,0.6),0.2); 
%          a(10,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%       if  tennodes19edges4==1
%          for i=I:I+4 
%            a(i,i)=-1; 
%          end;
%          a(1,3)=max(min(randn/10+0.4,0.6),0.2); 
%          a(10,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(2,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(3,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,10)=max(min(randn/10+0.4,0.6),0.2);
%          a(4,2)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,6)=max(min(randn/10+0.4,0.6),0.2);
%          a(6,9)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,1)=max(min(randn/10+0.4,0.6),0.2);
%          a(5,7)=max(min(randn/10+0.4,0.6),0.2);
%          a(8,5)=max(min(randn/10+0.4,0.6),0.2);
%          a(9,8)=max(min(randn/10+0.4,0.6),0.2);
%          
%          %two cycles
%          a(6,8)=max(min(randn/10+0.4,0.6),0.2);
%          a(10,4)=max(min(randn/10+0.4,0.6),0.2);
%          a(1,3)=max(min(randn/10+0.4,0.6),0.2);
%          a(7,5)=max(min(randn/10+0.4,0.6),0.2);
%          
%      end;
%      
%          
%      
%      
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     
%     
%     
%     
% 
% 
%   %  a(I,I+4)=max(min(randn/10+0.4,0.6),0.2);   %this makes and edge between ROI1 and ROI5 but in the opposite direction to avoid a cycle but to close the graph.
%                                                 %R2->R1;R3->R2;R4->R3;R5->R4;R5->R1;
%                                    
%    
%     
%     if Cyclic == 1        % cyclic connections
%       a(I,I+4)=0;         %make the direction between R1 and R5 zero.
%       a(I+4,I)=max(min(randn/10+0.4,0.6),0.2);  %and then make the direction between R5 and R1 positive to close the cycle.
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
     %assign values to the A connections (just positive values)
    a=Amatrix.*max(min(randn(size(Amatrix,1))/10+0.4,0.6),0.2);
    
    %include the -1 decays for each node
    a=a-eye(size(a));
    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%
    
    %Uncomment according to the case
    
    %used for structure2_contr_p2n6
    %a(3,4) = +0.2;
    %a(4,3) = -0.6;
    
    %used for structure2_contr_p6n2
    %a(3,4) = +0.6;
    %a(4,3) = -0.2;
    
    
    %used for Structure2_contr_c4
    %a(3,4)= +0.4;
    %a(4,3)= -0.4;
    
    %used for Structure2_amp_c4
    %a(3,4)= +0.4;
    %a(4,3)= +0.4;
    
    %a(6,7)=-a(6,7); %used for Structur3_amp_cont
    %a(4,5) = -a(4,5); %used for Structure3_cont_amp
    
    %a(4,3)=-a(4,3);  %used for Structure2_cont
    
    %a(4,5) = -a(4,5); %used for Structure1_cont
    %a(2,3) = -a(2,3); %used for Structure1_cont_b
    
    %a(4,3) = -a(4,3); %used for Structure5_cont
    
    %a(5,2) = -a(5,2); %used for Structure4_cont_amp
    %a(7,8) = -a(7,8); %usef for Structure4_amp_cont
    
    %Structure 4, heterogeneous subjects. 
    %For each subject add edges on top of the group structure
    %Add n random edges, n is randomly distributed, from 1 to a(numedges)/2
   % ne = randi(floor(size(a,1)/2)); 
   % u = find(a==0); %find candidates to assign edges, ie, zero entries of A
   % p = datasample(u,ne,'replace',false); %select the edges that will be added to the subject
    %assign values to the new edges
   % a(p)=max(min(randn(size(p,1),1)/10+0.4,0.6),0.2);
    
        
        
        
    
    
    
    

    c=eye(size(a,1));      %%%%%% external inputs. 1 input for each ROI

    %%%%%this is to simulate just 1 input in the first node
%     c = zeros(size(a,1));
%     c(1,1) = 1;

%     
%     
%     
%     if StrongLinks == 1        
%       a=a + 1.25 * a .* (a>0); % strong (mean 0.9)
%       if OneExternal > 0        % pretty much only one external input on node 1. And the other inputs very weak 0.1
%         c(2,2)=.1; c(3,3)=.1; c(4,4)=.1; c(5,5)=.1;
%       end
%     end
%     
%     
%     
%     
%     if SharedInputs == 1        % shared inputs with chance 20% ie. make the C matrix different from diagonal.
%       for i=I:I+4 
%           for j=I:I+4
%                 if j~=i && rand>0.8  %if it is not a auto-connection (i=j) and rand>0.8
%                     c(i,j)=0.3;    %make share inputs c(1,3) means that INPUT1 affects ROI3...
%                 end
%             end ; 
%       end
%     end
%   end
% 
%   %%%%% cross-subnet connections, MAKE CONNECTIONS ACROSS SUBNETS.  CHECK
%   %%%%% THAT THE ROI ARE SELECTED ASSUMING A NUMBER OF ROIS (N). IF WE SET
%   %%%%% N=40 THEN THE A(1,48) WILL NOT APPEAR....THIS HAS TO BE BUILD
%   %%%%% HAVING IN MIND THE NUMBER OF ROIS.
%   % if Nsubnets==2  
%    %  a(3,8)=max(min(randn/10+0.4,0.6),0.2);  %MAKE A CONNECTION BETWEEN ROI3 (SUBNET1) AND ROI8 (SUBNET2)
%    %end
%   
%   
%   
%   if Nsubnets==3                             %MAKE CONNECTION BETWEEN ROI3-ROI8; ROI3-ROI13; ROI8-13   WHEN NSUBNETS=3
%     a(3,8)=max(min(randn/10+0.4,0.6),0.2);
%     a(3,13)=max(min(randn/10+0.4,0.6),0.2);
%     a(8,13)=max(min(randn/10+0.4,0.6),0.2);
%   end
%   
%   if Nsubnets==10
%     a(3,8)=max(min(randn/10+0.4,0.6),0.2);
%     a(8,13)=max(min(randn/10+0.4,0.6),0.2);
%     a(13,18)=max(min(randn/10+0.4,0.6),0.2);
%     a(18,23)=max(min(randn/10+0.4,0.6),0.2);
%     a(3,23)=max(min(randn/10+0.4,0.6),0.2);
%     a(28,33)=max(min(randn/10+0.4,0.6),0.2);
%     a(33,38)=max(min(randn/10+0.4,0.6),0.2);
%     a(38,43)=max(min(randn/10+0.4,0.6),0.2);
%     a(43,48)=max(min(randn/10+0.4,0.6),0.2);
%     a(28,48)=max(min(randn/10+0.4,0.6),0.2);
%     a(3,28)=max(min(randn/10+0.4,0.6),0.2);
%   end
% 
%   
%   
%   
   N=size(a,1); Nreal=N;  % used to specify when there are additional nodes to be deleted later
%   if Modulation > 0    %check cases for those that have a modulation >0
%       if Modulation < 3     % original modulations used in netsim paper  (??????)
%                 a=[a a*0]; 
%                 a=[a;a*0];   
%                 a=a + 1.25 * a .* (a>0);   %why it changes the size of the matrix to 2*N
%                 N=N*2;   %double the size of the matrix
%                 c=eye(N);   
%         
%         
%                 for j=2:N/2  %the cycle runs from 2 to the original value of N 
%                     c(j,j)=.3;  %the value of the connection of the C matrix equals 0.3. All except the first ROI1-INPUT1
%                 end;
% 
% 
%                 if Modulation == 1
%                   DOMODext=zeros(N,N,N);  %define a 3D matrix to define modulation of intrinsic connections (a). This is THE b matrix. 
%                   DOMODext(1,2,6)=-0.12;  %the connection between ROI1-ROI2 is modulated by INPUT6
%                   DOMODext(2,3,7)=-0.12;  %the connection between ROI2-ROI3 is modulated by INPUT7
%                   DOMODext(3,4,8)=-0.12; 
%                   DOMODext(4,5,9)=-0.12; 
%                   DOMODext(1,5,10)=-0.12;
%                 end
% 
%                 
%       elseif Modulation < 5    % new modulations for low-TR HCP work
%                 a=a + 0.5 * a .* (a>0);  % slightly stronger links
%                 c=0.6 * eye(N);          % slightly weakened external inputs. 0.6 instead of 1.
%                 a(3,28)=-a(3,28);        % make link between two halves negative. I ASSUME HERE IT IS CONNECTING TWO SUBNETS.
%         
%                 if Modulation == 3
%                   a=[a zeros(N,4); zeros(4,N+4)]; %EXPAND THE A MATRIX WITH 4 NODES
%                   N=N+4; %WHY IT SUMS 4 NODES.
%                   c=0.6*eye(N);
%                   
%                   DOMODext=zeros(N,N,N); 
%                   modval=-0.085;
%                                                 %same modulation matrix as
%                                                 %above
%                   DOMODext(11,12,51)=modval; 
%                   DOMODext(12,13,51)=modval; 
%                   DOMODext(13,18,51)=modval;
%                   DOMODext(1,2,52)=modval; 
%                   DOMODext(2,3,52)=modval; 
%                   DOMODext(3,4,52)=modval; 
%                   DOMODext(4,5,52)=modval; 
%                   DOMODext(1,5,52)=modval;
%                   DOMODext(46,50,53)=modval; 
%                   DOMODext(47,48,53)=modval; 
%                   DOMODext(48,49,53)=modval; 
%                   DOMODext(49,50,53)=modval;
%                   DOMODext(43,48,54)=modval; 
%                   DOMODext(48,49,54)=modval;
%                 end
% 
%       else  % (Mod = 10)
%         DOMODint=zeros(N,N,N);  
%         DOMODint(6,7,5)=2; 
%         DOMODint(1,2,4)=2;
%       end
%   end
% 
%   
%   
%   
%   
%   %check conditions 101 onwards for value of the parameter ThreeNodes
%   if ThreeNodes > 0  %threenodes define graphs of only three nodes and test some cases.
%     std_BOLD_noise=0.1;  %define the noise
%     c=[1 0 0; 0 1 0; 0 0 1];  %each ROI is affected by an exclusive INPUT. C matrix.
%     
%     switch ThreeNodes
%       case {1}     % A -> B -> C
%         a=[-1 0.5 0;  0 -1 0.5;  0 0 -1];
%       case {2}     % A -> B <- C
%         a=[-1 0.5 0;  0 -1 0;  0 0.5 -1];
%       case {3}     % A <- B -> C
%         a=[-1 0 0;  0.5 -1 0.5;  0 0 -1];
%      
%       case {4}     % A -> B -> C   (negative connection to C)
%         a=[-1 0.5 0;  0 -1 -0.5;  0 0 -1];
%       case {5}     % A -> B <- C   
%         a=[-1 0.5 0;  0 -1 0;  0 -0.5 -1];
%       case {6}     % A <- B -> C   
%         a=[-1 0 0;  0.5 -1 -0.5;  0 0 -1];
%       
%       case {7}     % A -> B -> C   (only A has external input. SEE THE C MATRIX.)
%         c=[1 0 0; 0 0 0; 0 0 0];
%         a=[-1 0.5 0;  0 -1 0.5;  0 0 -1];
%       
%       case {8}     % A -> B <- C   (only A,C have external input)
%         c=[1 0 0; 0 0 0; 0 0 1];
%         a=[-1 0.5 0;  0 -1 0;  0 0.5 -1];
%     
%       case {9}     % A <- B -> C   (only B has external input)
%         c=[0 0 0; 0 1 0; 0 0 0];
%         a=[-1 0 0;  0.5 -1 0.5;  0 0 -1];
%       
%       case {10}     % A -> B -> C  (with a little negative from A to C)
%         a=[-1 0.5 -0.1;  0 -1 0.5;  0 0 -1];
%       
%       case {11}     % A -> B <- C
%         a=[-1 0.5 -0.1;  0 -1 0;  0 0.5 -1];
%       
%       case {12}     % A <- B -> C
%         a=[-1 0 -0.1;  0.5 -1 0.5;  0 0 -1];
%     end
%   end
% 
%   
%   
%   if NegativeConnections > 0
%     a(2,3)=-a(2,3);
%   end
%end   %THIS IS THE END OF THE DEFINITION OF THE GRAPH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we create the first half of the data. ie. n=nsub/2.
  aa(n,:,:)=a; 
  cc(n,:,:)=c;   %aa is a 3d matrix. the 1d contains the index of subject and 2d&3d contain connectivity matrix A. 
                 %cc is a 3d matrix. the 1d contains the index of subject and 2d&3d contain input connectivity matrix C.   
 
%%%%%% setup second-half of the subjects
%   if GroupTwoDiff==1     %CHECK CASE {21}
%     a(2,3)=a(2,3)*.5;  %DIVIDE THE VALUE OF A(2,3) CONNECTION BY HALF.
%   end
  
%   aa(n+Nsub/2,:,:)=a;
%   cc(n+Nsub/2,:,:)=c;

  
end

static=strcat('thesis_simulations/',outputfolder);
mkdir(static);
%save the A matrices of all the sessions in a 3D matrix
%first dimension is the number of sessions, 2nd and 3d is the graph.
save(strcat(static,'/Amatrices.mat'),'aa');


save(outname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nma_rsns_thesis;   % call the DCM-for-FMRI forward model to generate BOLD timeseries




%MAKE THE POOL WITH ALL THE DATA







save(outname);

%CHECK NMA_RSNS_RUBENNOTES FOR EXACT DETAILS OF THE FUNCTION.




 if Nreal < N     % remove the extra nodes that were used to create modulations
   N=Nreal;
   ts=ts(:,1:N);  %TS CONTIANS YS2 WHICH CONTAINS THE BOLD (Y) SIGNALS OF ALL THE ROIS.
   aa=aa(:,1:N,1:N);  %ALL THE SUBJECTS AND THE A INSTRINSIC CONNECTIVITY MATRIX
 end

%[vals i]=max(poop,[],1); std(i)*tr      %show HRF lag variability (need to set tr to 0.1 first)

%%%%% show arrow-plot of model
%grot=cell(1,N+P);
%for i=1:N, grot{i}='ellipse'; end;    for i=1:P, grot{i+N}='rect'; end;
%plotmodel([ [a+eye(size(a)) ; c'] zeros(N+P,P) ]' , 1:N+P, 'nodeshapes', grot);



 T=length(ts);
 TT=T/Nsub;
 tsn=ts; %MAKE A VERSION OF TS THAT WILL NOT BE FILTERED OR DEMEANED
% 
% 
% 
% 
% %%%% highpass and demean ts
 %[Bb,Ab] = butter(4, 0.03, 'high');      % what we used in the paper - too aggressive for low-TR
 [Bb,Ab] = butter(4, 0.03*tr/3, 'high');  % HPcutoff=200s regardless of TR   
% 
% 
% 

 %AUXILIARY TO SAVE THE BOLD SIGNALS OF EACH SUBJECT 
 for s=1:Nsub 
   ts1((s-1)*TT+1:s*TT,:)=demean(filter(Bb,Ab,ts((s-1)*TT+1:s*TT,:)),1);
   BOLDdemean1=ts1((s-1)*TT+1:s*TT,:);
   dlmwrite(strcat(static,'/BOLDdemefilt1_',num2str(s),'.txt'),BOLDdemean1,'delimiter','\t','precision','%.3f');
   %save(strcat(static,'/50_BOLDdemefilt1_',num2str(s),'.txt'),'BOLDdemean1','-ascii');
 end
 
 %THIS IS THE ORIGINAL FORM IN THE SMITH CODE. IT DEMEANED AND FILTER THE
 %FULL TS MATRIX FOR ALL THE SUBJECTS
%  for s=1:Nsub
%   ts((s-1)*TT+1:s*TT,:)=demean(filter(Bb,Ab,ts((s-1)*TT+1:s*TT,:)),1);
% end





%%%%%% prewhiten
%%%% ar1 = mean( sum( ts(2:T,:) .* ts(1:T-1,:) ) ./ sum(ts .* ts ) )
%%%% clear pwV;
%%%% for i = 1:TT
%%%%   for j = 1:TT
%%%%     pwV(i,j)=ar1^abs(i-j);
%%%%   end
%%%% end
%%%% hrfinv=inv(chol(pwV));
%%%% for s=1:Nsub
%%%%   ts((s-1)*TT+1:s*TT,:)=demean( hrfinv * ts((s-1)*TT+1:s*TT,:));
%%%% end



%CHECK CONDITION 11 AND 12
% if BadROIs == 1
%   MixingFraction=0.2;
%   [yy,ii]=sort(rand(1,N));
%   %ts= (1-MixingFraction)*ts + MixingFraction*ts(:,N:-1:1);
%   
%   ts1= (1-MixingFraction)*ts1 + MixingFraction*ts1(:,ii); %THIS IS FOR THE DEMEAN-FILTERED DATA
%   for s=1:Nsub 
%    BOLDdemfiltBADROI1=ts1((s-1)*TT+1:s*TT,:);
%    save(strcat(static,'/50_BOLDdemfiltBADROI1_',num2str(s),'.txt'),'BOLDdemfiltBADROI1','-ascii');
%   end
% 
%  tsn= (1-MixingFraction)*tsn + MixingFraction*tsn(:,ii); %THIS IS FOR THE ONLY NOISE DATA
%   for s=1:Nsub 
%    BOLDnoiseBADROI1=tsn((s-1)*TT+1:s*TT,:);
%    save(strcat(static,'/50_BOLDnoiseBADROI1_',num2str(s),'.txt'),'BOLDnoiseBADROI1','-ascii');
%  end
% 
% 
% 
% 
% %THE NUMBER OF SUBJECTS (NSUB) SHOULD ALWAYS BE EQUAL OR GRATER THAN THE NUMBER OF
% %ROIS (N)
% elseif BadROIs == 2
%   MixingFraction=0.2;
%   
%   for s=1:Nsub
%     for nn=1:N
%       ss=s+nn; if ss>Nsub, ss=ss-Nsub; end;
%       ts1((s-1)*TT+1:s*TT,nn)=(1-MixingFraction)*ts1((s-1)*TT+1:s*TT,nn) + MixingFraction*ts1((ss-1)*TT+1:ss*TT,nn);    
%     end
%     BOLDdemfiltBADROI2=ts1((s-1)*TT+1:s*TT,:);
%     save(strcat(static,'/50_BOLDdemfiltBADROI2_',num2str(s),'.txt'),'BOLDdemfiltBADROI2','-ascii');
%   end
%   
%   for s=1:Nsub
%     for nn=1:N
%       ss=s+nn; if ss>Nsub, ss=ss-Nsub; end;
%       tsn((s-1)*TT+1:s*TT,nn)=(1-MixingFraction)*tsn((s-1)*TT+1:s*TT,nn) + MixingFraction*tsn((ss-1)*TT+1:ss*TT,nn);    
%     end
%     BOLDnoiseBADROI2=tsn((s-1)*TT+1:s*TT,:);
%     save(strcat(static,'/50_BOLDnoiseBADROI2_',num2str(s),'.txt'),'BOLDnoiseBADROI2','-ascii');
%   end
%   
%   
%  
% end

%%%%%%%%%%%%%%%%%  estimate nulls
% 
% clear global cm opcm pcm names; global cm opcm pcm names ConfoundGlobalMean tr testrun NormICOV; cm=0; names=cell(1,1);
% net_measures(reshape(ts(:,2),TT,Nsub),N);   % using node 2 as it is more 'representative' than node 1
% save(outname);
% % figure; boxplot(cm);

%%%%%% optionally add in global mean timecourse confound
% if ConfoundGlobalMean > 0
%   ts1=add_confound(ts1,TT,ConfoundGlobalMean);
%   tsn=add_confound(tsn,TT,ConfoundGlobalMean);
%   
%   for s=1:Nsub 
%    BOLDdemfiltConfoundGlobalmean=ts1((s-1)*TT+1:s*TT,:);
%    save(strcat(static,'/50_BOLDdemfiltConfoundGlobalmean_',num2str(s),'.txt'),'BOLDdemfiltConfoundGlobalmean','-ascii');
%   end
% 
%  
%   for s=1:Nsub 
%    BOLDnoiseConfoundGlobalmean=tsn((s-1)*TT+1:s*TT,:);
%    save(strcat(static,'/50_BOLDnoiseConfoundGlobalmean_',num2str(s),'.txt'),'BOLDnoiseConfoundGlobalmean','-ascii');
%  end
% end




%%%%%%%%%%%%%%%%%  estimate real effects

% clear ozz zzz;
%  for s=1:Nsub
% %   s
%    net_measures(ts((s-1)*TT+1:s*TT,:),N);  %CHECK THIS FUNCTION IN THE FOLDER. IT IS RELATED TO THE CAUSAL SEARCH ALGORITHMS TESTED IN THE PAPER.
% %   ozz(s,:,:,:)=opcm;
% %   if Nsubnets > 3
%     save(outname);
%   end
% end
% Nmeasures=size(ozz,4);
% save(outname);

% use_null;
%save(outname);

%%%%%%% make matrix-image plots 
%%%%% one-group t-test on connections strengths
%[H,P,CI,STATS]=ttest(reshape(zzz,Nsub,N*N*Nmeasures));
%figure;  maxz=8;   sb=ceil(sqrt(Nmeasures));
%zzzz=reshape(STATS.tstat,size(opcm)); zzzz(isnan(zzzz))=0;
%for i=1:Nmeasures
%  subplot(sb,sb,i); imagesc(zzzz(:,:,i),[-maxz maxz]); set(gca,'XTick',[],'YTick',[]); xlabel(names(i));
%end
%%%%% mean x-x' causality differencing
%zzzz=reshape(mean(reshape(zzz,Nsub,N*N*Nmeasures)),size(opcm));
%figure;  maxz=3;   sb=ceil(sqrt(Nmeasures));
%for i=1:Nmeasures
%  subplot(sb,sb,i); imagesc(squeeze(zzzz(:,:,i))-squeeze(zzzz(:,:,i))',[-maxz maxz]); set(gca,'XTick',[],'YTick',[]); xlabel(names(i));
%end

