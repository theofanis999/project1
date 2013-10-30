%% Comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cellular simulator (cellsim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 22.05.2013 01:25
% started
% based on AAMC205.m
% (same as AAMC204.m but with 6 states, not 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 23.05.2013 01:27
% reached just before the...
% utility function calculation
% added simulation time measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 08.10.2013 01:37
% Included decision making and..
% utility function calculation
% Final calculations &
% a VERY good review still pending
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 09.10.2013 01:27
% Sanitized global & iteration initialisation
% Stopped before Adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12.10.2013 01:29
% Finished transmitting node loop...
% ... including nested receiving node loop
% Renewed ENERGYSPENTT calculation
% Started new node loop with policies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 16.10.2013 02:15
% Reviewed till result calculations:
% "Set (beta, mode) at the end of this time slot"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 22.10.2013 02:11
% The beta now gets correctly adapted
% The mode does not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 24.10.2013 02:06
% The mode does get adapted albeit very little
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 26.10.2013 03:38
% 1. From now action_choice_weight=<1 & is calc'd like this....
% action_choice_weight(i,pindex) = POLICYPOPUL(i,pindex)/t;
% 2. In the calculation of the UUNCH now....
% UUNCH(i,t) = (action_choice_weight(i,1))*METRIC7(1,i)+...
% ... rather than...
% UUNCH(i,t) = (action_choice_weight(i,1)/t)*METRIC7(1,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 29.10.2013 01:06
% Added information freshness measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 31.10.2013 00:37
% Different OST;
% if U>UUNCH
% start counting time
% when time τ > unknownfactor/e
% if AVEUT(τ) > AVEUT(t) for all t = 1...t-1
% then OST is fullfilled!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 31.10.2013 01:18
% Put metrics: "freshinfect", efficiency, infection cost, ave. utility
% Put plots
% Changed name to "cellsim"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial cleanup
% function []=AAMC(N,T,iterations,aggressiveness,noise_model,adaptivity,metric,orig_infected,expiry);
clc
clear
%% Simulation banner
yy=clock;
starthour=num2str(yy(4));
startminu=num2str(yy(5));
startseco=num2str(yy(6));
starthour,startminu,startseco
disp('===RUNNING===')
disp('Our policy-choice-based model')
%% Global initialization

% Basic parameters
N = 50; % number of nodes
T = 100; % duration
iterations = 3; %algorithm's iterations

netdens = 0.85; % network density

adaptive = 1; % value 1 means adaptive scheme, value 0 means static beta

orig_infected = 2;

aggressiveness = 3;   % Used in Utility function averaging
if aggressiveness == 1
    aggressistring = 'aggressive';
else if aggressiveness == 2
        aggressistring = 'conservative';
    else
        aggressistring = 'lazy';
    end
end

% b0 = BETAS(10); % original beta
numberofbetas = 25;
deltabeta = 0.02;

% Mode
MODES = [1 2 3 4 5 6];
startingmode = 3;

% Beta
index1=1:numberofbetas;
BETAS=0.20+index1*deltabeta;
deltaoverbeta = 0.1;
resetsize=7;
betasum = 0;

% Noise
loss = 0.8; % this is essentially BER (attention, not PER!!)
signal = 30e+015;
n0 = 15;
sigma = 7;
nb = 45; %amplitude of noise burst
trainprob = 0.2; %probability that a noise burst occurs
floor = 30; % floor in case of awgn
% noise_model = 'awgn';
noise_model = 'sos';

% Infections
SUCCESSFUL_INF_COUNT_GLOBAL = zeros(iterations,N);
SUCCESSFUL_INF_COUNT = zeros(N,T,iterations);

% Info freshness
expiry = T/4; % infecting information lifetime

%Channel state
STATES = [1 2 3 4 5 6];
startingchannelstate = 3;

% Energy cost
bitTxcost = 720; %bit transmission energy cost in nJ
bitRxcost = 110; %bit reception energy cost in nJ
instructioncost = 4;
pktsize = 50; % packet size in bits
maxcostsingle = 1.8*pktsize*bitTxcost; % this is max cost for a SINGLE transmission

ENERGYSPENT= zeros(N,T,iterations);    % energy spent UNTIL the specified timeslot
ENERGYSPENTTX= zeros(N,T,iterations);
ENERGYSPENTRX= zeros(N,T,iterations);

ENERGYSPENTT= zeros(N,T,iterations);    % energy spent during the ONE specified timeslot
ENERGYSPENTTXT= zeros(N,T,iterations);
ENERGYSPENTRXT= zeros(N,T,iterations);

ENORM = zeros(1,T);

% Other
weight = 1.0;
actualtime = 0;
berrors = 0;
bergood = 0;
unknownfactor = 15; % 100% arbitrary; only limitation: unknownfactor/2.718 > 1

changes_count = zeros(N,T,iterations);

% Channel noise model parameters acc. to Yannaki
GN(1)=7.9932;
GN(2)=3.4998;
GN(3)=1.6883;
GN(4)=0.6644;
GN(5)=0.3756;
GN(6)=0.09;
     
AN(1)=274.7229;
AN(2)=90.2514;
AN(3)=67.6181;
AN(4)=50.1222;
AN(5)=53.3987;
AN(6)=35.3508;
    
GPN(1)=0.7026;
GPN(2)=1.2865;
GPN(3)=2.4959;
GPN(4)=5.89;
GPN(5)=10.5896;
GPN(6)=39.6278;
GPN(7)=10; % programming trick; the value (10) is arbitrary but non-void ;-)
    
%Overhead for each encoding mode
OVERHEAD(1)= 2; % overhead factor= 1/(coding rate); code rate shows the useful data proportion, overhead shows the overhead(!) proportion
OVERHEAD(2)= 2;
OVERHEAD(3)= 1.3333; % 4/3;
OVERHEAD(4)= 1.7778; % 16/9
OVERHEAD(5)= 1.3333;
OVERHEAD(6)= 1.3333;
SNRLIMITS=GPN;

for index2 = 1:7
    NOISELIMITS(index2)=signal/SNRLIMITS(index2);
end

MODE_0 = zeros(1,iterations);   % starting mode average
DUPT=zeros(T,iterations);

POLICYPOPUL_ITER = zeros(iterations,4);  % because policyindex = 1:4; policy popularities
TC2FC = 0;
T2FC = 0;
E2FC = 0;
TC29C = 0;
T29C = 0;
E29C = 0;
%% Bursty noise
% Flapping function; it's a square function
FLAP = zeros(1,T);
FLAPJITTER = 1+rand(1,T);
 ft = 3; % flapping duration
 for f = 1:T
     if mod(f,25) == 0
         for ff = f:f+ft+5;  % guard time
             FLAP(ff+ round( FLAPJITTER(f) ) ) = 2;
         end
     else
         % do nothing
     end
 end
%% Import from SoS model
open transMatrix_nfsr_08.mat;
TM = ans.TM_ave;
clear ans;
%% Iteration loop
for iteration = 1:iterations
    iteration
%% Iteration-dependent initialization

% action_choice_weight = zeros(N,4)+0.25; % has been moved inside the temporal loop 
changes_count = zeros(N,T); % indices mean (node,chosen_policy);


% Noise

NOISE = zeros(N,N) + n0; % ATTENTION!!! DOES THIS GO INSIDE THE TEMPORAL LOOP?
SNR = zeros(N,N,T);
NOISET = zeros(N,N,T) + n0;
NOISEPURET = zeros(N,N,T) + n0;
RECENTNOISE = zeros(N,N,T) + n0;
AVENOISET = zeros(N,N,T) + n0;
S2 = zeros(N,T);
nodenoise = zeros(1,N);
nodenoise = zeros(T); % dummy initialization of nodenoise < ===== ATTENTION!! DIMESIONS??!
avenoise = zeros(1,N);
% PER = zeros(N);
% BER = zeros(N) + loss;
PER = zeros(N,N);
NOISEPURE = zeros(N,N) + n0;
PERCALC = zeros(4,4);

% For the initialization of noise for this interation; noise is NOISET(node,othernode,time)
for index3 = 1:N;
    for index4 = 1:N;
    randostate(index3,index4) = rand;
    if randostate(index3,index4) <= 0.0835 % so that (length(STATES)*randostate(node)) >= 0.5
        randostate(index3,index4) = 0.084; % int16(.) will safely give 1
    end
    end % index4
end % index3

% Beta
for index7 = 1:N
    randobeta(index7) = rand;
    if randobeta(index7) <= 0.02 % (length(BETAS)*randobeta(node))
        randobeta(index7) = 0.022;
    end
end % index7

B = BETAS(int16(length(BETAS)*randobeta)); % random original beta
% B = zeros(1,N) + 0.4;

BTEMP = B;    
BFACTOR = zeros(1,N);

for node = 1:N
    for t = 1:T
        BDET(node,t) = B(node); % BDET(= beta detailed) like B(node) but keeps temporal evolution too
    end
end

for a = 1:N;        % number of nodes
    for b = 1:4;    % number of possible actions
        BCALC(b,a) = B(a);
    end
end

UP = zeros(1,N); % number of times the beta has been forced up when minimum SNR-limit reached - ACCORDION-RELATED
DOWN = zeros(1,N); % number of times the beta has been forced down when maximum SNR-limit reached - ACCORDION-RELATED


% Mode
for index8 = 1:N
    randomode(index8) = rand;
    if randomode(index8) <= 0.0834 % so that (length(MODES)*randomode(node)) >= 0.5
        randomode(index8) = 0.0834;
    end
end % index8

MODE = MODES(int16(length(MODES)*randomode)); % random original mode
% MODE = zeros(1,N) + 3;
MODE_0(iteration) = mean(MODE);

MODETEMP = MODE; % arbitrary; value plays no role
MODEFACTOR = zeros(1,N);

MODEDET = zeros(N,T) + startingmode; % MODEDET(= mode detailed) like MODE(node) but keeps temporal evolution too  
MODECALC = zeros(4,N) + 3;

DOWNMODE = zeros(1,N);
UPMODE = zeros(1,N);

% A is the nodes' infection state vector; A is 1xN

A = zeros(1,N);
ATEMP = zeros(1,N);

% Simple case with 1 originally infected node

% A(1) = 1;
% ATEMP(1) = 1;

for index12 = 1:orig_infected
    index13 = max(1, round( abs(randn(1))*(N+1) ) ); % set random nodes to infected state
    A(index13) = 1;
end
ATEMP = A;

A_1 = sum(A);

% Information Freshness
% Freshness = 1 means just infected
% Freshness = 0 is not taken into account (not infected)

FRESH = zeros(1,N);
FRESHTEMP = zeros(1,N);
FRESH(1) = 1;
FRESHTEMP(1) = 1;

% Transmission and reception counts
TX = zeros(1,N);
RX = zeros(1,N);
TXTEMP = zeros(1,N);
RXTEMP = zeros(1,N);

% Errors (=corrupt packets)
ERR = zeros(1,N);
ERRTEMP = zeros(1,N);

ERRRATE = zeros(1,N);
ERRRATETEMP = zeros(1,N);

% Duplicates
DUP = zeros(1,N);
DUPTEMP = zeros(1,N);

DUPRATE = zeros(1,N);
DUPRATETEMP = zeros(1,N);
% DUPT=zeros(T,iterations);  % moved outside iterations loop!

% Infections
SUCCESSFUL_INF_COUNT = zeros(1,N);

% Mode lifetime; sliding window CAN BE CHOSEN to assume this value
WINDOW = zeros(1,N) + 1;

% Bit count
OVERHEADBITSENT= zeros(1,N);
BITSENT= zeros(1,N); % Number of bits sent by a node. So that when I sum up for all nodes I get the total bits sent.

% Energy cost
% ENERGYSPENT= zeros(N);
% ENERGYSPENTTX= zeros(N);
% ENERGYSPENTRX= zeros(N);

% Channel state
for index5 = 1:N;
    for index6 = 1:N;
        CHANSTATE(index5,index6) = STATES(int16(length(STATES)*randostate(index5,index6))); % random original channel state; initialize matrix with channel states (->SNR) as perceived by EACH NODE!!      
    end % index6
end % index5

CHANSTATETEMP = CHANSTATE;
CHANSTATET = zeros(T,N,N) + startingchannelstate;

UPIFSMC = zeros(N,N); % state increment count for channel status (input FSMC)
DOWNIFSMC = zeros(N,N); % state increment count for channel status (input FSMC)
SWITCHIFSMC = zeros(N,N) + UPIFSMC + DOWNIFSMC;
UPIFSMCTEMP = zeros(N,N);
DOWNIFSMCTEMP = zeros(N,N);

% Other
PERCEIVEDADJ = zeros(N,N); % PERCEIVED adjacency matrix! = KNOWN neighbors
POLICYPOPUL = zeros(N,4);  % there are 4 policies-possible actions

PINF = zeros(4,N);         % Four (4) policies to handle beta and mode
PINFAUX = zeros(4,N);
ENERGYSPENTCALC = zeros(4,N);
ENERGYSPENTTXCALC = zeros(4,N);
ENERGYSPENTRXCALC = zeros(4,N);

% SSS = zeros(N,4)+0.25;  % renamed, see below
action_choice_weight = zeros(N,4) + 0.25;  % indices mean (node,chosen_policy) % renamed from SSS(node,chosen_policy) 

% U.F.O. UFO
for t=1:T;
    g(t)=0.5*randn(1);
end

% Utility function
RECENTUT = zeros(N,T);
AVEUT = zeros(N,T);
U = zeros(N,T);
UTEMP = zeros(N,T);

% Measurements
DUPRATE_NO_TIE_INS = 0.01*abs(randn(N,T));
NO_TIE = 1+0.01*abs(randn(N,T));
DUPRATE_NO_TIE = zeros(N,T);

% Decision matrices initialization
% These values are indifferent; just for the matrices to be non-void
% These matrices are always changed; every node, timeslot & iteration(?)


%% Adjacency matrix
% Put it here for no mobility
ADJ = rand(N,N);
for m = 1:N
    for n = 1:N
        if ADJ(m,n)>netdens;
            ADJ(m,n) = 1;
            ADJ(n,m) = 1;
        else
            ADJ(m,n) = 0;
            ADJ(n,m) = 0;
        end
    end
end

% Neighbor count, attributed to beaconing
NC = sum(ADJ);
for node=1:N;
    NEIGHCOUNT(node) = NC(node);  % <-- ALL nodes, assumed known from beacons!
end

%% Temporal loop - Epidemic dissemination loop
for t = 1:T
    actualtime = actualtime + 0.001;
    action_choice_weight = zeros(N,4)+0.25;
    SUCCESSFUL_INF_COUNT_TEMP(node) = 0;  % calculated but not really used
    FRESH = FRESH + 1;   % info freshness ageing; if expiry threshold crossed, it is cured (see later)
    
    %% Markov Chain refresh - keep inside the big temporal loop to avoid giant matrices:
    %  CHANSTATE,DOWNIFSMC, UPIFSMC, & CHANSTATETEMP,DOWNIFSMCTEMP,UPIFSMCTEMP(???)
    
    % Refresh channel state markov chain for all LINKS (NOT nodes) start
    for index9 = 1:N
        for index10 = 1:N
            a = 1:6;             % possible channel states - modes
            weight = TM(CHANSTATE(index9,index10),:); %# corresponding weights
            results = 1;              % how many numbers to generate
            R = a( sum( bsxfun(@ge, rand(results,1), cumsum(weight./sum(weight))), 2) + 1 );
            CHANSTATETEMP(index9,index10) = R; % thought to make it CHANSTATETEMP(node20,node21,time) and move it out of the temporal loop
            if CHANSTATETEMP(index9,index10) > CHANSTATE(index9,index10)
                UPIFSMCTEMP(index9,index10) = UPIFSMCTEMP(index9,index10) + 1;
            end
            if CHANSTATETEMP(index9,index10) < CHANSTATE(index9,index10)
                DOWNIFSMCTEMP(index9,index10) = DOWNIFSMCTEMP(index9,index10) + 1;
            end
        end  % node21
    end % node20
    % No real reason to keep the CHANSTATETEMP, -TEMP, etc
    CHANSTATE = CHANSTATETEMP;
    DOWNIFSMC = DOWNIFSMCTEMP;
    UPIFSMC = UPIFSMCTEMP;
    SWITCHIFSMC = DOWNIFSMC + UPIFSMC;
    CHANSTATET(t,:,:) = CHANSTATE;

%% Transmitting node-related initializations

for node = 1:N  % (transmitting) node
    
%     Channel noise calculation
    for othernode = 1: N % othernode
        if CHANSTATE(node,othernode)~=7
            
            if strcmp(noise_model,'awgn') == 1
                NOISE(node,othernode) = floor+randn(1);
            end % if noise_model STRCMP awgn
            if strcmp(noise_model,'sos') == 1
                %             The following line gives sporadic negatives !!!!!!!!
                NOISEPURE(node,othernode) = (  NOISELIMITS(CHANSTATE(node,othernode)) + ( NOISELIMITS(CHANSTATE(node,othernode)+1)-NOISELIMITS(CHANSTATE(node,othernode))   )*randn(1,1)  ); % noise randmoly between limits
                NOISE(node,othernode) = (1+FLAP(t)) * NOISEPURE(node,othernode); % noise randmoly between limits
                %             NOISE(node,othernode) = NOISEPURE(node,othernode); % no BURSTY noise
            end % noise_model = sos
            
        else
            if CHANSTATE(node,othernode)==7
                NOISE(node,othernode) = NOISELIMITS(7)*(1+rand);
            end
        end
        
        NOISET(node,othernode,t)=NOISE(node,othernode);
        NOISEPURET(node,othernode,t)=NOISEPURE(node,othernode);
        SNR(node,othernode,t)= signal/(NOISET(node,othernode,t)); % Populate signal-to-noise matrix: a value for each node and timeSLOT.
        S2 = mean(SNR,2);  % SNR conceived by a *node*: mean over all its channels with neighbors
        S3 = squeeze(S2);  % auxiliary: S2 is Nx1xT; S3 is NxT, hence plottable!
        M = 2^(MODE(node));
        BER(node,othernode) = 0.2*exp(-3*SNR(node,othernode,t)/(2*(M-1)));  % Where is this equation from??? BER never used again
    end % othernode
    
%     PER calculation - use Giannaki!
    for othernode=1:N
        if othernode ~= node
            %         if MODE(node) == 0
            %             MODE(node) = 3;
            %         end
            if SNR(node,othernode,t)>GPN(MODE(node))
                PER(node,othernode) = AN(MODE(node))*exp(-GN(MODE(node))*SNR(node,othernode,t));
            elseif (SNR(node,othernode,t)<GPN(MODE(node)) && SNR(node,othernode,t)>0)
                PER(node,othernode) = 1;
            end
            
            if PER(node,othernode)>1
                %disp('PER ALARM !!!');
                %disp(PER(node));
                berrors = berrors +1;
                PER(node,othernode) = 1;
            else bergood = bergood+1;
            end % if PER
        end % if othernode
    end % othernode
       
end  % (transmitting) node
% disp('end node');
PER2 = mean(PER,2);
%% Transmitting node loop

for i = 1:N   % transmitting node loop
    % Attempt to cure
%     if A(i) == 1  % if node # i is infected
%         cure = rand;
%         if cure < B(i)*deltaoverbeta;
%             A(i) = 0;
%         end
%     end
    
    % Cure if info freshness crosses threshold (info expires)
    if FRESH(i) > expiry
        A(i) = 0;
    end
    
    % Attempt to infect (with banner)
    if A(i) == 1   % if node # i is infected
        % disp('===NOW ATTEMPT TO INFECT===')
        lucky = rand;
        if lucky < B(i) % here we make the beta experiment ONCE for each node: broadcast!
            %% Receiving node loop
            for j = 1:N
                %             if A(i) == 1   % if node # i is infected
                if ADJ(i,j) == 1 && (i ~= j)
                    % 	                   lucky = rand;
                    % 	                   if lucky < B(i)  % here we make the beta experiment on a p2p basis!
                    % 	                       TXTEMP(i)=TX(i)+1;
                    %                          ENERGYSPENT(i)=ENERGYSPENT(i)+pktsize*OVERHEAD(MODE(i))*bitTxcost; % Tx energy MUST be calculated ONCE for every broadcast!!
                    %                          BITSENT(i) = BITSENT(i) + pktsize + pktsize*OVERHEAD(MODE(i));
                    %                          OVERHEADBITSENT(i)=OVERHEADBITSENT(i)+pktsize*OVERHEAD(MODE(i));
                    RXTEMP(j)=RX(j)+1;
                    
                    %                          Add a new known neighbor; if it already known the element (j,i) is already =1
                    PERCEIVEDADJ(i,j) = 1;       % i.e. node j becomes aware of its neighbor i
                    %
                    % the following line is Rx energy spent
%                     ENERGYSPENT(j)  =ENERGYSPENT(j)  +pktsize*OVERHEAD(MODE(i))*bitRxcost+2*instructioncost; %packet size decided by sender's code rate!
                    ENERGYSPENTRXT(j,t,iteration)=ENERGYSPENTRXT(j,t,iteration)+pktsize*OVERHEAD(MODE(i))*bitRxcost+2*instructioncost; %packet size decided by sender's code rate!
                    channel = rand;
                    if channel > PER(i,j) %loss - PER(i,j) NOW! i=src, j=dst   ***PER(i,j) % if pkt received healthy
                        if A(j)==0      % if the receiver is susceptible...
                            ATEMP(j)=1; % ...then it becomes infected
                            FRESHTEMP(j) = FRESH(i) + 1; % assumes sender's info freshness
                            SUCCESSFUL_INF_COUNT_TEMP(j) = SUCCESSFUL_INF_COUNT(j)+1;
                        else % A(j)==0
                            if FRESH(j) > FRESH(i) + 1       % if the receiver is infected & its info less fresh than the sender's...
                                FRESHTEMP(j) = FRESH(i) + 1; % ...it is re-infected, and freshnesh is that of the sender's
                            else
                                DUPTEMP(j)=DUP(j)+1; % received info is a duplicate - rejected; freshness not affected
                            end
                        end % A(j)==0
                    else
                        ERRTEMP(j)=ERR(j)+1;
                    end % loss
                    % 	                   end % lucky % here we make the beta experiment on a p2p basis!
                end % ADJ
                %             end % A(i)
            end % j
            

            ENERGYSPENTTXT(i,t,iteration) = ENERGYSPENTTXT(i,t,iteration) + pktsize*OVERHEAD(MODE(i))*bitTxcost;
            ENERGYSPENTT(i,t,iteration)   = ENERGYSPENTTXT(i,t,iteration) + ENERGYSPENTRXT(i,t,iteration); 
            TXTEMP(i)=TX(i)+1;   % maybe do this like the ENERGYSPENT?            
            BITSENT(i) = BITSENT(i) + pktsize + pktsize*OVERHEAD(MODE(i));    % maybe do this like the ENRGYSPENTTXT?
            OVERHEADBITSENT(i)=OVERHEADBITSENT(i)+pktsize*OVERHEAD(MODE(i));  % maybe do this like the ENRGYSPENTTXT?
            
            AUXMAT1 = cumsum(ENERGYSPENTRXT,2);
            AUXMAT2 = cumsum(ENERGYSPENTTXT,2);
            
            for index99 = 1:t
                ENERGYSPENTRX(i,t,iteration) = AUXMAT1(i,t,iteration);
                ENERGYSPENTTX(i,t,iteration) = AUXMAT2(i,t,iteration);
            end
            
            ENERGYSPENT(i,t,iteration) = ENERGYSPENTRX(i,t,iteration) + ENERGYSPENTTX(i,t,iteration);
            
            
        end % lucky % here we make the beta experiment ONCE for each node: broadcast!
    end % A(i)
        
%% Calculate new (beta, mode) for the next time slot; time slot-dependent results

if RX(i) == 0
    BTEMP(i) = B(i);
end %RX(i) == 0

if RX(i)~=0  % This loop is for protection in the next two lines
    % Some measurements
    ERRRATE(i)=ERR(i)/(RX(i));
    DUPRATE(i)=DUP(i)/(RX(i));
    ERRRATETEMP(i)=ERRTEMP(i)/RXTEMP(i);
    DUPRATETEMP(i)=DUPTEMP(i)/RXTEMP(i);  
    
    if adaptive == 1
        
        %% Policies
        
        % Calculate some auxiliary quantities
        
        if S2(i,:,t)>S2(i,:,max(1,t-1))                     % LOWER noise
            BFACTOR(i) = - min(B(i)/deltabeta,resetsize);   % transmit less (AIMD) and...
            MODEFACTOR(i) = 1;                              % ...encode lighter
        else BFACTOR(i) = 1;                                % HIGHER noise, transmit more and...
            MODEFACTOR(i) = - 1;                              % ...encode stronger
        end % S2
        
        %     PER2 = mean(PER,2);
        
        % Set the sliding window according to the configured aggressiveness
        
        if aggressiveness == 1
            swindow = 1;
        else if aggressiveness == 2
                swindow = min(3,t-1);
            else
                swindow = min(WINDOW(i),t-1);
            end
        end
        
        % Some averaging calculations over the sliding window
        
        RECENTNOISE(i,t-swindow:t) = NOISET(i,t-swindow:t);
        AVENOISET(i,t) = mean(RECENTNOISE(i,t-swindow:t));
        RECENTUT(i,t-swindow:t) = U(i,t-swindow:t);
        AVEUT(i,t) = mean(RECENTUT(i,t-swindow:t));
        DUPRATE_NO_TIE(i,t) = mean(DUPRATE_NO_TIE_INS(i,t-fix(0.2*swindow):t));
        
        % Policies follow:
        
        % Policy #1 - keep mode, change beta
        
        BCALC(1,i) = B(i) + BFACTOR(i)*deltabeta;
        
        % PINF(1,i) = (0.01*abs(randn(1)) + 1)*(1-BCALC(1,i)/10)*BCALC(1,i)*(1-PER(i,j)^NEIGHCOUNT(i)); %reviewed with X"   ***PER(i,j)
        PINF(1,i) = (0.01*abs(randn(1)) + 1)*(1-BCALC(1,i)/10)*BCALC(1,i)*(1-PER2(i)^NEIGHCOUNT(i)); %reviewed with X"   ***PER(i,j)
        if PINF(1,i) < 0
            %         disp('PINF 1 < 0!'); % It never happens, ok!
            PINF(1,i) = 0; % Just set it to minimum(=0) so that the old value is not retained
        end
        ENERGYSPENTCALC(1,i)=(pktsize*OVERHEAD(MODE(i))*bitTxcost)*BCALC(1,i);
        
        % Policy #2 - change mode, keep beta
        
        if (MODE(i) == 6 && MODEFACTOR(i) > 0) || (MODE(i) == 1 && MODEFACTOR(i) < 0)
            MODECALC(2,i) = MODE(i);
        else MODECALC(2,i) = MODE(i) + MODEFACTOR(i);
        end
        
        PERCALC(2,i) = AN(MODECALC(2,i))*exp(-GN(MODECALC(2,i))*S2(i,:,t)); % THIS IS OK, THE NEW AN AND GN ARE USED ******************************
        
        if PERCALC(2,i) > 1
            PERCALC(2,i) = 1;
        end
        
        % PINF(2,i) = (0.01*abs(randn(1))+1) * ( (1-PERCALC(2,i)) *(1-B(i))^NEIGHCOUNT(i) ); % X" PINF calculation
        PINF(2,i) = (0.01*abs(randn(1))+1)*(1-B(i)/10)*B(i)*(1-PERCALC(2,i)^NEIGHCOUNT(i)); %reviewed with X"
        if PINF(2,i) < 0
            %         disp('PINF 2 < 0!');
            PINF(2,i) = 0; % Just set it to minimum(=0) so that the old value is not retained
        end
        ENERGYSPENTCALC(2,i)=(pktsize*OVERHEAD(MODECALC(2,i))*bitTxcost)*B((i));
        
        % Policy #3 - change mode, change beta
        
        BCALC(3,i) = B(i) + BFACTOR(i)*deltabeta;       % I CHECK BEFORE CHANGING B
        
        if (MODE(i) == 6 && MODEFACTOR(i) > 0) || (MODE(i) == 1 && MODEFACTOR(i) < 0)
            MODECALC(3,i) = MODE(i);
        else MODECALC(3,i) = MODE(i) + MODEFACTOR(i);
        end
        
        PERCALC(3,i) = AN(MODECALC(3,i))*exp(-GN(MODECALC(3,i))*S2(i,:,t)); % THIS IS OK, THE NEW AN AND GN ARE USED ******************************
        
        if PERCALC(3,i) > 1
            PERCALC(3,i) = 1;
        end
        
        % PINF(3,i) = (0.01*abs(randn(1)) + 1)* ( weight*(1-PERCALC(3,i))*(1-BCALC(3,i))^NEIGHCOUNT(i)); % X" --- the factor weight* is arbitrary weight to enhance the choice of this policy
        PINF(3,i) = (0.01*abs(randn(1)) + 1)*(1-BCALC(3,i)/10)*BCALC(3,i)*(1-PERCALC(3,i)^NEIGHCOUNT(i)); %reviewed with X"
        if PINF(3,i) < 0
            %         disp('PINF 3 < 0!');
            PINF(3,i) = 0; % Just set it to minimum(=0) so that the old value is not retained
        end
        ENERGYSPENTCALC(3,i)=(pktsize*OVERHEAD(MODECALC(3,i))*bitTxcost)*BCALC(3,i);
        
        % Policy #4 - keep mode, keep beta
        
        BCALC(4,i) = B(i);
        MODECALC(4,i) = MODE(i);
        AUXPERCALC = mean(PER,1);       % ***PER(i,j)
        PERCALC(4,i) = AUXPERCALC(i);   % ***PER(i,j)
        
        % PINF(4,i) = (0.01*abs(randn(1)) + 1)*((1-PER(i))*(1-B(i))^NEIGHCOUNT(i));  % X"
        PINF(4,i) = (0.01*abs(randn(1)) + 1)*(1-B(i)/10)*B(i)*(1-PER2(i)^NEIGHCOUNT(i)); %reviewed with X"
        if PINF(4,i) < 0
            %         disp('PINF 4 < 0!');
            PINF(4,i) = 0; % Just set it to minimum(=0) so that the old value is not retained
        end
        % I THINK THIS IS NEEDLESS THANKS TO THE USE OF PER2
        % for peer_of_i = 1:N
        %     PINFAUX(4,i) = PINF(4,i) + PINF(4,i,peer_of_i);
        % end % peer_of_i
        % PINF(4,i) = PINFAUX(4,i)/N;
        
        ENERGYSPENTCALC(4,i)=(pktsize*OVERHEAD(MODE(i))*bitTxcost)*B((i));
        
        % End of different options - policies
        
        if sum(MODE(i))>10
            disp('ALARM!!');
        end
        
        %% Utilities Calculation
        for policyindexx = 1:4
            %         METRIC4(policyindexx,i) = PINF(policyindexx,i) - ENERGYSPENT(i)/(maxcostsingle*t)-changes_count(i,max(1,t-1))/(max(1,t-1));
            %         METRIC45(policyindexx,i) = PINF(policyindexx,i) - (ENERGYSPENT(i)+ENERGYSPENTCALC(policyindexx,i))/(maxcostsingle*t)-changes_count(i,max(1,t-1))/(max(1,t-1));
            METRIC7(policyindexx,i) = DUPRATE_NO_TIE(i,t) - ENERGYSPENT(i,t,iteration)/(maxcostsingle*t)-changes_count(i,max(1,t-1))/(max(1,t-1));
            %         METRIC75(policyindexx,i) = DUPRATE_NO_TIE(i,t) - (ENERGYSPENT(i)+ENERGYSPENTCALC(policyindexx,i))/(maxcostsingle*t)-changes_count(i,max(1,t-1))/(max(1,t-1));
            %         METRIC8(policyindexx,i) = NO_TIE(i,t) - (ENERGYSPENT(i)+ENERGYSPENTCALC(policyindexx,i))/(maxcostsingle*t)-changes_count(i,max(1,t-1))/(max(1,t-1));
        end
        
        % Calculate expected Utility over all possible future choices
        UUNCH(i,t) = (action_choice_weight(i,1))*METRIC7(1,i) + (action_choice_weight(i,2))*METRIC7(2,i) + (action_choice_weight(i,3))*METRIC7(3,i) + (action_choice_weight(i,4))*METRIC7(4,i);
        %     UUNCH(i,t) = (action_choice_weight(i,1))*METRIC4(1,i) + (action_choice_weight(i,2))*METRIC4(2,i) + (action_choice_weight(i,3))*METRIC4(3,i) + (action_choice_weight(i,4))*METRIC4(4,i);
        
        
        % I THINK THIS CAN BE SIMPLIFIED AS DONE JUST AFTERWARDS
        %     for jjj = 1:4                  % policy index
        %         MMM(jjj) = METRIC7(jjj,i);  % Vector with policies of the current node
        %     end
        
        %     MMM = METRIC4(:,i);   % Vector with policies of the current node
        MMM = METRIC7(:,i);   % Vector with policies of the current node
        
        metric = max(MMM);      %I find the element with the highest metric FOR THIS NODE
        
        %% Action choice (OST)
        for pindex  = 1:4;              % yes, policy index
            %         if metric == METRIC4(pindex,i); % identifying that I am using the policy with the highest metric
            if metric == METRIC7(pindex,i); % identifying that I am using the policy with the highest metric
                % Start of optimal stopping condition
                
                
                if U(i,t)>UUNCH(i,t)           % utility: check that it is higher than the utility if state does not change = U-UNCHanged % this condition fires OST
                    
                    if t > unknownfactor/exp(1) % OST training time has finished?
                        if  AVEUT(i,t) == max(AVEUT(i,:))    % OST condition
                            (AVEUT(i,t)>max(AVEUT(i,max(1,t-1))));
                            % disp('optimal stopping condition fulfilled!!')
                            % policy adoption start
                            BTEMP(i) = BCALC(pindex,i);                             % policy adoption
                            MODETEMP(i) = MODECALC(pindex,i);                       % policy adoption
                            POLICYPOPUL(i,pindex) = POLICYPOPUL(i,pindex)+1;        % policy popularity update
                            action_choice_weight(i,pindex) = POLICYPOPUL(i,pindex)/t;
                        end  % AVEUT; OST condition
                        % policy adoption end
                        
                    end  % OST training time has finished!
                    
                else   % utility
                    % disp('no benefit'); % do nothing!
                end    % utility
                UTEMP(i,t) = metric;
                U(i,t)=UTEMP(i,t);      % UTILITY CALCULATED HERE!!
                % End of optimal stopping condition
                
                
                
                
                
            end  % if metric...
        end  % policy index
    end % if adaptive
    % End of decision making
    

    
    % Update mode swithcing counters
    
    if MODETEMP(i) > MODE(i)
        UPMODE(i) = UPMODE(i)+1;
    end
    if MODETEMP(i) < MODE(i)
        DOWNMODE(i) = DOWNMODE(i)+1;
    end 
end %RX(i)~=0  %

% Set the sliding window value for the next time slot
% This MUST go before you do MODE = MODETEMP
    if (MODETEMP(i) ~= MODE(i)) || (BTEMP(i) ~= B(i))
        WINDOW(i) = 1;% reset window size at node state change
        changes_count(i,t) = changes_count(i,max(1,t-1)) + 1;
        %                 disp('STATE CHANGE!!');
    else
        WINDOW(i) = WINDOW(i)+1;
        changes_count(i,t) = changes_count(i,max(1,t-1));
    end
end % i    % transmitting node loop
% disp('End of nodes loop');

% Set (beta, mode) at the end of this time slot; update matrices

A = ATEMP;
B = BTEMP;
TX = TXTEMP;
RX = RXTEMP;
ERR = ERRTEMP;
DUP = DUPTEMP;
SUCCESSFUL_INF_COUNT = SUCCESSFUL_INF_COUNT_TEMP;
MODE = MODETEMP;
ERRRATE = ERRRATETEMP;
DUPRATE = DUPRATETEMP;
DUPRATE_NO_TIE_INS(:,t) = (0.01*abs(randn(1))+1)*DUPRATE';
MODEDET(:,t)=MODE;
BDET(:,t)=B;
FRESH = FRESHTEMP;

% the most important ones

Ivstime(iteration,t) = sum(A)/N;
Bvstime(iteration,t) = mean(B);
Modevstime(iteration,t) = mean(MODE);
Chanstatevstime(iteration,t) = mean(CHANSTATET(t,:));  % (iter,t)!!
Freshvstime(iteration,t) = mean(FRESH);


% Other calculations

DUPRATET(t) = mean(DUPRATE);
ERRRATET(t) = mean(ERRRATE);
DUPRATEI(iteration,t) = DUPRATET(t);
ERRRATEI(iteration,t) = ERRRATET(t);

% Moved outside the temporal loop
%
% TXCOUNT(iteration) = sum(TX);
% RXCOUNT(iteration) = sum(RX);
% 
% UPI(iteration) = mean(UP);
% DOWNI(iteration) = mean(DOWN);
% 
% UPMODEI(iteration)=mean(UPMODE);
% DOWNMODEI(iteration)=mean(DOWNMODE);

AUXP = mean(POLICYPOPUL); % auxiliary matrix

for indexx=1:4;
    POLICYPOPUL_ITER(iteration,indexx) = AUXP(indexx);
end

changes_count_global(:,:,iteration) = changes_count;
SUCCESSFUL_INF_COUNT_GLOBAL(iteration,:) = SUCCESSFUL_INF_COUNT;

end % t
%% Iteration-dependent results

TXCOUNT(iteration) = sum(TX);
RXCOUNT(iteration) = sum(RX);

UPI(iteration) = mean(UP);
DOWNI(iteration) = mean(DOWN);

UPMODEI(iteration)=mean(UPMODE);
DOWNMODEI(iteration)=mean(DOWNMODE);


end % iteration
%% Global results

% Energy spent vs time - cumulative

Eavevstime = cumsum( mean(mean(ENERGYSPENTT,3),1) );

% Energy spent vs time - cumulative - normalised

for index_t = 1:T
    ENORM(index_t) = Eavevstime(index_t)/(maxcostsingle*index_t);
end

% plot(Eavevstime);

% Infection rate

Iavevstime = (mean(Ivstime,1));
for retrotime = 1:T
    if Iavevstime == 0.9
        T29C = retrotime
    end
    if Iavevstime == 1
        T2FC = retrotime
    end
end % retrotime

% Beta

Bavevstime = (mean(Bvstime,1));

% Mode

Modeavevstime = (mean(Modevstime,1));

% Information Feshness

Freshavevstime = (mean(Freshvstime,1));
% plot(Freshavevstime);

% channel state

Chanstateavevstime = (mean(Chanstatevstime,1));

% metrics

% efficiency

Efficiency = Iavevstime(T)/T2FC;

% "freshinfect

Freshinfect = zeros(1,T);

for t = 1:T
    Freshinfect(t) = Iavevstime(t)/Freshavevstime(t);
end

plot(Freshavevstime/3.5);
hold on
grid on
plot(Iavevstime,'--');
plot(Freshinfect,'r');
legend('Freshness','Infection','Freshinfect');

% average utility

UEFF=(1+mean(AVEUT,1));  % the factor 1 is arbitrary to keep utility positive for nicer plotting
plot(UEFF);

% infection cost

loglog(Iavevstime,Eavevstime,'r');


% End of results

xx=clock;
hour=num2str(xx(4));
minu=num2str(xx(5));
seco=num2str(xx(6));

% filenamepart = ['cellsim_' aggressistring];
% diary(strcat(filenamepart,'_metric7_',num2str(iterations),'_',date,'_',hour,'_',minu,'.csv'));

filename = strcat('cellsim_',aggressistring,'_metric7_',num2str(iterations),'_',date,'_',hour,'_',minu);
diary(strcat(filename,'.csv'));
diary on;

disp('===Input===');
disp('===Cellsim===');
disp('AMC and adaptive beta - OUR scheme!!');
disp('Nodes');
disp(N);
disp('Temporal duration');
disp(T);
disp('Iterations');
disp(iterations);
disp('aggressiveness');
disp(aggressiveness);
disp('Noise model');
disp('Adaptivity');
disp(adaptive);
disp('Metric');
disp('METRIC X');
disp('Number of originally infected nodes');
disp(A_1);
disp('===Output===');
disp('===Cellsim===');
disp('Infection rate I(T)');
disp(Iavevstime(T));
disp('E2FC');
disp(E2FC);
disp('T2FC');
disp(T2FC);
disp('TC2FC');
disp(TC2FC);
disp('E2FC');
disp(E2FC);
disp('T29C');
disp(T29C);
disp('TC29C');
disp(TC29C);
disp('Final Freshinfect');
disp(Freshinfect(T));
disp('Efficiency');
disp(Efficiency);
disp('=========');
hour,minu,seco
disp('The total simulation duration was...')
xx-yy
disp('===OUR SCHEME w/ METRIC X ===');
disp('===Cellsim===');
disp('== END ==')
diary off;
% filename = ['cellsim_' num2str(iterations) '_' date '_' hour '_' minu '.mat'];
filenamemat = (strcat(filename,'.mat'));
save(filenamemat);