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

%% Initial cleanup
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
numberofbetas = 25;
deltabeta = 0.02;

for inindex=1:numberofbetas
    BETAS(inindex)=0.20+inindex*deltabeta;
end
BETAS;

MODES = [1 2 3 4 5 6];
STATES = [1 2 3 4 5 6];

N = 50; % number of nodes
% b0 = BETAS(10); % original beta
T = 100; % number of temporal iterations
I = 1; % number of originally infected nodes
loss = 0.8; % this is essentially BER (attention, not PER!!)
signal = 30e+015;
n0 = 15;
sigma = 7;
deltaoverbeta = 0.1;
netdens = 0.85;
% delta = beta/5
betasum = 0;
pktsize = 50; % packet size in bits
bitTxcost = 720; %bit transmission energy cost in nJ
bitRxcost = 110; %bit reception energy cost in nJ
instructioncost = 4;
nb = 45; %amplitude of noise burst
trainprob = 0.2; %probability that a noise burst occurs
iterations = 5; %algorithm's iterations
startingmode = 3;
resetsize=7;
weight = 1.0;
startingchannelstate = 3;
actualtime = 0;
maxcostsingle = 1.8*pktsize*bitTxcost; % this is max cost for a SINGLE transmission
berrors = 0;
bergood = 0;
   
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
for dummynindex1 = 1:7
        SNRLIMITS(dummynindex1)=GPN(dummynindex1);
end
for dummynindex2 = 1:7
    NOISELIMITS(dummynindex2)=signal/SNRLIMITS(dummynindex2);
end

%% Bursty noise
% Flapping function; it's a square function
 ft = 3; % flapping duration
 for f = 1:T
     if mod(f,25) == 0
         for ff = f:f+ft+5  % guard time
             FLAP(ff) = 2;
         end
     else
         % do nothing
     end
 end
%% Import from SoS model
open transMatrix_nsr_08.mat;
TM = ans.TM_ave;
clear ans;
%% Iteration loop
for iteration = 1:iterations
    iteration
%% Iteration-dependent initialization

TC2FC(iteration)=0;
T2FC(iteration)=0;
E2FC(iteration)=0;
    
TC29C(iteration)=0;
T29C(iteration)=0;
E29C(iteration)=0;

for policyindex=1:4;
    POLICYCHOSENI(iteration,policyindex)=0;
end

% action_choice_weight(node,chosen_policy)
% renamed from SSS(node,chosen_policy)
action_choice_weight = zeros(N,4)+0.25;
changes_count = zeros(N,T);

% initialization of noise for this interation; noise is NOISET(node,othernode,time)
for iindexxx = 1:N;
    for jjindexxx = 1:N;
    randostate(iindexxx,jjindexxx) = rand;
    if randostate(iindexxx,jjindexxx) <= 0.0835 % so that (length(STATES)*randostate(node)) >= 0.5
        randostate(iindexxx,jjindexxx) = 0.084; % int16(.) will safely give 1
    end
    end % jjndexxx
end % iindexxx


NOISE = zeros(N,N) + n0;
SNR = zeros(N,N,T);
NOISET = zeros(N,N,T) + n0;
NOISEPURET = zeros(N,N,T) + n0;
RECENTNOISE = zeros(N,N,T) + n0;
AVENOISET = zeros(N,N,T) + n0;
S2 = mean(SNR,2);

for iindexx = 1:N;
    for jjindexx = 1:N;
        CHANSTATE(iindexx,jjindexx) = STATES(int16(length(STATES)*randostate(iindexx,jjindexx))); % random original channel state; initialize matrix with channel states (->SNR) as perceived by EACH NODE!!      
    end % jjindexx
end % iindexx

CHANSTATETEMP = CHANSTATE;
    
DUPT=zeros(T,iterations);

A = zeros(N);
ATEMP = zeros(N);
A(1) = 1;
ATEMP(1) = 1;


for node10 = 1:N
    randobeta(node10) = rand;
    if randobeta(node10) <= 0.02 % (length(BETAS)*randobeta(node))
        randobeta(node10) = 0.022;
    end
end % node10 ...is an index only

for node11 = 1:N
    randomode(node11) = rand;
    if randomode(node11) <= 0.0834 % so that (length(MODES)*randomode(node)) >= 0.5
        randomode(node11) = 0.0834;
    end
end % node11 ...is an index only
TX = zeros(N);
RX = zeros(N);
TXTEMP = zeros(N);
RXTEMP = zeros(N);
ERR = zeros(N);
ERRTEMP = zeros(N);
DUP = zeros(N);
DUPTEMP = zeros(N);
BTEMP = zeros(N);
nodenoise = zeros(N);
avenoise = zeros(N);
WINDOW = zeros(N) + 1;
% B = zeros(N) + b0; 
BFACTOR = zeros(N);
ENERGYSPENT= zeros(N);
ENERGYSPENTTX= zeros(N);
ENERGYSPENTRX= zeros(N);
OVERHEADBITSENT= zeros(N);
UP = zeros(N); % number of times the beta has been forced up when minimum SNR-limit reached - ACCORDION-RELATED
DOWN = zeros(N); % number of times the beta has been forced down when maximum SNR-limit reached - ACCORDION-RELATED
DOWNMODE = zeros(N);
UPMODE = zeros(N);
DUPRATE = zeros(N);
ERRRATE = zeros(N);
ERRRATETEMP = zeros(N);
DUPRATETEMP = zeros(N);
% PER = zeros(N);
BITSENT= zeros(N); % Number of bits sent by a node. So that when I sum up for all nodes I get the total bits sent.
%BER = zeros(N) + loss;

for node = 1:N
    B(node) = BETAS(int16(length(BETAS)*randobeta(node))); % random original beta
    MODE(node)= MODES(int16(length(MODES)*randomode(node))); % random original mode
end

MODETEMP = MODE; % arbitrary; value plays no role

CHANSTATET = zeros(T,N,N) + startingchannelstate;
UPIFSMC = zeros(N,N); % state increment count for channel status (input FSMC)
DOWNIFSMC = zeros(N,N); % state increment count for channel status (input FSMC)
SWITCHIFSMC = zeros(N,N) + UPIFSMC + DOWNIFSMC;
UPIFSMCTEMP = zeros(N,N);
DOWNIFSMCTEMP = zeros(N,N);
PER = zeros(N,N);
NOISEPURE = zeros(N,N) + n0;
PERCEIVEDADJ = zeros(N,N); % PERCEIVED adjucency matrix! = KNOWN neighbors
POLICYCHOSEN = zeros(N,4); % there are 4 policies-possible actio

for t=1:T;
    nodenoise(t) = 0; % dummy initialization of nodenoise
    g(t)=0.5*randn(1);
end

% NOISET = zeros(N,T) + n0;
% NOISEPURET = zeros(N,T) + n0;
% RECENTNOISE = zeros(N,T) + n0;
% AVENOISET = zeros(N,T) + n0;
RECENTUT = zeros(N,T);
AVEUT = zeros(N,T);
U = zeros(N,T);
    
MODEDET = zeros(N,T) + startingmode; % MODEDET(= mode detailed) like MODE(node) but keeps temporal evolution too
% CHANSTATET = zeros(N,T) + startingchannelstate;
% BDET = zeros(N,T) + b0; % BDET(= beta detailed) like B(node) but keeps temporal evolution too

for node = 1:N
    for t = 1:T
        BDET(t,node) = B(node); % BDET(= beta detailed) like B(node) but keeps temporal evolution too
    end
end

% Decision matrices initialization
% These values are indifferent; just for the matrices to be non-void
% These matrices are always changed; every node, timeslot & iteration

MODECALC = zeros(4,4) + 3;
PERCALC = zeros(4,4);

for a = 1:N;        % number of nodes
    for b = 1:4;    % number of possible actions
        BCALC(b,a) = B(a);
    end
end
%% Adjacency matrix
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

%% Markov Chain refresh - moved out of the big temporal loop
% Refresh channel state markov chain for all LINKS (NOT nodes) start    
for t = 1:T
    for node20 = 1:N
        for node21 = 1:N
            a = 1:6;             %# possible channel states - modes
            weight = TM(CHANSTATE(node20,node21),:); %# corresponding weights
            results = 1;              %# how many numbers to generate
            R = a( sum( bsxfun(@ge, rand(results,1), cumsum(weight./sum(weight))), 2) + 1 );
            CHANSTATETEMP(node20,node21) = R; % make it CHANSTATETEMP(node20,node21,time) and move it out of the temporal loop
            if CHANSTATETEMP(node20,node21) > CHANSTATE(node20,node21)
                UPIFSMCTEMP(node20,node21) = UPIFSMCTEMP(node20,node21) + 1;
            end
            if CHANSTATETEMP(node20,node21) < CHANSTATE(node20,node21)
                DOWNIFSMCTEMP(node20,node21) = DOWNIFSMCTEMP(node20,node21) + 1;
            end
        end  % node21
    end % node20
    
    CHANSTATE = CHANSTATETEMP;
    DOWNIFSMC = DOWNIFSMCTEMP;
    UPIFSMC = UPIFSMCTEMP;
    % SWITCHIFSMC = DOWNIFSMC + UPIFSMC;
    CHANSTATET(t,:,:) = CHANSTATE;
    % Refresh channel state markov chain for all LINKS (NOT nodes) end
end

%% Attempt to infect banner

% disp('===NOW ATTEMPT TO INFECT===')

%% Temporal loop - Epidemic dissemination loop
    for t = 1:T
        actualtime = actualtime + 0.001;
% %% Markov Chain refresh - Time slot-dependent initialization
% % Refresh channel state markov chain for all LINKS (NOT nodes) start    
% 
% for node20 = 1:N
%     for node21 = 1:N
%         a = 1:6;             %# possible channel states - modes
%         weight = TM(CHANSTATE(node20,node21),:); %# corresponding weights
%         results = 1;              %# how many numbers to generate
%         R = a( sum( bsxfun(@ge, rand(results,1), cumsum(weight./sum(weight))), 2) + 1 );
%         CHANSTATETEMP(node20,node21) = R; % make it CHANSTATETEMP(node20,node21,time) and move it out of the temporal loop
%         if CHANSTATETEMP(node20,node21) > CHANSTATE(node20,node21)
%             UPIFSMCTEMP(node20,node21) = UPIFSMCTEMP(node20,node21) + 1;
%         end
%         if CHANSTATETEMP(node20,node21) < CHANSTATE(node20,node21)
%             DOWNIFSMCTEMP(node20,node21) = DOWNIFSMCTEMP(node20,node21) + 1;
%         end
%     end  % node21
% end % node20
% 
% CHANSTATE = CHANSTATETEMP;
% DOWNIFSMC = DOWNIFSMCTEMP;
% UPIFSMC = UPIFSMCTEMP;
% % SWITCHIFSMC = DOWNIFSMC + UPIFSMC;
% CHANSTATET(t,:,:) = CHANSTATE;
% % Refresh channel state markov chain for all LINKS (NOT nodes) end

%% Transmitting node-related initializations

for node = 1:N
    
    for othernode = 1: N % othernode
        if CHANSTATE(node,othernode)~=7
            NOISEPURE(node,othernode) = (  NOISELIMITS(CHANSTATE(node,othernode)) + ( NOISELIMITS(CHANSTATE(node,othernode)+1)-NOISELIMITS(CHANSTATE(node,othernode))   )*randn(1,1)  ); % noise randmoly between limits
            NOISE(node,othernode) = (1+FLAP(t)) * NOISEPURE(node,othernode); % noise randmoly between limits
            %                     NOISE(node,othernode) = NOISEPURE(node,othernode); % no BURSTY noise
        else
            if CHANSTATE(node,othernode)==7
                NOISE(node,othernode) = NOISELIMITS(7)*(1+rand);
            end
        end
        
        NOISET(node,othernode,t)=NOISE(node,othernode);
        NOISEPURET(node,othernode,t)=NOISEPURE(node,othernode);
        SNR(node,othernode,t)= signal/(NOISET(node,othernode,t)); % Populate signal-to-noise matrix: a value for each node and timeSLOT.
        S2=mean(SNR,2);
        M = 2^(MODE(node));
        BER(node,othernode) = 0.2*exp(-3*SNR(node,othernode,t)/(2*(M-1)));
        
    end % othernode
    
    % use Giannaki for PER calculation:
    for othernode=1:N
        if SNR(node,othernode,t)>GPN(MODE(node))
            PER(node,othernode) = AN(MODE(node))*exp(-GN(MODE(node))*SNR(node,othernode,t));
        elseif (SNR(node,othernode,t)<GPN(MODE(node)) && SNR(node,othernode,t)>0)
            PER(node,othernode) = 1;
        end
        
        if PER(node,othernode)>1
            %disp('***PER ALARM !!!***');
            %disp(PER(node));
            berrors = berrors +1;
            PER(node,othernode) = 1;
        else bergood = bergood+1;
        end % if
    end % othernode
end

%% Transmitting node loop

for i = 1:N
    %attempt to cure
    cure = rand;
    if cure < B(i)*deltaoverbeta;
        A(i) = 0;
    end
    lucky = rand;
    if lucky < B(i) % here we make the beta experiment ONCE for each node: broadcast!
        %% Receiving node loop
        for j = 1:N
            if A(i) == 1
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
                    ENERGYSPENT(j)  =ENERGYSPENT(j)  +pktsize*OVERHEAD(MODE(i))*bitRxcost+2*instructioncost; %packet size decided by sender's code rate!
                    ENERGYSPENTRX(j)=ENERGYSPENTRX(j)+pktsize*OVERHEAD(MODE(i))*bitRxcost+2*instructioncost; %packet size decided by sender's code rate!
                    channel = rand;
                    if channel > PER(i,j) %loss - PER(i,j) NOW! i=src, j=dst   ***PER(i,j)
                        if A(j)==0
                            ATEMP(j)=1;
                        else
                            DUPTEMP(j)=DUP(j)+1;
                        end
                    else
                        ERRTEMP(j)=ERR(j)+1;
                    end % loss
                    % 	                   end % lucky % here we make the beta experiment on a p2p basis!
                end % ADJ
            end % A(i)
        end % j
        
        TXTEMP(i)=TX(i)+1;
        ENERGYSPENT(i)  =ENERGYSPENT(i)  +pktsize*OVERHEAD(MODE(i))*bitTxcost; % energy calculated ONCE for every broadcast!!
        ENERGYSPENTTX(i)=ENERGYSPENTTX(i)+pktsize*OVERHEAD(MODE(i))*bitTxcost;
        BITSENT(i) = BITSENT(i) + pktsize + pktsize*OVERHEAD(MODE(i));
        OVERHEADBITSENT(i)=OVERHEADBITSENT(i)+pktsize*OVERHEAD(MODE(i));
    end % lucky % here we make the beta experiment ONCE for each node: broadcast!
    
    %Calculate next beta & encoding mode
    if RX(i)~=0
        %BTEMP(i)= b0*(1+ERR(i)/RX(i)-DUP(i)/RX(i));
        ERRRATE(i)=ERR(i)/RX(i);
        DUPRATE(i)=DUP(i)/RX(i);
        
        ERRRATETEMP(i)=ERRTEMP(i)/RXTEMP(i);
        DUPRATETEMP(i)=DUPTEMP(i)/RXTEMP(i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % set MODE according to SNR limits...
        % and change beta only when SNR-thresholds are crossed
        
        % No Giannakis mode choosing -------------------------
        
        % No accordion action -------------------------------
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              

%% Utility calculation

%% Action choice (OST)


    end %RX(i)
    if RX(i)==0
        %BTEMP(i) = b0; %whichever you like most. the second is more 'natural'
        BTEMP(i) = B(i);
    end
    
end  % i % end of transmitting node loop
%% Time slot-dependent results
    end % t
%% Iteration-dependent results
end % iteration
%% Global results

xx=clock;
hour=num2str(xx(4));
minu=num2str(xx(5));
seco=num2str(xx(6));
diary(strcat('cellsim_',date,'_',hour,'_',minu,'.csv'));
diary on;

disp('===Input===');
disp('AMC and adaptive beta - OUR scheme!!');

disp('===Output===');

disp('=========');
disp('===OUR SCHEME W/ METRIC METRIC7===');
disp('== END ==')
save cellsim.mat
hour,minu,seco
disp('The total simulation duration was...')
xx-yy