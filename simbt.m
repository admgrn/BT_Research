GEN_TRACE = 1;

BTRatio = 0.33; %Probability of BT transmission
BT_LEN = 625;

if GEN_TRACE == 1

BT_ROUND_SLOT_N = 4; 
BT_MASTER_L = 2;
BT_SLAVE_L  = 0.3;

TotalRoundNumArray = [500,250,166]*4; % with these numbers, we probe for 500*2*625=0.625 seconds
TotalRoundNum = TotalRoundNumArray(BT_ROUND_SLOT_N/2);

BT_SLOT_NUM = BT_ROUND_SLOT_N * TotalRoundNum;
TotalSimTime = BT_SLOT_NUM * BT_LEN;
BTTranRound = rand(1,TotalRoundNum); BTTranRound = ceil(BTTranRound-1+BTRatio);


% for each simu time instant, bluetooth using the channel or not
BT_busy = zeros(1,TotalSimTime);
for btr = 1:TotalRoundNum
    if BTTranRound(btr) == 1
        btrstart = (btr-1)*BT_ROUND_SLOT_N*BT_LEN+1;
        BT_busy(btrstart:btrstart+ceil(BT_MASTER_L*BT_LEN)) = 1;
        tempp = btrstart+(BT_ROUND_SLOT_N-1)*BT_LEN;
        BT_busy(tempp+1:tempp+ceil(BT_SLAVE_L*BT_LEN)) = 1;
    end
end

%BT_busy = zeros(1,TotalSimTime);

RTS_t = 60; 
CTS_t = RTS_t; 
ACK_t = 40; 
SIFS = 10; 
DIFS = 28;

DUMMY_t = 50; %TODO: check this constant

randlossp = 0.2; 
randlossp1good2bad = 0.02;
busylossp = 0.8;

INIT_OFFSET = 300 + 5*BT_LEN;

SampleRes = []; aampleTime = []; SampleDirty = [];

SamplesPlot = zeros(1,TotalSimTime); % for debugging

lastSampleTime = 0; thisSampIdx = 1;
while 1    %if BT_busy(samp0time) + BT_busy(samp0time+RTS_t) == 0
    
    thisSampleTime = ceil(rand(1)*BT_LEN*4) + lastSampleTime;
    
    lastSampleTime = thisSampleTime + RTS_t+CTS_t+2*SIFS+DUMMY_t+SIFS+ACK_t+DIFS;
    if lastSampleTime > TotalSimTime break; end
    
    samp0time = thisSampleTime; 
    samp1time = samp0time + RTS_t+CTS_t+2*SIFS;
    
    SampleTime(thisSampIdx) = samp0time; 
    SampleTime(thisSampIdx+1) = samp1time;
   
    RES0 = 0; RES1 = 0;
    thisdirty = 0;
    
    sample0busy = (BT_busy(samp0time) + BT_busy(samp0time+RTS_t)) > 0;
    sample1busy = (BT_busy(samp1time) + BT_busy(samp1time+DUMMY_t)) > 0;
    
    if sample0busy == 0
        RES0 = rand(1) > randlossp;
    else
        RES0 = rand(1) > busylossp;
    end
    
    if RES0 % will send sample 1 only if sample 0 is good
        if sample0busy == 0 && sample1busy == 0
            tempp = rand(1); 
            if tempp > randlossp1good2bad 
                RES1 = 1;
            else
                RES1 = 0;
                thisdirty = 1; 
            end
        end
        if sample0busy == 0 && sample1busy == 1
            RES1 = rand(1) > busylossp;
        end
        if sample0busy == 1 && sample1busy == 0
            RES1 = rand(1) > randlossp;
        end
        if sample0busy == 1 && sample1busy == 1
            tempp = rand(1); 
            if tempp > randlossp1good2bad 
                RES1 = 1;
            else
                RES1 = 0;
                thisdirty = 1; 
            end
        end
    end
    
    %if thisdirty
    if (RES0 == 1 && RES1 == 0) && ~(sample0busy == 0 && sample1busy == 1)
        XXX = 0;
    end
    
    SampleRes(thisSampIdx) = RES0; 
    SampleRes(thisSampIdx+1) = RES1; 
    SampleDirty(thisSampIdx+1) = thisdirty;

    %for debugging begin
    tempp = 1 + (RES0 == 1 && RES1 == 0) + thisdirty;
    SamplesPlot(samp0time:samp0time+RTS_t-1) = tempp; % for debugging
    SamplesPlot(samp1time:samp1time+DUMMY_t-1) = tempp; % for debugging
    %for debugging end
    
    thisSampIdx = thisSampIdx + 2;
end

SampleTime = SampleTime + INIT_OFFSET;
end %GEN_TRACE

%------------------------------------------------------------------



%------------------------------------------------------------------
% collecting data begin
numbertotalsamples = length(SampleRes)/2;
number10 = 0;
found10time = [];
for j=1:numbertotalsamples
    if SampleRes(j*2-1) == 1 && SampleRes(j*2) == 0
        number10 = number10 + 1;
        found10time(number10) = SampleTime(j*2);
    end
end

cleanT = mod(found10time - CTS_t - 2*SIFS,BT_LEN)+1;
T = BT_LEN; D = DUMMY_t + CTS_t+2*SIFS;

% cleanT(j): the time right after the RTS of a 10 sample. the boundary may start as early as here, i.e., as long as it does not overlap with the RTS to make the 1 in the 10 sample. 
% plus D, is the time the boundary must start before it, i.e., must overlap with the dummy to make the 0 in the 10 sample.

cleanF = zeros(1,T);
for j=1:length(cleanT)
    thisT = cleanT(j);
    if thisT < T-D
        thisbgn = thisT; thisend = thisT+D-1;
        cleanF(thisbgn:thisend) = cleanF(thisbgn:thisend) + 1;
    else
        thisbgn = thisT; thisend = T;
        cleanF(thisbgn:thisend) = cleanF(thisbgn:thisend) + 1;
        thisbgn = 1; thisend = D - (T-thisT) ;
        cleanF(thisbgn:thisend) = cleanF(thisbgn:thisend) + 1;
    end
end

if 1
    %plot(cleanF, 'g');
    %xlabel('Time (microseconds)');
    %ylabel('Number of Collisions');
    %title('Failure of Wi-fi Data Frames in the Presence and Absence of Bluetooth')
    %legend('Sum of Clean & Dirty Frames','Dirty Frames', 'Clean Frames')
end

% collecting data end
%------------------------------------------------------------------

%------------------------------------------------------------------
% detect BT bgn

[maxval, maxloc] = max(cleanF);

BELIEF_RAND10_max_p = 0.01; % NOTE: magic number here, based on the calculation that if there are BT, schedule of 6, 1/3 transmission prob, D = 100, prob of hitting the boundary is about (2/6)*(1/3)*(100/600) = 0.02. 
morethanrandom01 = 0;
sample10ratio = number10 / numbertotalsamples;
if sample10ratio > BELIEF_RAND10_max_p && maxval > mean(cleanF) * 2.5
    % NOTE: will check only if the number of 10 samples is more than normal which will also guarantee enough number of samples. if not, the random nature of the samples may pass the next check
    % NOTE: magic number 2.5 here
    fprintf(1,'10 sample ratio %f, maxval %d, avg %f, We think there is BT!\n', sample10ratio, maxval, mean(cleanF));
else
    fprintf(1,'10 sample ratio %f, maxval %d, avg %f, We think there is NOOOOOOOOOO BT!\n', sample10ratio, maxval, mean(cleanF));
end
% detect BT end
%------------------------------------------------------------------


%------------------------------------------------------------------
% finding boundary bgn

RIGHT_len = CTS_t + 2*SIFS;
shouldtake = zeros(1,length(found10time));

estmidoint = maxloc;
rangemargin = 10; % NOTE: magic number here
for j=1:length(found10time)
    dummytime = mod(found10time(j), BT_LEN);
    shouldtakethisone = 0;  
    if dummytime < estmidoint
        if estmidoint - dummytime <= DUMMY_t + rangemargin 
            shouldtakethisone = 1;
        end
    else
        if estmidoint + BT_LEN- dummytime <= DUMMY_t + rangemargin 
            shouldtakethisone = 1;
        end
    end    
    if dummytime > estmidoint
        if dummytime - estmidoint <= RIGHT_len + rangemargin 
            shouldtakethisone = 1;
        end
    else 
        if dummytime + BT_LEN- estmidoint <= RIGHT_len + rangemargin 
            shouldtakethisone = 1;
        end
    end
    shouldtake(j) = shouldtakethisone;
end
totaltakenum = sum(shouldtake);
measureddirtyp = (length(found10time) - totaltakenum) / (BT_LEN-D); % probability that a random us landed a dirty sample

EXPECTED_R_LEN = 6; % NOTE: if we send 1000 probes, the overlapping is expected to be 6 us
target_p = 0.005; % NOTE: prob that the output does not contain the actual boundary
one_p = EXPECTED_R_LEN * measureddirtyp;
for stepdown=1:4 %NOTE: max go down by 4
    if power(one_p,stepdown) < target_p 
        break;
    end
end

increment_array = [-1,1];
    
outputrange = zeros(1,2);
for j=1:2
    thisinc = increment_array(j);
    currpnt = maxloc;  currval = maxval;
    currdist = 0;
    while 1
        currpnt = currpnt + thisinc; 
        currdist = currdist + 1;
        if currpnt == 0
            currpnt = BT_LEN;
        end
        if currpnt == BT_LEN + 1
            currpnt = 1;
        end
        currval = cleanF(currpnt);
        % NOTE: to apply this, has to assume that there is at most one dirty sample landed in the actual intersection. So, has to have enough number of clean samples.
        if currval < (maxval- stepdown) 
            fprintf(1,'broke from the loop when the val is %d, measureddirtyp %f, stepdown %d\n', currval,measureddirtyp, stepdown);
            break;
        end
    end
    outputrange(j) = currpnt;
end
 
outputrange

% finding boundary end
%------------------------------------------------------------------

%------------------------------------------------------------------
% finding schedule begin, include length (6,4,2) and the beginning slot

% TODO: handling the case when the actual boudary is very close to our
% preset boundary

MAX_SLOT_NUM = 6;
scores = zeros(2,MAX_SLOT_NUM); % scores(1,1:4), scores for 4, scores(2,1:6), for 6
slenarray = [4,MAX_SLOT_NUM];
for j=1:2
    slen = slenarray(j);
    for jj=1:slen
        thisestbgn = estmidoint + (jj-1) * BT_LEN;
        for jjj=1:length(found10time)
            if shouldtake(jjj) == 1 
                thisslotidx = found10time(jjj) - thisestbgn;
                thisslotidx = mod(thisslotidx,(slen*BT_LEN));
                thisslotidx = round(thisslotidx/BT_LEN);
                if thisslotidx ~=0 && thisslotidx ~= slen - 1 && thisslotidx ~= slen 
                    scores(j,jj) = scores(j,jj)+1;
                end
            end
        end
    end
    for jj=slen+1:MAX_SLOT_NUM
        scores(j,jj) = 1000000;
    end
end

scores

minest4 = min(scores(1,:)); 
minest6 = min(scores(2,:)); 
MAGIC_C = 0.3; % NOTE: if 2, there will be equal number of 10 samples in each slot. if the actual schedule is 6, 2/3 of them will reegister a point. if 4, 1/2. Given there are random variations, 0.3 turns out to be a good number

if minest4 > totaltakenum * MAGIC_C && minest6 > totaltakenum * MAGIC_C
    fprintf(1,'we think the schdule is 2!!!\n');
else
    if minest4 < 0.5 * minest6 % NOTE: 2 because a dirty sample is twice more likely to register a point with 6 than with 4
        fprintf(1,'we think the schdule is 4!!!\n');
    else
        fprintf(1,'we think the schdule is 6!!!\n');
    end
end

% finding schedule end
%------------------------------------------------------------------


fid = fopen('output.txt','w');

fprintf(fid,'%d\n',numbertotalsamples);
for i=1:length(found10time)
    fprintf(fid,'%d\n',found10time(i));
end

fclose(fid);
