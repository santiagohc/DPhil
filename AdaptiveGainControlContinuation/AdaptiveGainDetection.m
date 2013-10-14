function[] = AdaptiveGainDetection(NumTrialsTot, AllMeans, AllConcentrations, AllContrasts, AllBoundaryAngles)
%[data] = AdaptiveGainSensitivity()
%
%...
tic
NameFileToSave = 'Prueba';

%% Define Initial variables if empty

%
BlueCircleRight = 0; %1 if you want the blue circle response cue to appear on the right side

% a DEFAULT set of MEAN deviations from the test value for each STREAM of stimuli
if isempty(AllMeans)
    AllMeans = [-20 -10 10 20];
end

%a DEFAULT set of CONCENTRATIONS for the stream of stimuli...
if isempty(AllConcentrations)
    AllConcentrations = [5];
end

%a DEFAULT set of CONTRAST  values for the stimuli...
if isempty(AllContrasts)
    AllContrasts = [30];
end

%a DEFAULT number of TRIALS per BLOCK
if isempty(NumTrialsTot)
    NumTrials = 30;%numel(AllMeans)*numel(AllVariances)*numel(AllContrasts);
else
    NumTrials = NumTrialsTot;
end

%a DEFAULT number of TRIALS per BLOCK
if isempty(AllBoundaryAngles)
    AllBoundaryAngles = 1:360;%[90];
end

NumCombinations = numel(AllMeans)*numel(AllConcentrations)*numel(AllContrasts);

% Compute the optimal amount of trials if they don't surpass 300 else make it 300
% if isempty(NumTrialsTot) && NumCombinations < 300 && NumCombinations > 20
%     NumTrials = NumCombinations;
% elseif ~isempty(NumTrialsTot) && NumCombinations > 300
%     NumTrials = 300;
% elseif ~isempty(NumTrialsTot) && NumCombinations < 20
%     NumTrials = 20;
% end




%% Define other initial parameters

%Define Number of samples(stimuli) that each stream will have
NumMaxSamplesStream = 10;

%Define Number of samples(stimuli) that each stream will have
NumMinSamplesStream = 1;

%Define the initial total reward to zero
TotalReward = 0;
TotalRewardLimited = 0;

%Duration of each mask and sample
DurMask = 0.15;
InterMaskDur = 0.05; % sampITI
TimeStartSamp = 100;
TimeEndSamp = 250;
InterSampInterval = TimeStartSamp/1000;
SampDur = (TimeEndSamp - TimeStartSamp)/1000;
data.sampdur = SampDur;
data.intersampinterval = InterSampInterval;

%Duration of Gap after Stimuli presentation and before Response Cue presentation
DurGap = .2 + rand(1,NumTrials)*.3;

maxMuDevCrit = .5;%Critical maximum acceptable value of MEAN in stream
maxKappaDevCrit = .5;%Critical maximum acceptable value of VARIANCE in stream

data.maxmudevcrit = maxMuDevCrit;
data.maxkappadevcrit = maxKappaDevCrit;

% Pseudorandom string of 1's and 2's that says whether there id 1 or 2 masks
NumMasksIni = ceil(rand(1,NumTrials)*2);

% Pseudorandom string of 1's and 2's that says whether there id 1 or 2 masks
NumMasksEnd = ceil(rand(1,NumTrials)*2);

%Define whether it is an Integration-of-Mean trial or a
%Sensitivity-to-Angle trial Currently set to 2:1 proportion
TypeOfTrial = (Shuffle(rem((1:(ceil(NumTrials/3))*3),3))+1)<3;
TypeOfTrial =  TypeOfTrial(1:NumTrials); % zeros(NumTrials,1);%

%Random vector of Sample-Numbers, we need to take the those values only for
%the Sensitivity-to-Angle trials
IsTypeTwoWhichSampleSize = Shuffle(rem((1:(ceil(sum(~TypeOfTrial)/(NumMaxSamplesStream-NumMinSamplesStream)))*NumMaxSamplesStream-NumMinSamplesStream+1),NumMaxSamplesStream-NumMinSamplesStream+1))+NumMinSamplesStream;%(ceil(rand(1,NumTrials)*(NumMaxSamplesStream-NumMinSamplesStream+1))) + NumMinSamplesStream-1;
NumSamplesOddTrial = IsTypeTwoWhichSampleSize(1:sum(~TypeOfTrial));
NumSamplesEachTrial = NumMaxSamplesStream*ones(1,NumTrials);
NumSamplesEachTrial(~TypeOfTrial) = NumSamplesOddTrial;

%Random Vector of the Desired Means for each trial
DesiredMeansEachTrial = Shuffle(rem((1:(ceil(NumTrials/numel(AllMeans)))*numel(AllMeans)),numel(AllMeans)))+1;%Shuffle(repmat(1:numel(AllMeans),1,ceil(NumTrials/numel(AllMeans))));
DesiredMeansEachTrial = AllMeans(DesiredMeansEachTrial(1:NumTrials));

% InitialDificultyFactor = 20;
% SignMeanEachTrial = -2 * (DesiredMeansEachTrial < 0) + 1;

%Random Vector of the Desired Concentrations for each trial
DesiredConcentrationEachTrial = Shuffle(rem((1:(ceil(NumTrials/numel(AllConcentrations)))*numel(AllConcentrations)),numel(AllConcentrations)))+1;% level of across-element CONCENTRATION in the stream
DesiredConcentrationEachTrial = AllConcentrations(DesiredConcentrationEachTrial(1:NumTrials));

%Random Vector of the Desired Contrasts for each trial
DesiredContrastsEachTrial = Shuffle(rem((1:(ceil(NumTrials/numel(AllContrasts)))*numel(AllContrasts)),numel(AllConcentrations)))+1;% level of across-element CONCENTRATION in the stream
DesiredContrastsEachTrial = AllContrasts(DesiredContrastsEachTrial(1:NumTrials));

% Values that the angle at which the pink-blue boundary lies can take.
TrialBoundaryAngles = Shuffle(rem((1:(ceil(NumTrials/numel(AllBoundaryAngles)))*numel(AllBoundaryAngles)),numel(AllBoundaryAngles)))+1;
TrialBoundaryAngles = AllBoundaryAngles(TrialBoundaryAngles(1:NumTrials));
% Define the reference angle for each trial
% TrialBoundaryAngles = ceil(361*rand(NumTrials,1))-1;% One degree resolution

%Random Vector of the Position of the Blue Response Probe for each trial
AppearBlueCircleRight = Shuffle(rem((1:(ceil(NumTrials/2))*2),2))+1;% One means right
AppearBlueCircleRight = BlueCircleRight*ones(NumTrials,1);% Fixates the Blue Response Probe to the right

%Random Vector of the Presence or absence of a Gabor for detection
AppearGaborInNoise = Shuffle(rem((1:(ceil(NumTrials/2))*2),2));% One means Present
ContrastDetection = 30;
NoiseFactorDetection = 128;

%% Open Screen Window

[w, rect] = Screen('OpenWindow', 0, 255,[],32,2);%[0 0 900 600],32,2);%[10 10 800 600],32,2);% open SCREEN of PSYCHTOOLBOX with certain dimensions
HideCursor

%% Define some parameters

DiamCircle = 300;%Radius of th CENTRAL CIRCLE (try making it an even number)
RadiusCircle = DiamCircle/2;
ExtraWidthPercentage = 0.1;
SizeCue = DiamCircle * (1 + ExtraWidthPercentage);
SizeCue = floor(SizeCue/2)*2;
FixationPointSize = 5;
GaborFactor = 7;
fixColour = 100;
GaborFreq = 0.02;


ScreenSize = rect(3:4);%width and height of SCREEN
CentreOfScreen = [ScreenSize(1)/2 ScreenSize(2)/2];% rect([3 4]) are the width and height of the screen.
CentralRectangle = [CentreOfScreen CentreOfScreen] + [-DiamCircle/2 -DiamCircle/2 (DiamCircle/2)-1 (DiamCircle/2)-1];%Rectangle where STIMULUS CIRCLE will be inscribed
RewardBox = floor(CentralRectangle + [-DiamCircle*.05 0 DiamCircle*.05 0]) + [0 (DiamCircle + floor(.40*DiamCircle)) 0 (floor(.50*DiamCircle))];%Reward box will start at a distance of 15% of the screen from the end of the Central Circle
RewardBoxFilling = RewardBox; RewardBoxFilling(3)= RewardBoxFilling(1);%Give height and and position but not width to the filling

CentralCircle = CreateCircleINT(DiamCircle,3)*-128;%CreateCircle(size,width,falloff)
patch=CreateCircularApertureINT(DiamCircle);
IsCentralCircle = patch >.80;
CueRectangle = [CentreOfScreen CentreOfScreen] + [-SizeCue/2 -SizeCue/2 (SizeCue/2)-1 (SizeCue/2)-1];%Rectangle where STIMULUS CIRCLE will be inscribed
HalfDifCueCircle = (SizeCue - DiamCircle)/2;
ResponseRectRight = [CentreOfScreen CentreOfScreen] + [+DiamCircle/2 -DiamCircle/4 (DiamCircle) (DiamCircle/4)];%Rectangle where STIMULUS CIRCLE will be inscribed]
ResponseRectLeft  = [CentreOfScreen CentreOfScreen] + [-(DiamCircle) -DiamCircle/4 -(DiamCircle/2) (DiamCircle/4)];%Rectangle where STIMULUS CIRCLE will be inscribed
ResponseRectRight = ResponseRectRight + [DiamCircle*.1 0 DiamCircle*.1 0];
ResponseRectLeft  = ResponseRectLeft  - [DiamCircle*.1 0 DiamCircle*.1 0];
CentreCueRectangle = false(SizeCue,SizeCue);
CentreCueRectangle(1+HalfDifCueCircle: end-HalfDifCueCircle,1+HalfDifCueCircle: end-HalfDifCueCircle) = true;
CentreCueCircle = false(SizeCue,SizeCue);
CentreCueCircle(CentreCueRectangle) = IsCentralCircle;
CentralCircleResp = CentralCircle;
ProgressBar = RewardBox([1 2 1 4]);
ProgressBarWhole = RewardBox;
ProgressBarLength = RewardBox(3)-RewardBox(1);
[x,y] = meshgrid((1:DiamCircle)-(DiamCircle+1)/2);
FinalCueCompRGB = nan(DiamCircle, DiamCircle,3);
TextRectangle1 = ProgressBarWhole + [0  CentreOfScreen(2)  0 CentreOfScreen(2) ];
TextRectangle1Height = abs(TextRectangle1(2) - TextRectangle1(4));
TextRectangle2 = TextRectangle1 + [0  CentreOfScreen(2)+ 2*TextRectangle1Height  0 CentreOfScreen(2)+ 2*TextRectangle1Height ];
PresentTextPosition = [ CentreOfScreen(2)+ 2*TextRectangle1Height];
AbsentTextPosition = 0;

FixationPoint = [CentreOfScreen CentreOfScreen] + [-FixationPointSize -FixationPointSize FixationPointSize FixationPointSize];




% Create a progress bar
% h=waitbar(0,'building sequence...');
% set(findobj(h,'type','patch'), ...
%     'edgecolor','b','facecolor','b')





try %Welcome Screen until click
    
    %     SizeWelcomeFont = floor(ScreenSize(2)*.05);
    
    if NumTrials < 108;
        Screen(w,'FillRect',255);  % black screen
        Screen('TextSize', w , floor(ScreenSize(2)*.05) );% Set the font size
        Screen('TextFont', w, 'Calibri');% Set the font of text
        Screen(w,'DrawText','PRACTICE',(CentreOfScreen(1)-50),(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
        Screen(w,'DrawText','Press any button when ready...',(CentreOfScreen(1)-200),(CentreOfScreen(2))+100,0);%Initial Instructions
    else
        Screen(w,'FillRect',255);  % black screen
        Screen('TextSize', w , floor(ScreenSize(2)*.05) );% Set the font size
        Screen('TextFont', w, 'Calibri');% Set the font of text
        Screen(w,'DrawText','ACTUAL EXPERIMENT',(CentreOfScreen(1)-150),(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
        Screen(w,'DrawText','Press any button when ready...',(CentreOfScreen(1)-250),(CentreOfScreen(2))+100,0);%Initial Instructions
        
    end
    Screen(w,'Flip');%Print to SCREEN
    
    WaitPressButton = 0;
    %Wait for Subject to press a Button to initiate the task...
    while WaitPressButton==0;
        [kdown secs keyCode] = KbCheck;  % check for key press
        if kdown == 1
            WaitPressButton = 1;
        end
    end
    
    %Have a black space after KeyPress before the task begins
    Screen(w,'FillRect',128);  % Gray screen
    Screen(w,'Flip');        % write to screen
    WaitSecs(1);
    
    
catch
    Screen('closeall');% Close Screen
    ShowCursor;%
    rethrow(lasterror);%InformError
    %     save datatmp data;%Save data acquired so far
    
end




%% Define variables to be filled for each trial and sample
% Keep Record of Trial Duration
TrialsLengthTime = nan(NumTrials,1);
ResponseRightSide = nan(NumTrials,1);
CorrectResponse = nan(NumTrials,1);

AllAngles = nan(NumMaxSamplesStream,NumTrials);
AllAnglesDegreesObserved = nan(NumMaxSamplesStream,NumTrials);
RelativeMeanEachTrial = nan(1,NumTrials);
ObservedConcentrationEachTrial = nan(1,NumTrials);
DesiredMeansCorrectedEachTrial = DesiredMeansEachTrial;
DifficultyFactorTrial10 = nan(1,NumTrials);
TotalRewardEachTrial = nan(1,NumTrials);


NoGaborRGBAllTrials = cell(1,1);
GaborRGBAllTrials = cell(NumMaxSamplesStream,1);

IsBlueAngle = false(NumTrials,1);
ErrorAngle = nan(1,NumTrials);
TrueAngle = nan(1,NumTrials);
ReportedAngle = nan(1,NumTrials);
NTrCompDiff = 5;

try
    
    mask=SCreateGaborINT(DiamCircle,GaborFactor*DiamCircle,0,GaborFreq,0,30)+ SCreateGaborINT(DiamCircle,GaborFactor*DiamCircle,90,GaborFreq,0,30);
    DifficultyFactor10 = 0;
    DifficultyFactor20 = 0;
    for Trial = 1: NumTrials
        
        if rem(Trial,100) == 1 && Trial > 1
            TotalRewardLimited = 50;
            BlockNumber = ceil(Trial/100);
            Screen(w,'FillRect',255);  % black screen
            Screen(w,'DrawText',strcat('Block', num2str(BlockNumber)-1, ' ', 'Completed'),(CentreOfScreen(1)-150),(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
            Screen('TextFont', w, 'Calibri');% Set the font of text
            Screen('TextSize', w , floor(ScreenSize(2)*.05) );% Set the font size
            Screen(w,'DrawText','Press any button when ready to continue',(CentreOfScreen(1)-250),(CentreOfScreen(2))+100,0);%Initial Instructions
            Screen(w,'Flip');        % write to screen
            pause
        end
        
        if rem(Trial,100) == 1
            BlockNumber = ceil(Trial/100);
            Screen(w,'FillRect',128);  % black screen
            Screen(w,'DrawText',strcat('Block', num2str(BlockNumber)),(CentreOfScreen(1)-50),(CentreOfScreen(2))-100,0);%Inform Subject of PRACTICE session
            Screen('TextFont', w, 'Calibri');% Set the font of text
            Screen('TextSize', w , floor(ScreenSize(2)*.05) );% Set the font size
            Screen(w,'DrawText','Get Ready...',(CentreOfScreen(1)-70),(CentreOfScreen(2)),0);%Initial Instructions
            Screen(w,'Flip');        % write to screen
            pause(2)
        end
        
        Screen(w,'FillRect',128);  % Gray screen
        Screen(w,'Flip');        % write to screen
        %% Compute the values for the trial
        
        if Trial > min([10 NTrCompDiff])
            LastFive10Ang = find((TypeOfTrial(1:Trial)==1 & abs(DesiredMeansEachTrial(1:Trial)) == 10),1+NTrCompDiff,'last');
            LastFive10Ang = LastFive10Ang(1:end-1);
            LastFive20Ang = find((TypeOfTrial(1:Trial)==1 & abs(DesiredMeansEachTrial(1:Trial)) == 20),1+NTrCompDiff,'last');
            LastFive20Ang = LastFive20Ang(1:end-1);

            if length(LastFive10Ang) >= NTrCompDiff & TypeOfTrial(Trial)
                RecentMeanCorrect10 = mean(CorrectResponse(LastFive10Ang));
                if RecentMeanCorrect10 >= .90
                    DifficultyFactor10 = min([DifficultyFactor10+2 9]);
                elseif RecentMeanCorrect10 >= .65
                    DifficultyFactor10 = min([DifficultyFactor10+1 9]);
                else
                    DifficultyFactor10 = max([DifficultyFactor10-1 -12]);
                end
            end
            if length(LastFive20Ang) >= NTrCompDiff & TypeOfTrial(Trial)
                RecentMeanCorrect20 = mean(CorrectResponse(LastFive20Ang));
                if RecentMeanCorrect20 >= .90
                    DifficultyFactor20 = min([DifficultyFactor20+2 9]);
                elseif RecentMeanCorrect20 >= .80
                    DifficultyFactor20 = min([DifficultyFactor20+1 9]);
                else
                    DifficultyFactor20 = max([DifficultyFactor20-1 -12]);
                end
            end
        end
        
        
        DifficultyFactorTrial10(Trial) = DifficultyFactor10;
        DifficultyFactorTrial20(Trial) = DifficultyFactor20;
        if abs(DesiredMeansEachTrial(Trial)) == 10
            DesiredMeansCorrectedEachTrial(Trial) = DesiredMeansCorrectedEachTrial(Trial)*(1 - DifficultyFactor10/10);
        elseif abs(DesiredMeansEachTrial(Trial)) == 20
            DesiredMeansCorrectedEachTrial(Trial) = DesiredMeansCorrectedEachTrial(Trial)*(1 - DifficultyFactor20/10);
        end
        mu = DesiredMeansCorrectedEachTrial(Trial);
        
        %             mu = InitialDesiredMean * SignMeanEachTrial(Trial);
        %         mu = DesiredMeansEachTrial(Trial);
        kappa = DesiredConcentrationEachTrial(Trial);
        
        numSamplesThisTrial = NumSamplesEachTrial(Trial);
        [NewValues, NewMu, NewKappa] = DrawArrayCircularINT(mu, kappa,numSamplesThisTrial, maxMuDevCrit, maxKappaDevCrit);
        AllAngles(1:numSamplesThisTrial,Trial) = NewValues;
        AllAnglesDegreesObserved(1:numSamplesThisTrial,Trial) = rem(NewValues + TrialBoundaryAngles(Trial),360);
        
        
        RelativeMeanEachTrial(1,Trial) = rem(NewMu,360);
        PositiveMeanTrial = rem(NewMu + 360*5,360);
        ObservedConcentrationEachTrial(1,Trial) = NewKappa;
        % Create the cue with respect to the reference angle of each trial
        [CueBluePink, SizeCue] = CreateBluePinkCueINT(DiamCircle,ExtraWidthPercentage, TrialBoundaryAngles(Trial));
        
        if Trial == 1
            %         CueBluePinkAllTrials = nan(SizeCue,SizeCue,NumTrials);
        end
        %     CueBluePinkAllTrials(:,:,Trial) = CueBluePink;
        NoGaborR = CueBluePink(:,:,1);
        NoGaborG = CueBluePink(:,:,2);
        NoGaborB = CueBluePink(:,:,3);
        NoGaborR(CentreCueCircle)= 128+CentralCircle(IsCentralCircle);
        NoGaborG(CentreCueCircle)= 128+CentralCircle(IsCentralCircle);
        NoGaborB(CentreCueCircle)= 128+CentralCircle(IsCentralCircle);
        NoGaborRGB(:,:,1) = NoGaborR;
        NoGaborRGB(:,:,2) = NoGaborG;
        NoGaborRGB(:,:,3) = NoGaborB;
        NoGaborRGBAllTrials = NoGaborRGB;
        
        MaskR = CueBluePink(:,:,1);
        MaskG = CueBluePink(:,:,2);
        MaskB = CueBluePink(:,:,3);
        MaskR(CentreCueCircle) = 128+(mask(IsCentralCircle)*255)+CentralCircle(IsCentralCircle);
        MaskG(CentreCueCircle) = 128+(mask(IsCentralCircle)*255)+CentralCircle(IsCentralCircle);
        MaskB(CentreCueCircle) = 128+(mask(IsCentralCircle)*255)+CentralCircle(IsCentralCircle);
        MaskRGB(:,:,1) = MaskR;
        MaskRGB(:,:,2) = MaskG;
        MaskRGB(:,:,3) = MaskB;
        
        
        for CurrSamp = 1:NumSamplesEachTrial(Trial)
            
            if CurrSamp == 1
                
                %Generate NoisyPatch even if integration trial
%                 Noise = (rand(300)-.5)*NoiseFactorDetection;
                [Noise] = CreateSmoothedNoise(300,10,30);
                
            end
            
            if TypeOfTrial(Trial) == 0
                if CurrSamp == NumSamplesEachTrial(Trial)
                    AllAnglesDegreesObserved(CurrSamp,Trial) = round(rand(1)*360);
                    GaborSamp = SCreateGaborINT(DiamCircle,GaborFactor*DiamCircle,AllAnglesDegreesObserved(CurrSamp,Trial),GaborFreq,0,60)*255;%ContrastDetection
                    NoisyPatch = GaborSamp;
                    if AppearGaborInNoise(Trial)
                        NoisyPatch(IsCentralCircle) = ...
                            (NoisyPatch(IsCentralCircle)+Noise(IsCentralCircle));
                    else
                        NoisyPatch(IsCentralCircle) = Noise(IsCentralCircle);
                    end
                    GaborSamp = NoisyPatch;
                else
                    GaborSamp = SCreateGaborINT(DiamCircle,GaborFactor*DiamCircle,AllAnglesDegreesObserved(CurrSamp,Trial),GaborFreq,0,DesiredContrastsEachTrial(1,Trial))*255;
                    
                end
            else
                GaborSamp = SCreateGaborINT(DiamCircle,GaborFactor*DiamCircle,AllAnglesDegreesObserved(CurrSamp,Trial),GaborFreq,0,DesiredContrastsEachTrial(1,Trial))*255;            
            end
            
            
            GaborSquare = 128+((GaborSamp .* patch) + CentralCircle);
            
            %             For computing the contrast:
            %             GforRMS = GaborSamp(IsCentralCircle);
            %             RMS = rms(GforRMS(:))
            
            
            GaborSampR = CueBluePink(:,:,1);
            GaborSampG = CueBluePink(:,:,2);
            GaborSampB = CueBluePink(:,:,3);
            
            GaborSampR(CentreCueCircle)= GaborSquare(IsCentralCircle);
            GaborSampG(CentreCueCircle)= GaborSquare(IsCentralCircle);
            GaborSampB(CentreCueCircle)= GaborSquare(IsCentralCircle);
            
            GaborRGB(:,:,1) = GaborSampR;
            GaborRGB(:,:,2) = GaborSampG;
            GaborRGB(:,:,3) = GaborSampB;
            GaborRGBAllTrials{CurrSamp,1} = GaborRGB;
        end % End of samples for current trial
        
        %% Store important variables in DATA structure
        data.allanglesobserved = AllAnglesDegreesObserved;
        
        %Check what category each trial belongs to: Blue = 1 Pink = 0
        IsBlueAngle(Trial) = (PositiveMeanTrial > 0) & (PositiveMeanTrial < 90) |...
            (PositiveMeanTrial > 180) & (PositiveMeanTrial < 270);
        

        
        
        TrialStartTime = GetSecs;% Record the beginning of the Trial
        %Print the initial Setting of Trial and wait for trial to begin
        Screen(w,'FillRect',128);
        Screen(w,'FillRect',[0 0 255],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
        Screen(w, 'FillOval', 0, FixationPoint);
        Screen(w,'FrameRect',0,RewardBox,3);%Add Rectangle Box
        Screen(w,'PutImage',NoGaborRGB,CueRectangle);%Add Circle in the CENTRAL RECTANGLE
        Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
        Screen(w,'Flip');%Print to Screen
        WaitSecs(1);%Give time for subjects to identify beginning of trial
        
        OnsetTimeTrial = GetSecs;%Record the time of Trial Onset
        
        %Show one or two masks
        
        %         for j = 1: NumMasksIni(Trial)
        IniTimeMask = GetSecs;
        ElapsedTimeMask = 0;
        %             if j == 2
        %                 WaitSecs(InterMaskDur)
        %             end
        Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
        Screen(w,'Framerect',0,RewardBox,3);
        Screen(w,'PutImage',MaskRGB,CueRectangle);
        Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
        Screen(w,'Flip');%Print to screen after an interval after Mask...
        FlagMaskEnd = (SampDur*1000);
        while ElapsedTimeMask < FlagMaskEnd
            ElapsedTimeMask = ceil(1000*(GetSecs - IniTimeMask));
        end
        Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
        Screen(w,'Framerect',0,RewardBox,3);
        Screen(w,'PutImage',NoGaborRGB,CueRectangle);%Add Circle in the CENTRAL RECTANGLE
        %    Screen(w,'FillOval',255,tinyrect);
        Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
        Screen(w,'Flip');%Print to screen at a point relative to Trial Onset
        %         end
        
        
        AfterMaskIni = GetSecs;
        AfterMaskTime = 0;
        
        while AfterMaskTime < InterSampInterval
            AfterMaskTime = GetSecs-AfterMaskIni;
        end
        
        
        ElapsedMilisecs = 0;
        CurrentSample = 1;
        RunningSamp = 1;
        
        while  RunningSamp <= NumMaxSamplesStream %Run until signal EndOfTrial
            
            if ElapsedMilisecs == 0;
                SampleStartTime = GetSecs;% Record the beginning of the sample presentation
            end
            
            ElapsedMilisecs = ceil(1000*(GetSecs - SampleStartTime));%Compute the Elapsed Time of Trial
            
            if ElapsedMilisecs < TimeStartSamp % Leave the screen without Stimulus for 150ms
                ElapsedMilisecs = ceil(1000*(GetSecs - SampleStartTime));%Compute the Elapsed Time of Trial);
                Screen(w,'FillRect',128);  % blank screen
                Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                %category
                Screen(w,'Framerect',0,RewardBox,3);
                Screen(w,'PutImage',NoGaborRGB,CueRectangle); %128+CentralCircle,CentralRectangle);
                Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
                Screen(w,'Flip');
                
            elseif ElapsedMilisecs >= TimeStartSamp && ElapsedMilisecs < TimeEndSamp
                
                ElapsedMilisecs = ceil(1000*(GetSecs - SampleStartTime));%Compute the Elapsed Time of Sample Presentation);
                Screen(w,'FillRect',128);  % blank screen
                GaborRGBSamp = GaborRGBAllTrials{CurrentSample,1}(:,:,:);
                Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                Screen(w,'Framerect',0,RewardBox,3);
                Screen(w,'PutImage',GaborRGBSamp,CueRectangle);
                Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
                Screen(w,'Flip');
%                 pause
                
            elseif ElapsedMilisecs >= TimeEndSamp
                %pause
                ElapsedMilisecs = 0;%Compute the Elapsed Time of Trial);
                Screen(w,'FillRect',128);  % blank screen
                Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                Screen(w,'Framerect',0,RewardBox,3);
                Screen(w,'PutImage',NoGaborRGB,CueRectangle); %128+CentralCircle,CentralRectangle);
                Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
                Screen(w,'Flip');
                RunningSamp = RunningSamp +1;
                if NumSamplesEachTrial(Trial) == CurrentSample
                    
                    break
                end
                CurrentSample = CurrentSample + 1;
            end
            
        end
        
        AfterMaskIni = GetSecs;
        AfterMaskTime = 0;
        
        while AfterMaskTime < InterSampInterval
            AfterMaskTime = GetSecs - AfterMaskIni;
        end
        
        
        if TypeOfTrial(Trial)
            %         %Show one or two masks
            %         for j = 1: NumMasksEnd(Trial)
            IniTimeMask = GetSecs;
            ElapsedTimeMask = 0;
            %             if j == 2
            %                 WaitSecs(InterMaskDur)
            %             end
            Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
            %category
            Screen(w,'Framerect',0,RewardBox,3);
            %             Screen(w,'FillOval',255,tinyrect);
            Screen(w,'PutImage',MaskRGB,CueRectangle);
            Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
            Screen(w,'Flip');%Print to screen after an interval after Mask...
            FlagFinishMask = (SampDur + InterSampInterval) *1000;
            while ElapsedTimeMask < FlagFinishMask
                ElapsedTimeMask = ceil(1000*(GetSecs - IniTimeMask));
            end
            Screen(w,'FillRect',128);  % Gray screen
            Screen(w,'PutImage',NoGaborRGB,CueRectangle);
            Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
            %    Screen(w,'FillOval',255,tinyrect);
            Screen(w,'Flip');%Print to screen at a point relative to Trial Onset
            %         end
            
        else
            
        end
        %% End Of Stimulus Presentations Beginning of Variable Gap
        EndOfStimPresentation = GetSecs;
        ElapsedTimeGap = 0;
        while ElapsedTimeGap < (DurGap(Trial))
            ElapsedTimeGap = GetSecs - EndOfStimPresentation;
        end
        
        %% End of Gap
        
        
        if TypeOfTrial(Trial) == 1 %% If Integration
            
            TimeEndOfTrial = GetSecs; % Register the end of trial
            
            Screen(w,'FillRect',128);  % blank screen
            if AppearBlueCircleRight(Trial) == 1
                Screen(w, 'FillOval', [0   128 255], ResponseRectRight);
                Screen(w, 'FillOval', [255 0 128], ResponseRectLeft);
            else
                Screen(w, 'FillOval', [255 0 128], ResponseRectRight);
                Screen(w, 'FillOval', [0   128 255], ResponseRectLeft);
            end
            %category
            Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
            Screen(w,'PutImage',NoGaborRGB,CueRectangle);
            Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
            Screen(w,'Framerect',0,RewardBox,3);
            
            
            Screen(w,'Flip');
            
            WaitPressButton = 0;
            
            ReactionTime = 9999;
            
            %Wait for Subject to press a Button to report response...
            while WaitPressButton==0;
                
                [kdown secs keyCode]=KbCheck;  % check for key press
                if kdown==1;
                    if keyCode(27)==1 || keyCode(41)==1;  % if escape key
                        Screen('closeall');% Close Screen
                        ShowCursor;%
                        rethrow(lasterror);%InformError
                        save datatmp data;%Save data acquired so far
                        break    % exit program
                    end
                    
                    if keyCode(39)==1 || keyCode(79)==1;
                        ReactionTime = GetSecs-TimeEndOfTrial;  % log RT
                        ResponseRightSide(Trial) = 1;
                        WaitPressButton = 1;
                        disp(['key: RIGHTDSIDE,',num2str(ReactionTime*1000),' ms']);
                    elseif keyCode(37)==1 || keyCode(80)==1;
                        ReactionTime = GetSecs-TimeEndOfTrial;  % log RT
                        ResponseRightSide(Trial) = 0;
                        WaitPressButton = 1;
                        disp(['key: LEFTSIDE,',num2str(ReactionTime),' ms']);
                    end
                end
            end
            
            if RelativeMeanEachTrial(Trial) ~=0
                if AppearBlueCircleRight(Trial) == ResponseRightSide(Trial)
                    WasResponseBlue = true;
                else
                    WasResponseBlue = false;
                end
                if WasResponseBlue == IsBlueAngle(Trial)
                    CorrectResponse(Trial) = 1;
                    
                    Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                    Screen(w,'Framerect',0,RewardBox,3);
                    Screen(w,'FillRect',[0 0 255],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
                    Screen(w,'PutImage',NoGaborRGB,CueRectangle);
                    Screen(w, 'FillOval', [0 255 0], FixationPoint);
                    Screen(w,'Flip');WaitSecs(0.1);
                    disp(['mean=',num2str(RelativeMeanEachTrial(Trial)),' correct']);
                    TotalReward=TotalReward+(CorrectResponse(Trial))*5;
                    TotalRewardLimited = TotalRewardLimited +(CorrectResponse(Trial))*5;

                    TotalRewardEachTrial(Trial) = TotalReward;
                    
                else
                    CorrectResponse(Trial) = 0;
                    Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                    Screen(w,'Framerect',0,RewardBox,3);
                    Screen(w,'FillRect',[0 0 255],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
                    Screen(w,'PutImage',NoGaborRGB,CueRectangle);
                    Screen(w, 'FillOval', [255 0 0], FixationPoint);
                    Screen(w,'Flip');
                    WaitSecs(0.1);
                    disp(['mean=',num2str(RelativeMeanEachTrial(Trial)),' error']);
                    TotalReward = TotalReward -(CorrectResponse(Trial))*10;
                    TotalRewardLimited = TotalRewardLimited - (CorrectResponse(Trial))*10;
                    TotalRewardEachTrial(Trial) = TotalReward;
                                        
                    
                end
                
                
                Screen(w,'Framerect',0,RewardBox,3);
                TotalRewardLimited = max([0 TotalRewardLimited]);
                TotalRewardLimited = min([TotalRewardLimited RewardBox(3)]);
                RewardBoxFilling(3) = min([RewardBox(3) RewardBox(1)+TotalRewardLimited]);
                
                if CorrectResponse(Trial) == 1
                    Screen(w,'Fillrect',[0 255 0],RewardBoxFilling);
                else
                    Screen(w,'FillRect',[255 0 0],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
                    
                end
                Screen(w,'PutImage',NoGaborRGB,CueRectangle);
                Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
                Screen(w,'Flip');
                WaitSecs(1);
                
                
            end
            
            
            
        else %If TrialType ~=1 If DETECTION
            
            TimeEndOfTrial = GetSecs; % Register the end of trial
            ResponseSectionDetection = 0;
            
            Screen(w,'FillRect',128);  % blank screen
            if AppearBlueCircleRight(Trial) == 1
                Screen('TextSize', w , floor(ScreenSize(2)*.06) );% Set the font size
                Screen('TextFont', w, 'Calibri');% Set the font of text
                Screen(w,'DrawText','PRESENT',(CentreOfScreen(1))+ 300,(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
                Screen(w,'DrawText','ABSENT',(CentreOfScreen(1))- 350,(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
            else
                Screen('TextSize', w , floor(ScreenSize(2)*.06) );% Set the font size
                Screen('TextFont', w, 'Calibri');% Set the font of text
                Screen(w,'DrawText','PRESENT',(CentreOfScreen(1)) + 250,(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
                Screen(w,'DrawText','ABSENT',(CentreOfScreen(1)) - 310,(CentreOfScreen(2)),0);%Inform Subject of PRACTICE session
            end
            %category
            Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
            Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
            Screen(w,'Framerect',0,RewardBox,3);
            
            
            Screen(w,'Flip');
            
            WaitPressButton = 0;
            
            ReactionTime = 9999;
            
            %Wait for Subject to press a Button to report response...
            while WaitPressButton==0;
                
                [kdown secs keyCode]=KbCheck;  % check for key press
                if kdown==1;
                    if keyCode(27)==1 || keyCode(41)==1;  % if escape key
                        Screen('closeall');% Close Screen
                        ShowCursor;%
                        rethrow(lasterror);%InformError
                        save datatmp data;%Save data acquired so far
                        break    % exit program
                    end
                    
                    if keyCode(39)==1 || keyCode(79)==1;
                        ReactionTime = GetSecs-TimeEndOfTrial;  % log RT
                        ResponseRightSide(Trial) = 1;
                        WaitPressButton = 1;
                        disp(['key: RIGHTDSIDE,',num2str(ReactionTime*1000),' ms']);
                    elseif keyCode(37)==1 || keyCode(80)==1;
                        ReactionTime = GetSecs-TimeEndOfTrial;  % log RT
                        ResponseRightSide(Trial) = 0;
                        WaitPressButton = 1;
                        disp(['key: LEFTSIDE,',num2str(ReactionTime),' ms']);
                    end
                end
            end
            
            if  ResponseRightSide(Trial)
                WasResponsePresent = true;
            else
                WasResponsePresent = false;
            end
            if WasResponsePresent == AppearGaborInNoise(Trial)
                CorrectResponse(Trial) = 1;
                
                Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                Screen(w,'Framerect',0,RewardBox,3);
                Screen(w,'FillRect',[0 0 255],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
                Screen(w, 'FillOval', [0 255 0], FixationPoint);
                Screen(w,'Flip');WaitSecs(0.1);
                disp(['mean=',num2str(RelativeMeanEachTrial(Trial)),' correct']);
                TotalReward=TotalReward+(CorrectResponse(Trial))*5;
                TotalRewardLimited = TotalRewardLimited +(CorrectResponse(Trial))*5;
                
                TotalRewardEachTrial(Trial) = TotalReward;
                
            else
                CorrectResponse(Trial) = 0;
                Screen(w,'Fillrect',[0 0 255],RewardBoxFilling);
                Screen(w,'Framerect',0,RewardBox,3);
                Screen(w,'FillRect',[0 0 255],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
                Screen(w,'PutImage',NoGaborRGB,CueRectangle);
                Screen(w, 'FillOval', [255 0 0], FixationPoint);
                Screen(w,'Flip');
                WaitSecs(0.1);
                disp(['mean=',num2str(RelativeMeanEachTrial(Trial)),' error']);
                TotalReward = TotalReward -(CorrectResponse(Trial))*10;
                TotalRewardLimited = TotalRewardLimited - (CorrectResponse(Trial))*10;
                TotalRewardEachTrial(Trial) = TotalReward;
                
                
            end
            
            data.RecordOfNoisyPatch{Trial,1} = GaborSquare;
            
            Screen(w,'Framerect',0,RewardBox,3);
            TotalRewardLimited = max([0 TotalRewardLimited]);
            TotalRewardLimited = min([TotalRewardLimited RewardBox(3)]);
            RewardBoxFilling(3) = min([RewardBox(3) RewardBox(1)+TotalRewardLimited]);
            
            if CorrectResponse(Trial) == 1
                Screen(w,'Fillrect',[0 255 0],RewardBoxFilling);
            else
                Screen(w,'FillRect',[255 0 0],RewardBoxFilling);%Add REWARD FILLING rectangle, originally set to no width
                
            end
            Screen(w,'PutImage',NoGaborRGB,CueRectangle);
            Screen(w, 'FillOval', [fixColour fixColour fixColour], FixationPoint);
            Screen(w,'Flip');
            WaitSecs(1);
            
            
            
            
            
        end
        
        
        Screen(w,'FillRect',128);
        Screen(w,'Flip');
        WaitSecs(.5);
        
        
        TrialsLengthTime(Trial) = GetSecs - TrialStartTime;
        
        data.difficultyfactor10 = DifficultyFactorTrial10;

        data.difficultyfactor20 = DifficultyFactorTrial20;

        data.desiredmeansdiffcorrected = DesiredMeansCorrectedEachTrial;
        
        %Response right (1) or Left (0)
        data.responserightside = ResponseRightSide;
        
        %Response correct (1) or Incorrect(0)
        data.correctresponse = CorrectResponse;
        
        %Presence or Absence of Gabor (1) or (0)
        data.appeargaborinnoise = AppearGaborInNoise;
    
        % Integration or Discrimination
        data.typeoftrial(Trial) = TypeOfTrial(Trial);
        
        % Number of Samples for each trial
        data.numsampsize(Trial) = NumSamplesEachTrial(Trial);
        
        %Desired Mean for each trial
        data.desiredmean(Trial) = DesiredMeansEachTrial(Trial);
        
        %Observed Mean for each trial without adding reference
        data.relativemean(Trial) = RelativeMeanEachTrial(Trial);
        
        %Observed Mean for each sample in Degrees having added the reference
        data.observedmeansample(Trial) = AllAnglesDegreesObserved(Trial);
        
        %Desired Concentration for each trial
        data.desiredconcentration(Trial) = DesiredConcentrationEachTrial(Trial);
        
        %Observed Concentration for each trial
        data.observedconcentration(Trial) = ObservedConcentrationEachTrial(Trial);
        
        %Contrast for each trial (useless...)
        data.desiredcontrasts(Trial) = DesiredContrastsEachTrial(Trial);
        
        % Boundary angles for each trial
        data.boundaryangles(Trial) =  TrialBoundaryAngles(Trial);
        
        % Position of Blue response probe 1 = right
        data.appearblueright(Trial) = AppearBlueCircleRight(Trial);
        
        % Category of each trial: 1 = Blue mean 0 = Pink mean
        data.isblueangle(Trial) = IsBlueAngle(Trial);
        
        %Number of masks at the beginning of trial
        data.nummasksini(Trial) = NumMasksIni(Trial);
        
        %Number of masks at the end of trial
        data.nummasksend(Trial) = NumMasksEnd(Trial);
        
        %Before Response Variable Gap for each trial
        data.durgap(Trial) = DurGap(Trial);
        
  

        data.reactiontime(Trial) = ReactionTime; 
        data.totreward = TotalRewardEachTrial;
        data.trialslengthtime(Trial) = TrialsLengthTime(Trial);
        
    end
    
    Screen(w,'FillRect',128);
    Screen(w,'DrawText','Block Completed',CentreOfScreen(1)-50,CentreOfScreen(2)+50,0);
    Screen(w,'DrawText','Saving Data...',CentreOfScreen(1)-50,CentreOfScreen(2)-50,0);% draw instuctions
    Screen(w,'Flip');
    save(NameFileToSave, 'data') ;
    OverallTime = toc/60;
    
catch
    ShowCursor;%
    Screen('closeall');% Close Screen
    OverallTime = toc

%     save datatmp data;
    rethrow(lasterror);%InformError
    
end


Screen('closeall');% Close Screen
ShowCursor;%

end


%% Function called for creating Gabor Patches note that NormalPDF is in capitals

function [patch] = SCreateGaborINT(siz,envelopedev,angle,frequency,phase,contrast)
%  [patch] = CreateGabor(siz,envelopedev,angle,frequency,phase,[contrast])


if nargin < 6, contrast = 1; end
if nargin < 5, error('Not enough input arguments.'); end

siz = floor(siz/2)*2;
[xint,yint] = meshgrid((1:siz)-(siz+1)/2);%xint and yint are the internal versions
patch = 0.5*contrast*cos(2*pi*(frequency*(sin(pi/180*angle)*xint+cos(pi/180*angle)*yint)+phase));
patch = patch.*NormalPDF(sqrt(xint.^2+yint.^2),0,envelopedev);

end

function [patch] = CreateCircleINT(siz,width,falloff)
%  [patch] = CreateCircle(siz,width,[falloff])

if nargin < 3, falloff = 2; end
if nargin < 2, error('Not enough input arguments.'); end

sigmoidfun = @(x,lims)lims(1)+diff(lims)./(1+exp(-x));

siz = floor(siz/2)*2;
[x,y] = meshgrid((1:siz)-(siz+1)/2);

coef = log(1/0.01-1)*2/falloff;

patch = sigmoidfun(coef*(sqrt(x.^2+y.^2)-siz/2),[1,0]);
patch = min(patch,sigmoidfun(coef*(sqrt(x.^2+y.^2)-(siz/2-width)),[0,1]));

end

function [patch] = CreateCircularApertureINT(siz,falloff)
%  [patch] = CreateCircularApertureINT(siz,[falloff])

if nargin < 2, falloff = 2; end
if nargin < 1, error('Not enough input arguments.'); end

sigmoidfun = @(x,lims)lims(1)+diff(lims)./(1+exp(-x));

siz = floor(siz/2)*2;
[x,y] = meshgrid((1:siz)-(siz+1)/2);

coef = log(1/0.01-1)*2/falloff;

patch = sigmoidfun(coef*(sqrt(x.^2+y.^2)-siz/2),[1,0]);

end

function [CueBluePink, SizeCue] = CreateBluePinkCueINT(DiamCircle, ExtraWidthPercentage, Angle)

NewDiamCircle = DiamCircle * (1 + ExtraWidthPercentage);
siz = floor(NewDiamCircle/2)*2;
NewDiamCircle = siz;
angle = (45+Angle)*pi/180;

CorrectionX = sin(angle);
CorrectionY = cos(angle);

Circulo = CreateCircularApertureINT(NewDiamCircle);

[x,y] = meshgrid((1:siz)-(siz+1)/2);

x = Circulo.*x;
y = Circulo.*y;

z = abs((-x*CorrectionX) + (-y*CorrectionY));
minZ = abs(min(min(z)));
maxZ = max(max(z));
z = (z + (minZ))/ (maxZ + minZ);


CueBluePink = zeros(NewDiamCircle,NewDiamCircle,3);
CueBluePink (:,:,1) = 256*z;
CueBluePink (:,:,2) = 128*rot90(z,1);
CueBluePink (:,:,3) = 128*rot90(z,1) + 128*Circulo;

SumCue = mean(CueBluePink,3);
IsZeroCue = SumCue <128/2;
logic3d = false(NewDiamCircle,NewDiamCircle,3);

for i = 1:3
    logic3d(:,:,i) = IsZeroCue;
end

CueBluePink(logic3d) = 128;
SizeCue = length(CueBluePink);

end

function[NewValues, NewMu, NewKappa] = DrawArrayCircularINT(mu, kappa,numElements, maxMuDevCrit, maxKappaDevCrit)
% mu is the mean angle in degrees of the desired distribution
% kappa is the concentration (analogous to 1/variance) of the desired
% distribution
%numElements is the number of elements that the distribution should have
%maxMuDevCrit is the maximum deviation of the final mu value from the
%desired value. Note, this value is in degrees
%maxKappaDevCrit is the maximum deviation of the final kappa value from the
%desired value. Note, this value is in radians

InRadMu = pi*mu/180;
InRadmaxMuDevCrit = pi*maxMuDevCrit/180;



w = ones(numElements,1);
dim = 1;

SatisfactoryValues = 0;
while SatisfactoryValues == 0
    alpha = Scirc_vmrndINT(InRadMu, kappa, numElements);
    
    if  numElements < 3 % Compute a Fake estMu and kappa
        if numElements < 3
            estMu = Scirc_meanINT(alpha,w,dim);
        else
            [estMu, ~, ~] = Scirc_meanINT(alpha,w,dim);
        end
        estKappa = std(kappa);
        InRadmaxMuDevCrit = InRadmaxMuDevCrit + .2;
        maxKappaDevCrit = maxKappaDevCrit + .2;
        SatisfactoryMu = abs(Scirc_distINT(InRadMu, estMu)) < InRadmaxMuDevCrit;
        SatisfactoryKappa = true;
        
    else
        [estMu, ~, ~] = Scirc_meanINT(alpha,w,dim);
        estKappa = Scirc_kappaINT(alpha,w);
        SatisfactoryMu = abs(Scirc_distINT(InRadMu, estMu)) < InRadmaxMuDevCrit;
        SatisfactoryKappa = abs(kappa - estKappa) < maxKappaDevCrit;
        
    end
    
    if SatisfactoryMu && SatisfactoryKappa
        SatisfactoryValues = 1;
    end
end

NewValues = (180*alpha/pi);
% NewValues(NewValues < 0) = NewValues(NewValues < 0) + 360;
NewMu =  (180*estMu/pi);

NewKappa = estKappa;
end

function alpha = Scirc_vmrndINT(theta, kappa, n)

%Modified by Santiago from Circular Statistics Toolbox
%alpha = circ_vmrnd(theta, kappa, n)
%   Simulates n random angles from a von Mises distribution, with preferred
%   direction thetahat and concentration parameter kappa.
%
%   Input:
%     [theta    preferred direction, default is 0]
%     [kappa    width, default is 1]
%     [n        number of samples, default is 10]
%
%     If n is a vector with two entries (e.g. [2 10]), the function creates
%     a matrix output with the respective dimensionality.
%
%   Output:
%     alpha     samples from von Mises distribution
%
%
%   References:
%     Statistical analysis of circular data, Fisher, sec. 3.3.6, p. 49
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu


% default parameter
if nargin < 3
    n = 10;
end

if nargin < 2
    kappa = 1;
end

if nargin < 1
    theta = 0;
end

if numel(n) > 2
    error('n must be a scalar or two-entry vector!')
elseif numel(n) == 2
    m = n;
    n = n(1) * n(2);
end

% if kappa is small, treat as uniform distribution
if kappa < 1e-6
    alpha = 2*pi*rand(n,1);
    return
end

% other cases
a = 1 + sqrt((1+4*kappa.^2));
b = (a - sqrt(2*a))/(2*kappa);
r = (1 + b^2)/(2*b);

alpha = zeros(n,1);
for j = 1:n
    while true
        u = rand(3,1);
        
        z = cos(pi*u(1));
        f = (1+r*z)/(r+z);
        c = kappa*(r-f);
        
        if u(2) < c * (2-c) || ~(log(c)-log(u(2)) + 1 -c < 0)
            break
        end
        
        
    end
    
    alpha(j) = theta +  sign(u(3) - 0.5) * acos(f);
    alpha(j) = angle(exp(1i*alpha(j)));
end

if exist('m','var')
    alpha = reshape(alpha,m(1),m(2));
end
end

function [mu ul ll] = Scirc_meanINT(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
    t = Scirc_confmeanINT(alpha,0.05,w,[],dim);
    ul = mu + t;
    ll = mu - t;
end
end

function kappa = Scirc_kappaINT(alpha,w)
%
% kappa = circ_kappa(alpha,[w])
%   Computes an approximation to the ML estimate of the concentration
%   parameter kappa of the von Mises distribution.
%
%   Input:
%     alpha   angles in radians OR alpha is length resultant
%     [w      number of incidences in case of binned angle data]
%
%   Output:
%     kappa   estimated value of kappa
%
%   References:
%     Statistical analysis of circular data, Fisher, equation p. 88
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


alpha = alpha(:);

if nargin<2
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) > size(w,1)
        w = w';
    end
end

N = length(alpha);

if N>1
    R = Scirc_rINT(alpha,w);
else
    R = alpha;
end

if R < 0.53
    kappa = 2*R + R^3 + 5*R^5/6;
elseif R>=0.53 && R<0.85
    kappa = -.4 + 1.39*R + 0.43/(1-R);
else
    kappa = 1/(R^3 - 4*R^2 + 3*R);
end

if N<15 && N>1
    if kappa < 2
        kappa = max(kappa-2*(N*kappa)^-1,0);
    else
        kappa = (N-1)^3*kappa/(N^3+N);
    end
end
end

function t = Scirc_confmeanINT(alpha, xi, w, d, dim)
%
% t = circ_mean(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- d yields upper/lower (1-xi)% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 5
    dim = 1;
end

if nargin < 4 || isempty(d)
    % per default do not apply correct for binned data
    d = 0;
end

if nargin < 3 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
    xi = 0.05;
end

% compute ingredients for conf. lim.
r = Scirc_rINT(alpha,w,d,dim);
n = sum(w,dim);
R = n.*r;
c2 = chi2inv((1-xi),1);

% check for resultant vector length and select appropriate formula
t = zeros(size(r));

for i = 1:numel(r)
    if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
        t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
    elseif r(i) >= .9
        t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
    else
        t(i) = NaN;
        warning('Requirements for confidence levels not met.');
    end
end

% apply final transform
t = acos(t./R);

end

function r = Scirc_rINT(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

if nargin < 3 || isempty(d)
    % per default do not apply correct for binned data
    d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
    c = d/2/sin(d/2);
    r = c*r;
end

end


function r =  Scirc_distINT(x,y)
%
% r = circ_dist(alpha, beta)
%   Pairwise difference x_i-y_i around the circle computed efficiently.
%
%   Input:
%     alpha      sample of linear random variable
%     beta       sample of linear random variable or one single angle
%
%   Output:
%     r       matrix with differences
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 3/19/2009
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


if size(x,1)~=size(y,1) && size(x,2)~=size(y,2) && length(y)~=1
  error('Input dimensions do not match.')
end

r = angle(exp(1i*x)./exp(1i*y));
end