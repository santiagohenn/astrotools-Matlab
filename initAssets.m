%% Variables
RT = 6378.137; % Radio Terrestre (máximo) [Km]
uT = 3.986004418E14; % Earth's Standard Gravitational parameter

% [p1 p2 p3 ... pn     iteration
%  p1 p2 p3 ... pn     current move
%  p1 p2 p3 ... pn]    best found so far

limInf=0; limSup=0;
rng('shuffle'); fprintf('Loading Assets..')

%% STK root init
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('UTCG');
root.UnitPreferences.SetCurrentUnit('Distance','Km');
scenario.SetTimePeriod('20 Jul 2020 16:00:05.000','30 Jul 2020 16:00:05.000'); root.Rewind;
root.ExecuteCommand('SubObjUnload *');
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); % .SetCurrentUnit('JDate');
scTime=scenario.StopTime;

% Elements order: [Eccentricity, Semimajor axis, Inclination, Longitude of the ascending node, Argument of periapsis, True anomaly]
%% Mission Satellite/s
nSM=2;
parametrosSat(1,:,1)=[0.00001,7007.7,97.8892,193.4537,90,270];     
parametrosSat(1,:,2)=[0.00001,7007.7,97.8892,193.4537,90,90];

% Elements order: [Eccentricity, Semimajor axis, Inclination, Longitude of the ascending node, Argument of periapsis, True anomaly]
%% Relay satellite/s
nRelay=1;
parametrosSat(1,:,3)=[0,7378.137,45,193.4537,180,0];
% parametrosSat(1,:,5)=[0,1000,80,0,0,0];
% parametrosSat(1,:,6)=[0,1000,80,0,0,0];

% Parameters order: % [Index, Latitude, Longitude, Altitude]
%% Ground Assets
nGA=2; parametrosGA = zeros(3,4,nGA);
parametrosGA(1,:,1)=[1 -31.5263 -64.4674 0];
parametrosGA(1,:,2)=[1 -54.5104 -67.1164 0];

%% Populate in STK

for nAsset=1:nSM+nRelay
    if nAsset<=nSM
        parametrosTemp=parametrosSat(1,:,nAsset);
        sat = scenario.Children.New('eSatellite',strcat('MS_',num2str(nAsset)));
    else
        parametrosTemp=parametrosSat(1,:,nAsset);
        sat = scenario.Children.New('eSatellite',strcat('RELAY_',num2str(nAsset-nSM)));
    end
keplerian = sat.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical');
sat.Propagator.step=60; % 60 segundos step 
%keplerian.SizeShapeType = 'eSizeShapeAltitude';
keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
%keplerian.LocationType = 'eLocationTrueAnomaly';
keplerian.LocationType = 'eLocationMeanAnomaly';
keplerian.Orientation.AscNodeType = 1;  % 1=RAAN, 0=LAN
eccentricity = parametrosTemp(1);
%semiMajorAxis = parametrosTemp(2);
keplerian.SizeShape.Eccentricity = parametrosTemp(1);
keplerian.SizeShape.SemiMajorAxis = parametrosTemp(2);
%keplerian.SizeShape.PerigeeAltitude = semiMajorAxis*(1-eccentricity);
%keplerian.SizeShape.ApogeeAltitude = semiMajorAxis*(1+eccentricity);
keplerian.Orientation.Inclination = parametrosTemp(3);
keplerian.Orientation.AscNode.Value = parametrosTemp(4);
keplerian.Orientation.ArgOfPerigee = parametrosTemp(5);
% keplerian.SizeShape.MeanMotion=parametrosTemp(6);
keplerian.Location.Value = parametrosTemp(6);
% keplerian.Orientation.MeanAnomaly = parametrosTemp(6);
sat.Propagator.InitialState.Representation.Assign(keplerian);
sat.Propagator.Propagate;
accessConstraints = sat.AccessConstraints;
minAngle = accessConstraints.AddConstraint('eCstrGroundElevAngle');
minAngle.EnableMin = true;
minAngle.Min = 5;   % Degrees
minGrazingAlt = accessConstraints.AddConstraint('eCstrGrazingAlt');
minGrazingAlt.EnableMin = true;
minGrazingAlt.Min = 100; % Km
end

%% Ground Assets
    for j=1:nGA
        facility = scenario.Children.New('eFacility',strcat('GA_',num2str(j)));
        facility.Position.AssignGeodetic(parametrosGA(1,2,j),parametrosGA(1,3,j),parametrosGA(1,4,j));
    end
    
%% Initial parameters

parametrosSat(3,:,:)=parametrosSat(1,:,:);
parametrosGA(3,:,:)=parametrosGA(1,:,:);

%%
fprintf('..Assets Loaded.'); fprintf('\n');

% SABIA MAR
%    Mean Motion (rev/day) = 14.56416259        a=(uT/(((2*pi*14.56416259)/(86400))^2))^(1/3)                                     
%    Mean Motion Dot       = 0.00000000                                              
%    Mean Motion Dot Dot   = 0.00000000e+000                                         
%    Bstar                 = 0.00000000e+000                                         
%    Inclination           = 98.22020000                                             
%    RAAN                  = 75.49950000                                             
%    Eccentricity          = 0.00001000                                              
%    Arg of Periapsis      = 90.00000000                                             
%    Mean Anomaly          = 89.87940000                                             
%    Rev Number            = 0                                                       
%    Element Set Number    = 1
%    Semi-major axis       = 707.6 + RT

% SAOCOM 1-A
%    Mean Motion (rev/day) = 14.82154588                                             
%    Mean Motion Dot       = 0.00000000                                              
%    Mean Motion Dot Dot   = 0.00000000e+000                                         
%    Bstar                 = 0.00000000e+000                                         
%    Inclination           = 97.88920000                                             
%    RAAN                  = 193.45370000                                            
%    Eccentricity          = 0.00001000                                              
%    Arg of Periapsis      = 90.00000000                                             
%    Mean Anomaly          = 270.12210000                                            
%    Rev Number            = 0                                                       
%    Element Set Number    = 625.4

% SAOCOM 1-B
%    Mean Motion (rev/day) = 14.82154588                                             
%    Mean Motion Dot       = 0.00000000                                              
%    Mean Motion Dot Dot   = 0.00000000e+000                                         
%    Bstar                 = 0.00000000e+000                                         
%    Inclination           = 97.88920000                                             
%    RAAN                  = 193.45370000                                            
%    Eccentricity          = 0.00001000                                              
%    Arg of Periapsis      = 90.00000000                                             
%    Mean Anomaly          = 89.87790000                                             
%    Rev Number            = 0                                                       
%    Element Set Number    = 625.4
    