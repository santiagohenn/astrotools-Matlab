clear; clc;
app = actxserver('STK10.application');
%STKXApplication = actxserver('STKX11.application');

% app.NoGraphics = 1;
% STKXApplication.NoGraphics = 1;

%app.Visible
%root = actxserver('AgStkObjects11.AgStkObjectRoot');
root = app.Personality2;

format long
scenario = root.Children.New('eScenario','DEFAULT');
%scenario = root.LoadScenario('C:\Users\Usuario\Desktop\STK Files\Anything\Anything.sc');
%scenario = root.CurrentScenario;
%root.UnitPreferences.Item('DateFormat').SetCurrentUnit('UTCG');
%root.ExecuteCommand('VO * SnapFrame SetValues Format jpeg');
%root.ExecuteCommand('VO * SnapFrame SetValues AntiAlias On 4');
scenario.SetTimePeriod('10 Mar 2020 16:00:05.000','11 Mar 2020 16:00:05.000');
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec');
root.Rewind;

%%

    nFT = 12;
    gridSpace = 360/12;

    for j=1:nFT
    facility = scenario.Children.New('ePlace',strcat('Eq_',num2str(j)));
    facility.Position.AssignGeodetic(0,(j-1)*gridSpace,0);
    end
    
    nFT = 6;
    gridSpace = 360/6;

    for j=1:nFT
    facility = scenario.Children.New('ePlace',strcat('P_',num2str(j)));
    facility.Position.AssignGeodetic(60,(j-1)*gridSpace,0);
    end
    
%%
        RE = 6378135;
        Hmax = 700000;
        th=deg2rad(40);
        etaMax = asin((RE*cos(th))/(RE+Hmax))
        etaMax = rad2deg(etaMax);
        
%         planesG = [3     5    5     5     5     5     5     5];
%         satsG =   [3     2    2     3     3     4     6     8];
%         incG =    [85.96 78.0 84.77 84.17 86.54 88.95 88.27 89];

        planes = 5;
        sats = 8;

		planePhase = 360/planes;
		satsPhase = 360/sats;

    for nPlane = 1:planes
        
        for nSat = 1:sats
            
            sat = scenario.Children.New('eSatellite',strcat('s',num2str(nSat),'_p',num2str(nPlane)));

            keplerian = sat.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical');
            sat.Propagator.step=120; % 60 segundos step 
            %keplerian.SizeShapeType = 'eSizeShapeAltitude';
            keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
            %keplerian.LocationType = 'eLocationTrueAnomaly';

            keplerian.LocationType = 'eLocationMeanAnomaly';
            keplerian.Orientation.AscNodeType = 1;  % 1=RAAN
            keplerian.SizeShape.SemiMajorAxis = (6378135+700000)/1000;
            keplerian.SizeShape.Eccentricity = 0;
            keplerian.Orientation.Inclination = 89;%82.38;
            keplerian.Orientation.AscNode.Value = (nPlane-1)*planePhase;
            keplerian.Orientation.ArgOfPerigee = 0;
            keplerian.Location.Value = (nSat-1)*satsPhase-(nPlane-1)*(planePhase/sats);
            sat.Propagator.InitialState.Representation.Assign(keplerian);
            sat.SetPropagatorType(1);   % 1 = J2 perturbation, 2 = J4, 7 = Two body.
            sat.Propagator.Propagate;
            accessConstraints = sat.AccessConstraints;
            minAngle = accessConstraints.AddConstraint('eCstrGroundElevAngle');
            minAngle.EnableMin = true;
            minAngle.Min = 5; % Grados
            minGrazingAlt = accessConstraints.AddConstraint('eCstrGrazingAlt');
            minGrazingAlt.EnableMin = true;
            minGrazingAlt.Min = 100; % Km
            sensor = sat.Children.New('eSensor',strcat('sensor_',num2str(nSat),'_plane_',num2str(nPlane)));
            sensor.CommonTasks.SetPatternSimpleConic(etaMax,1);
        end
        
    end

%%

scenario = root.CurrentScenario;
sat = scenario.Children.Item('SAOCOM_1A_43641');
accessConstraints = sat.AccessConstraints;
minAngle = accessConstraints.AddConstraint('eCstrGroundElevAngle');
minAngle.EnableMin = true;
minAngle.Min = 5; % Degrees

%%

timeElapsed = 0;
scenario.SetTimePeriod('20 Mar 2020 11:00:00.000','18 Jul 2020 11:00:00.000');

    tic
    sat.Propagator.Propagate;
    gs = scenario.Children.Item('Facility1');
    access = gs.GetAccess('Satellite/SAOCOM_1A_43641');
    access.ComputeAccess();
    intervalCollection = access.ComputedAccessIntervalTimes;
    timeElapsed = timeElapsed + toc

%%
    clc
    L = deg2rad(60);
    lambda = deg2rad(20);
    
    P = 0.258;

    result = 2*atan((sqrt(-sin(lambda)^2*sec(L)^2+tan(L)^2+pi^2*P^2)-pi*P)/(sin(lambda)*sec(L)+tan(L)));
    
    incDeg = rad2deg(result);
    inc = result;
    
    result2 = (1/pi)*((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)));
    

%%
    clc
    clear 
    series = zeros(1,120);
    
    figure
    
    for la = 10:10:30
    lambda = deg2rad(la);
        for L = 10:89
        inc = 2*atan(sin(lambda)*sec(deg2rad(L)));
        series(1,L) = rad2deg(inc);
        end
    plot(10:1:89,series(1,10:89)); hold on;
    end
    
    xlabel('Greatest latitude to cover');
    ylabel('Orbit inclination needed');
    legend('\lambda = 10','\lambda = 20','\lambda = 30','\lambda = 16','\lambda = 18','\lambda = 20','\lambda = 22','\lambda = 24','\lambda = 26','\lambda = 28','\lambda = 30');
    
    grid on
    
    % CORROBORARRR
    
    %%
    
    L = deg2rad(60);
    lambda = deg2rad(17.51);
    figure
    
    for i = 1:90
        
        inc = deg2rad(i);
        RL = real(acos((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)))/pi);
        Ro = real(1 - 2*acos((sin(lambda))/(inc))/pi);
        Ropt = RL-Ro;
        plot(i,Ropt,'k.'); hold on;
        
    end
    


%%
    clear
    clc
    figure
    L = deg2rad(60);
    lambda = deg2rad(22.13);
    
    for i = 1:90
        
        inc = deg2rad(i);
        RL = real(acos((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)))/pi);
        Ro = real(1 - (2*acos((sin(lambda))/sin(inc)))/pi);
        Ropt = RL-Ro;
        plot(i,Ropt,'k.'); hold on;
        
    end
   
    for i = 1:90
        
        inc = deg2rad(i);
        Ro = real(acos((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)))/pi) - real(1 - (2/pi)*acos((sin(lambda))/sin(inc)));
        Ropt = RL-Ro;
        plot(i,Ropt,'r.'); hold on;
        
        if(RL>Ro && Ro~=1)
            inc0 = i-1;
            inc1 = i;
            wdt = 0;
            
            while (abs(inc0-inc1) >= 0.01)
                
                incx = (inc1+inc0)/2;
                plot(incx,Ropt,'r.'); hold on;
                incRad = deg2rad(incx);
                RL = real(acos((-sin(lambda)+cos(incRad)*sin(L))/(sin(incRad)*cos(L)))/pi);
                Ro = real(1 - (2/pi)*acos((sin(lambda))/sin(incRad)));
                Ropt = RL-Ro;

                if(Ropt==0) 
                     break;
                end

                if(Ropt<0)
                    inc0 = incx;
                elseif (Ropt>0)
                    inc1 = incx;
                end
            
                if (wdt>1000)
                    break;
                end

            end
            break;
        end
        
    end
    
        grid on
        
        RL
        Ro
        incx

%%
    clc
    series = zeros(1,90);
    % inc = deg2rad(51.587);
    %lambda = deg2rad(30);
    %%figure
    
    RE = 6378135;
    Hmax = 700000;
    th=deg2rad(15);
    etaMax = asin((RE*cos(th))/(RE+Hmax))
    lambda = 90-5-rad2deg(etaMax)
    lambda = deg2rad(lambda);
    
%%
    figure 
    for i = 20:10:90
        inc = deg2rad(i);
    
    for la = 0:1:89
    L = deg2rad(la);
    phi1 = (-sin(lambda) + cos(inc)*sin(L))/(sin(inc)*cos(L));
        if (L<inc-lambda)
        	phi2 = (sin(lambda) + cos(inc)*sin(L))/(sin(inc)*cos(L));
            series(1,la+1)=(acos(phi1)-acos(phi2))*100/pi;
        elseif ((inc-lambda <= L) && (L <= inc + lambda))
            series(1,la+1)=acos(phi1)*100/pi;
        elseif (L > inc+lambda)
            series(1,la+1)=0;
        end
    end
    
        plot(0:1:89,series(1,1:90),'LineWidth',1.2); hold on; %,'color',rand(1,3)
    
    end
    
    title({'Percentage of orbits with coverage for a single point','in a given latitude'});
    xlabel('Latitude [degrees]');
    ylabel('percentage of orbits with coverage [0 - 1]');
    %legend({'i = 20','i = 25','i = 30','i = 35','i = 40','i = 45','i = 50','i = 55','i = 60','i = 65','i = 70'},'NumColumns',2);
    legend({'i = 20','i = 30','i = 40','i = 50','i = 60','i = 70'},'NumColumns',2);
    
    grid on
    
    % CORROBORARRR
    
%%

    figure;

    RE = 6378135;
    Hmax = 1000000;
    th = deg2rad(20);
    L = deg2rad(60);
    seriesTotal = zeros(90,3);
    seriesRo = zeros(4,90);
    seriesRL = zeros(4,90);
    h = 1;
    
    color = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
    
    for nPlot=1:3
    
        L = deg2rad(nPlot*25);
        h = 1;
        
    for Hmax = 400000:200000:1000000
        
    etaMax = asin((RE*cos(th))/(RE+Hmax));
    lam = 90-th-rad2deg(etaMax);
    lambda = deg2rad(lam);
    
    for i = 1:90
        inc = deg2rad(i);
        
        RL = real(acos((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)))/pi);
        Ro = real(1 - 2*acos((sin(lambda))/(inc))/pi);
        
        %seriesRo(h,i) = Ro;
        %seriesRL(h,i) = RL;
        seriesTotal(i,:) = [i RL Ro];
        
        if(RL>Ro && Ro~=1)
            inc0 = i-1;
            inc1 = i;
            
            while (abs(inc0 - inc1) >= deg2rad(0.1))

                incx = (inc1+inc0)/2;
                RL = real(acos((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)))/pi);
                Ro = real(1 - 2*acos((sin(lambda))/(inc))/pi);
                Ropt = RL-Ro;

                if(Ropt==0) 
                     break;
                end

                if(Ropt<0)
                    inc0 = incx;
                elseif (Ropt>0)
                    inc1 = incx;
                end

            end
            
            incx
            %plot(i,RL,'kx'); hold on;
            break;
        end
        
    end
        subplot(1,3,nPlot)
        p1 = plot(seriesTotal(1:i,1),seriesTotal(1:i,2),'Color',color(h,:),'LineStyle','--'); hold on;
        p2 = plot(seriesTotal(1:i,1),seriesTotal(1:i,3),'Color',color(h,:)); hold on;
        h = h+1;
        
    end
        hdl = zeros(4,1);
        hdl(1) = plot(1,seriesRo(1,1),'Color',color(1,:));
        hdl(2) = plot(1,seriesRo(2,1),'Color',color(2,:));
        hdl(3) = plot(1,seriesRo(3,1),'Color',color(3,:));
        hdl(4) = plot(1,seriesRo(4,1),'Color',color(4,:));
        legend(hdl,'H=400Km','H=600Km','H=600Km','H=1000Km');
        
        grid on
        %gtext('RL'); hold on;
        %gtext('Ro'); hold on;
        title(strcat('Lat= ',num2str(nPlot*25),' °'));
        xlabel('Inclination (degrees)');
        ylabel('Percentage of orbits with access');
    end

    %% MATLAB SOLVE
    
    clc
    syms x;
    % inc = deg2rad(47);    49.5188
    beta = deg2rad(60);
    alpha = deg2rad(12.13);
    assume(x>0 & x<pi);
    %assume(alpha>1 & alpha<pi/2);
    %assume(beta>0 & beta<pi/2);
    vpasolve(real((1/pi)*acos((-sin(alpha)+cos(x)*sin(beta))/(sin(x)*cos(beta)))) - real(1 - (2/pi)*acos(sin(alpha)/sin(x))) == 0,x)
    
    %% MATLAB SOLVE
    
    beta = deg2rad(60);
    alpha = deg2rad(12.13);
    solu = deg2rad(49.5217);
    Ro = 1 - (2/pi)*acos(sin(alpha)/sin(solu))
    RL = (1/pi)*acos((-sin(alpha)+cos(solu)*sin(beta))/(sin(solu)*cos(beta)))

    %% BISECTION SOLVE
    
    clear
    clc
    figure
    L = deg2rad(60);
    lambda = deg2rad(12.138928409178732);
    
    for i = 1:90
        
        inc = deg2rad(i);
        RL = real(acos((-sin(lambda)+cos(inc)*sin(L))/(sin(inc)*cos(L)))/pi);
        Ro = real(1 - (2*acos((sin(lambda))/sin(inc)))/pi);
        Ropt = RL-Ro;
        plot(i,Ropt,'k.'); hold on;
        
    end
        
        inc0 = 1;
        inc1 = 89;
            
    while (abs(inc0-inc1) >= 0.001)

        incx = (inc1+inc0)/2;
        incRad = deg2rad(incx);
        %RL = real(acos((-sin(lambda)+cos(incRad)*sin(L))/(sin(incRad)*cos(L)))/pi);
        %Ro = real(1 - (2/pi)*acos((sin(lambda))/sin(incRad)));
        
        RL = acos((-sin(lambda)+cos(incRad)*sin(L))/(sin(incRad)*cos(L)))/pi;
        Ro = 1 - (2/pi)*acos((sin(lambda))/sin(incRad));
        
        Ropt = RL-Ro;

        plot(incx,Ropt,'rX'); hold on;

        if(Ropt==0) 
             break;
        end

        if(Ropt<0)
            inc0 = incx;
        elseif (Ropt>0)
            inc1 = incx;
        end

    end
    
        grid on
        ylabel('Difference in percentage of access');
        xlabel('Inclination (degrees)');
        
        RL
        Ro
        incx