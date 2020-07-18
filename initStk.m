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
