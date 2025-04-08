% The purpose of this script is to collect environmental parameters, system parameters, 
% and signal-to-noise ratio data of the optical system for data fitting using a neural network.
% Modify the telescope type and the file path for saving the collected data.

t_start=clock;
% init ZOS Connect
if ~exist('args', 'var') 

    args = [];
end
% Initialize the OpticStudio connection
TheApplication = InitConnection();
if isempty(TheApplication)
    % failed to initialize a connection
    r = [];
end




%% Define parameters
temobs = Observer();
workingfolder = strcat('workingfolder');
telescopefolder = strcat('telescopetemplate');
TelEfficiency = 0.78;
temobs.TelEfficiency = TelEfficiency;
SNRLimit = 10;          
overallfov = 10000;     
Overlap = 5;
stackerror = 0.3;
      
SkyMag_ori = 22.3;
temobs.SkybackNoise = SkyMag_ori;

wavelength = 0.55;
temobs.Wavelength = wavelength;
width_height = [1,1];   
T_transform = 10;   

telescopelist = {'RC-CORRECT', };
Cameralist = {'CCD','CMOS'};

%% ---------
telType = 1; 
camType = 1;

% Telescopic aperture occlusion ratio
radio = [0.795, ];
D_radio = radio(telType); 

Fov_list = [2.0, ];
Fov = Fov_list(telType);                    
temobs.Fov = Fov;
fprintf('field is %g \n', Fov) 

telescopename = telescopelist{telType};
temobs = temobs.definecam(Cameralist{camType});


path1 = pwd;
[pathstr,name,suffix] = fileparts(telescopename);
outfile = strcat(path1,'\workingfolder\',name,'.zmx');
if contains(pathstr, ':') 
    orgfile = strcat(telescopename);
    fprintf('file path for work：%s \n',orgfile);
else
    orgfile = strcat(path1, '\telescopetemplate\', name, '.zmx');
    fprintf('file path origin：%s  \n',orgfile);
end
fprintf('camera：%s \n',Cameralist{camType});


Aperture_list = linspace(0.5, 1.5,11);
seeingfwhm_list = linspace(1.0, 2.0, 6);
magnitude_list = linspace(21, 21, 1);
SkyMag_list = linspace(22.3, 22.3, 1);  
zenith_list = linspace(25, 35, 6);    

FovX_list = linspace(0, 0, 1);  
FovY_list = linspace(0, (temobs.Fov/2), 6);
Exposure_list = linspace(0, 10, 301);

datasize = (size(Aperture_list,2))*(size(seeingfwhm_list,2))*(size(magnitude_list,2))*(size(SkyMag_list,2))*...
           (size(zenith_list,2))*(size(FovX_list,2))*(size(FovY_list,2))*(size(Exposure_list,2));
data_Aperture = zeros(1,datasize);
data_seeingfwhm = zeros(1,datasize);
data_magnitude = zeros(1,datasize);
data_SkyMag = zeros(1,datasize);
data_zenith = zeros(1,datasize);
data_FovX = zeros(1,datasize);
data_FovY = zeros(1,datasize);
data_Exposuretime = zeros(1, datasize);
data_result = zeros(1,datasize);

%% circle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global signal_dataload;
signal_dataload = 0;
if signal_dataload == 1
    data_before = importdata('？？？？？？？？？？？');
    data_Aperture = data_before(1, :);
    index_start = find(data_Aperture==0,1,'first');
    index_start = (index_start-1)/(size(Exposure_list,2));
    disp(index_start);

    data_seeingfwhm = data_before(2,:);
    data_magnitude = data_before(3,:);
    data_SkyMag =  data_before(4,:);
    data_zenith = data_before(5,:);
    data_FovX = data_before(6,:);
    data_FovY = data_before(7,:);
    data_Exposuretime = data_before(6,:);
    data_result = data_before(9,:);
    disp('The data has been processed')
end


hWaitBar = waitbar(0, 'calculating...');
totalIterations = length(Aperture_list) * ...
                length(seeingfwhm_list) * ...
                length(magnitude_list) * ...
                length(SkyMag_list) * ...
                length(zenith_list) * ...
                length(FovX_list) * ...
                length(FovY_list);
completedIterations = 0;

global num
global circle
num = 1;
circle = 1;
for Aperture = Aperture_list
    for seeingfwhm = seeingfwhm_list
        for magnitude_ori = magnitude_list
                for zenith = zenith_list
                        for FovY = FovY_list

                            if signal_dataload == 1
                                if circle <= index_start
                                    num = circle*(size(Exposure_list,2));
                                    circle = circle + 1;
                                    num = num + 1;
                                    continue
                                else   
                                    signal_dataload = 0;
                                    fprintf('circle is %g and num is %g \n', circle, num)
                                end
                            end
                            
                            [SkyMag_,magnitude] = Extinction(SkyMag_ori, magnitude_ori, zenith, wavelength);
                            temobs.SkybackNoise = SkyMag_;
                            [temobs] = temobs.rescal_lens(TheApplication, args, Aperture, workingfolder, ...
                                telescopefolder, telescopename, Fov);

                            temobs.CamType = Cameralist{camType};
                            if strcmpi(temobs.CamType, 'CCD')
                                camtype_num = 1; 
                                temobs.Camsize = 9*1e-6;
                                temobs.CamEff = 0.97;
                                temobs.CamDarkCurrent = 0.0011;
                                temobs.CamReadOutNoise = 1.0;
                                temobs.CamWellDepth = 160000;
                                
                            else
                                camtype_num = 2;
                                temobs.Camsize = 9*1e-6;
                                temobs.CamEff = 0.95;
                                temobs.CamDarkCurrent = 0.0011;
                                temobs.CamReadOutNoise = 2.4;
                                temobs.CamWellDepth = 97000;
                            end

                            [temobs.CamPrice, temobs.CamPixelNumber] = CalCamPriceZOS(TheApplication, args, camtype_num, temobs.Camsize*1000);
                            overalleff = temobs.CamEff*temobs.TelEfficiency;

                            psfmesssage = temobs.getpsf(TheApplication, args, 0.0, FovY, wavelength);
                            data = psfmesssage{1};
                            pixelscale = psfmesssage{2};
                            psfdata = temobs.onerealisticpsf(data, pixelscale, wavelength, temobs.Efocal, seeingfwhm);

                            for Exposure_time = Exposure_list
                                exposuretime = Exposure_time*60;
                                [imagecube,photons] = getimagecube_replace_single(temobs, psfdata,pixelscale,magnitude,wavelength,...
                                                temobs.Camsize,temobs.CamReadOutNoise,temobs.CamDarkCurrent,exposuretime,overalleff,...
                                                stackerror,D_radio);

                                [snr,distance_half,signalphoton,noisephoton] = ...
                                                calsnr11(imagecube{1}{3},imagecube{1}{2},psfdata,pixelscale,temobs.Camsize);

                                data_Aperture(1,num) = Aperture;
                                data_seeingfwhm(1,num) = seeingfwhm;
                                data_magnitude(1,num) = magnitude_ori;
                                data_SkyMag(1,num) = SkyMag_ori;
                                data_zenith(1,num) = zenith;
                                data_FovX(1,num) = 0;
                                data_FovY(1,num) = FovY;
                                data_Exposuretime(1,num) = Exposure_time;
                                data_result(1,num) = snr;

                                num = num + 1;
                            end

                            fprintf('circle is %g / %g \n', circle, totalIterations)

                            waitbar(circle / totalIterations, hWaitBar, ...clear
                                sprintf('wait bar: %.2f%%', (circle / totalIterations) * 100));

%                             if mod(circle, 200) == 0
%                                 resultvector_temp = [data_Aperture;data_seeingfwhm;data_magnitude;data_SkyMag;
%                                                     data_zenith;data_FovX;data_FovY;data_Exposuretime;data_result];
%                                 % 保存结果数组信息
%                                 filename = sprintf('./00_result/fittingDataResult-Tel%d-temp.mat', telType);
%                                 save(filename,'resultvector_temp');
%                                 
%                                 pause(0.02);
%                             end

                            circle = circle + 1;
                            % -------------------------------------------------
                        end 

                end

        end
    end
end
resultvector = [data_Aperture;data_seeingfwhm;data_magnitude;data_SkyMag;data_zenith;data_FovX;data_FovY;data_Exposuretime;data_result];

% close wait bar
close(hWaitBar);
% save
filename = sprintf('./00_result/net_fit_data/fittingDataResult-Tel%d.mat', telType);
save(filename,'resultvector');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% final close ZOS Connect
CleanupConnection(TheApplication);

t_end = clock;
elapsedTime = etime(t_end, t_start);
fprintf('time used: %.2f s\n', elapsedTime);
%% ZOS FUNCTION
function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
% Note -- uncomment the following line to use a custom NetHelper path
% NetHelper = 'C:\Users\shinelong\Documents\Zemax\ZOS-API\Libraries\ZOSAPI_NetHelper.dll';
% This is the path to OpticStudio
NET.addAssembly(NetHelper);

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\OpticStudio\');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% Attempt to create a Standalone connection

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
% This is an issue with how MATLAB interfaces with .NET, and the only
% current workaround is to use 32- or 64-bit versions of both applications.
app = TheConnection.CreateNewApplication();
if isempty(app)
   HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
    app = [];
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MException('zosapi:HandleError', error);
throw(ME);
end

function  CleanupConnection(TheApplication)
% Note - this will close down the connection.

% If you want to keep the application open, you should skip this step
% and store the instance somewhere instead.
TheApplication.CloseApplication();
end






