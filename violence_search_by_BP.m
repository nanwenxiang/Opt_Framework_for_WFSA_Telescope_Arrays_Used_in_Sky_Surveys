% Use a BP neural network to replace the part of the interactive computation of SNR in the original exposure meter.

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





for telType = 1:1  % telType circle, if there are more than one tel
    temobs = Observer();
    % Hardware defintion
    %workingfolder is defined here
    workingfolder = strcat('workingfolder');
    %folder for telescope samples
    telescopefolder = strcat('telescopetemplate');
    % telescope efficiency
    TelEfficiency = 0.78;
    temobs.TelEfficiency = TelEfficiency;

    SNRLimit = 10;          % SNR limit
    overallfov = 10000;     % field request
    %Overlap area (Defined by Professor Shang 20210326 with 5 col and 5 lin)
    Overlap = 5;
    stackerror = 0.3;
    seeingfwhm = 1.5;     %Seeing
    magnitude_ori = 21;         
    SkyMag_ori = 22.3;
    temobs.SkybackNoise = SkyMag_ori;
    zenith = 30;
    wavelength = 0.55;

    temobs.Wavelength = wavelength;
    [SkyMag,magnitude] = Extinction(SkyMag_ori, magnitude_ori, zenith, wavelength);
    temobs.SkybackNoise = SkyMag;
    width_height = [1,1];   
    
    T_transform = 10;         
    
%     telescopelist = {'MSCH-F2-7deg','PFC-F123-48deg',...
%                 'PFC-F176-42deg',...
%                 'PFC-F28-2deg',...
%                 'PFC-F2-6degree',...
%                 'RC-F3-4deg',...
%                 'RC-F4-3deg',...
%                 'TMA-F3-4deg'};
    telescopelist = {'RC-CORRECT', };
    Cameralist = {'CCD','CMOS'};
%     telType = 2; 
    camType = 1;

%     radio = [0.870, 0.77973, 0.835621, 0.8978, 0.75154, 0.72195, 0.874649, 0.705832];
    radio = [0.795, ];
    D_radio = radio(telType);        

    % full field, not half
%     Fov_list = [7, 4.8, 4.2, 2, 6, 4, 3, 4];
    Fov_list = [2.0, ];
    Fov = Fov_list(telType);     

    FovY_list = linspace(0,(Fov/2),6);
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
    

    theNumOfAperture_down = 0.5;
    theNumOfAperture_up = 1.5;
    theNumOfExposuretime_down = 0;
    theNumOfExposuretime_up = 15;
    SamplingOfEXPOs = 1/(2*60);
    
    
    % mintotalprice 
    global mintotalprice;
    mintotalprice = 9.0e+99;
    
    % circle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    position = {};
    %% load net with file path name
    FileName = ['./00_result/net_fit_data/Tel' num2str(telType) 'fitting_net.mat']; 
    load(FileName, 'net')
    

    for Aperture =  theNumOfAperture_up:-0.05:theNumOfAperture_down     
        fprintf('now, aperture is %g \n', Aperture)   
        [temobs] = temobs.rescal_lens(TheApplication, args, Aperture, workingfolder, ...
        telescopefolder, telescopename, Fov);
        temobs.CamType = Cameralist{camType};
        if strcmpi(temobs.CamType, 'CCD')
            camtype_num = 1;
            fprintf('camera is CCD \n');   
            temobs.Camsize = 9*1e-6;
            temobs.CamEff = 0.97;
            temobs.CamDarkCurrent = 0.0011;
            temobs.CamReadOutNoise = 1.0;
            temobs.CamWellDepth = 160000;
        else
            camtype_num = 2;
            fprintf('camera is CMOS \n');   
            temobs.Camsize = 9*1e-6;
            temobs.CamEff = 0.95;
            temobs.CamDarkCurrent = 0.0011;
            temobs.CamReadOutNoise = 2.4;
            temobs.CamWellDepth = 97000;
        end
        [temobs.CamPrice, temobs.CamPixelNumber] = CalCamPriceZOS(TheApplication, args, camtype_num, temobs.Camsize*1000);
        fprintf('camera price is %g \n', temobs.CamPrice)   
        overalleff = temobs.CamEff*temobs.TelEfficiency;
    
        for Exposuretime = theNumOfExposuretime_down:SamplingOfEXPOs:theNumOfExposuretime_up
            
		    exposuretime = Exposuretime*60;
            fprintf('exposure time is %g s \n', exposuretime) 
    
            datasize = size(FovY_list,2);
            % BP NN ----------------------------------------------------------

            data_Aperture = zeros(1,datasize);
            data_seeingfwhm = zeros(1,datasize);
            data_magnitude = zeros(1,datasize);
            data_SkyMag = zeros(1,datasize);
            data_zenith = zeros(1,datasize);
            data_FovX = zeros(1,datasize);
            data_FovY = zeros(1,datasize);
            data_Exposuretime = zeros(1, datasize);
           
            num = 1;
    
            for FovY = FovY_list
                FovX = 0;
                data_Aperture(1,num) = Aperture;
                data_seeingfwhm(1,num) = seeingfwhm;
                data_magnitude(1,num) = magnitude;
                data_SkyMag(1,num) = SkyMag;
                data_zenith(1,num) = zenith;
                data_FovX(1,num) = FovX;
                data_FovY(1,num) = FovY;
                data_Exposuretime(1,num) = Exposuretime;
    
                num = num + 1;
            end
    
    
            input_train = [data_Aperture;data_seeingfwhm;data_magnitude;data_SkyMag;data_zenith;
                            data_FovX;data_FovY;data_Exposuretime];
    

            snrmat = sim(net, input_train);
            [minsnr,I] = min(snrmat);
    %         fprintf('minimum SNR：%g \n', minsnr);
    
            if minsnr >= SNRLimit
	            stacknumber = (30*60)/(exposuretime+T_transform);
                [Ntel,allover,radio_Field] = temobs.CalNumOfEnough(Overlap,overallfov,width_height);
                temobs.Ntel = Ntel;
                Alltel = ceil(Ntel/stacknumber);
                totalprice = temobs.CalTotalPrice(Alltel,Aperture,exposuretime);

                position{end+1,1} = Aperture;
                position{end,2} = exposuretime;
                position{end,3} = Alltel;
                position{end,4} = totalprice;            
                position{end,5} = allover;
                position{end,6} = radio_Field;
    %             position{end,7} = snrmat;
                position{end,7} = minsnr;

                if totalprice <= mintotalprice   
                    mintotalprice = totalprice;
                    mintelnumber = Alltel;
                    minexposuretime = exposuretime;
                    mintAperture = Aperture;
                    minsnr = min(snrmat);
                    minallover = allover;                        
                    minfacter = [mintotalprice,mintelnumber,minexposuretime,...
                                 mintAperture,minsnr,minallover];
                end

                theNumOfExposuretime_down = Exposuretime;
                break;
            end
    
        end
    end
    
    
    

    % figure
    resultPic = figure();

    title([['Tel：',telescopename],newline,['Cam：',temobs.CamType]]);                    
    xlabel('Ntel');
    zlabel('Exposure Time (s)');
    ylabel({'Aperture (m)'});
    hold on;
    
    allnum = size(position);
    
    prices = cell2mat(position(:, 4)); 
    minPrice = min(prices); 
    maxPrice = max(prices); 
   
    scatter3(cell2mat(position(:, 3)), cell2mat(position(:, 1)), cell2mat(position(:, 2)), ...
        30, prices, 'filled'); 

    c = colorbar;
    c.Label.String = 'Price';
    caxis([minPrice maxPrice]); 
    colormap("autumn"); 

    for K = 1:allnum(1)
        if position{K, 4} == mintotalprice
            str = sprintf('%g', mintotalprice);
            strposition = sprintf('(%g,%g,%g)', position{K, 3}, position{K, 1}, position{K, 2});
            arrow = '\leftarrow';
            text(position{K, 3}, position{K, 1}, position{K, 2}, ...
                {[], [arrow, 'min total price:', str], ['     position:', strposition]}, ...
                'color', 'red ', 'HorizontalAlignment', 'left');
    
    %         arrow = '\rightarrow'; 
    %         text(position{K, 3}, position{K, 1}, position{K, 2}, ...
    %             {[], ['mintotalprice:', str, arrow], ['position:', strposition]}, ...
    %             'color', 'red ', 'HorizontalAlignment', 'right');
        end
    end
    
    % plot line
    for n = 1:allnum(1)-1
        plot3([position{n, 3}, position{n+1, 3}], ...
              [position{n, 1}, position{n+1, 1}], ...
              [position{n, 2}, position{n+1, 2}], 'Color', 'r', 'LineWidth', 1.0); 
    end

    grid on;        
    x = view(3);
    try
	    currentFolder = pwd;   % get the present file path
    
        outfile = sprintf('./00_result/project_BP-Tel%d.fig', telType);  % save as .fig 
        saveas(gcf, outfile);  
        saveas(gcf, sprintf('./00_result/project_BP-Tel%d.png', telType));  % save as PNG 
    catch
	    disp('error')
    end
    hold off;
    
    
    
    
    
    
    
    
    % figure 2
    fig = figure;
    set(fig, 'Position', [50, 50, 625, 450]);
    FontSize = 14;
    prices = cell2mat(position(:, 4)); 
    minPrice = min(prices); 
    maxPrice = max(prices); 
      
    yyaxis left; 
    hold on; 
    h1 = plot(cell2mat(position(:, 1)), cell2mat(position(:, 2)), 'r-', 'LineWidth', 1); % 连接线
    scatter(cell2mat(position(:, 1)), cell2mat(position(:, 2)), 30, prices, 'filled', 'MarkerEdgeColor', 'k'); % 数据点
    colormap("autumn"); 
    caxis([minPrice maxPrice]); 
    xlabel('Aperture (m)', 'FontSize', FontSize);
    ylabel('Exposure Time (s)', 'FontSize', FontSize);
    title([['Tel：',telescopename],newline,['Cam：',temobs.CamType]], 'FontSize', FontSize);
    grid on;
    set(gca, 'YColor', 'red', 'FontSize', 10);
    currentYLimitsLeft = ylim; 
    minY = currentYLimitsLeft(1);
    maxY = currentYLimitsLeft(2);
    midY = (minY + maxY) / 2; 
    yticks([minY, midY, maxY]); 


    yyaxis right; 
    h2 = plot(cell2mat(position(:, 1)), cell2mat(position(:, 3)), 'b-', 'LineWidth', 1);
    scatter(cell2mat(position(:, 1)), cell2mat(position(:, 3)), 30, prices, 'filled', 'MarkerEdgeColor', 'k');
    colormap("autumn"); 
    caxis([minPrice maxPrice]); 
    ylabel('Ntel', 'FontSize', FontSize);
    set(gca, 'YColor', 'blue', 'FontSize', 10); 
    currentYLimitsRight = ylim; 
    minYRight = currentYLimitsRight(1);
    maxYRight = currentYLimitsRight(2);
    midYRight = (minYRight + maxYRight) / 2;
    yticks([minYRight, midYRight, maxYRight]); 

    maxValue = max(cell2mat(position(:, 3)));
    ylim([0, maxValue * 1.1]); 
    

    c2 = colorbar; 
    c2.Label.String = 'Price';
    c2.Label.FontSize = 12; 
    c2.Position = [0.85, 0.15, 0.02, 0.7]; 
    
    [minPriceValue, minIndex] = min(prices);
    minXValue = cell2mat(position(minIndex, 1)); 
    
    str = sprintf('min price = %g', mintotalprice); 
    xline(minXValue, 'k--', 'LineWidth', 1); 
    text(minXValue-0.2, midYRight, str, 'FontSize', FontSize, 'Color', 'k', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'FontSize', 10);
    text(minXValue, minY - 0.1, sprintf('%.2f', minXValue), ...
        'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center');
    
    
    rightYValue = cell2mat(position(minIndex, 3)); 
    yyaxis right; 
    yline(rightYValue, 'b--', 'LineWidth', 1.5);
    text(max(cell2mat(position(:, 1))) + 0.1, rightYValue, sprintf('%.2f', rightYValue), ...
        'Color', 'b', 'FontSize', 13, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    
    leftYValue = cell2mat(position(minIndex, 2)); 
    yyaxis left; 
    yline(leftYValue, 'r--', 'LineWidth', 1.5); 
    text(min(cell2mat(position(:, 1))) - 0.1, leftYValue, sprintf('%.2f', leftYValue), ...
        'Color', 'r', 'FontSize', 13, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    
    hold off; 
    

    legend([h1, h2], {'Exposure Time - Aperture', 'Ntel - Aperture'}, 'Location', 'northeast');

    ax = gca; 
    ax.Position = [0.1, 0.15, 0.65, 0.72]; 

    newFilenameFig = sprintf('./00_result/BP_two_dimensional_%s-Tel%d.fig', Cameralist{camType}, telType);  
    saveas(fig, newFilenameFig);  % save as .fig 
    
    newFilenamePng = sprintf('./00_result/BP_two_dimensional_%s-Tel%d.png', Cameralist{camType}, telType);
    saveas(fig, newFilenamePng);  % save as PNG 
    
    
    
    
    
    
    
    
    
    
    
    
    
    %     figure 3
    %     fig = openfig(outfile, 'reuse');
    %     disp(['User selected ', outfile]);
    %     
    %     set(fig, 'Position', [1252,-5,560,420]);  
    %     set(fig, 'InnerPosition', [1252,-5,560,420]);   
    %     set(fig, 'OuterPosition', [1244,-13,576,514]);  
    %     
    %     
    %     ax = gca;
    %     
    %     dataObjs = findobj(ax, '-property', 'ZData');
    %     if isempty(dataObjs)
    %         error('No 3D data found in the figure.');
    %     end
    %     
    %     xData = get(dataObjs, 'XData');
    %     yData = get(dataObjs, 'YData');
    %     zData = get(dataObjs, 'ZData');
    %     
    %     if ~iscell(xData)
    %         xData = {xData};
    %         yData = {yData};
    %         zData = {zData};
    %     end
    %     
    %     allZData = [];
    %     for k = 1:length(zData)
    %         if ~isempty(zData{k})  
    %             allZData = [allZData; zData{k}(:)];  
    %         end
    %     end
    %     if isempty(allZData)
    %         error('No data found in zData to compute the minimum value.');
    %     end
    %     minZ = 0.0; 
    %     
    %     colors = lines(length(xData));
    %     figure(fig);
    %     hold on;
    %     for k = 1:length(xData)
    %         if k <= 23
    %             plot3(xData{k}, yData{k}, minZ * ones(size(zData{k})), 'Color', [1, 0, 1], 'LineStyle', '--', ...
    %                 'Marker', '.', 'MarkerFaceColor', [1, 0, 1], 'MarkerEdgeColor', [1, 0, 1]);
    %         else
    %             plot3(xData{k}, yData{k}, minZ * ones(size(zData{k})), 'Color', [0, 0, 1], 'LineStyle', '-', ...
    %                 'Marker', '.', 'MarkerFaceColor', [1, 0, 1], 'MarkerEdgeColor', [1, 0, 1]);
    %         end
    %     end
    %     
    %     disp(length(yData));
    %     for k = 1:length(yData)
    %     %     plot3(xData{k}, 1.2 * ones(size(yData{k})), zData{k}, 'Color', colors(k, :), 'LineStyle', '--');
    %         if k <= 23
    %             plot3(xData{k}, theNumOfAperture_up * ones(size(yData{k})), zData{k}, 'Color', [1, 0, 1], 'LineStyle', '-', ...
    %                 'Marker', '.', 'MarkerFaceColor', [1, 0, 1], 'MarkerEdgeColor', [1, 0, 1]);
    %         else
    %             plot3(xData{k}, theNumOfAperture_up * ones(size(yData{k})), zData{k}, 'Color', [0, 0, 1], 'LineStyle', '--', ...
    %                 'Marker', '.', 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', [0, 0, 1] );
    %         end
    %     end
    %     hold off;
    %     
    %     xlabel(ax.XLabel.String);
    %     ylabel(ax.YLabel.String);
    %     zlabel(ax.ZLabel.String);
    %     
    %     set(findobj(ax, 'Type', 'text'), 'FontSize', 14); 
    %     set(ax, 'FontSize', 13);
    %     set(ax.XLabel, 'FontSize', 14); 
    %     set(ax.YLabel, 'FontSize', 14); 
    %     set(ax.ZLabel, 'FontSize', 14); 
    %     title(ax.Title.String, 'FontSize', 14); 
    %     disp('Projection complete.');
    %     newFilename = sprintf('./00_result/proBP_Modified_CCD-Tel%d.fig', telType);  
    %     saveas(fig, newFilename); 
    %     saveas(fig, sprintf('./00_result/proBP_Modified_CCD-Tel%d.png', telType));  

end




    



%% final close ZOS Connect
CleanupConnection(TheApplication);

t_end = clock;
elapsedTime = etime(t_end, t_start);
fprintf('time cost: %.2f s\n', elapsedTime);

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




