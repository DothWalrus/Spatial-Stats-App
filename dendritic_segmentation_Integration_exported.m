classdef dendritic_segmentation_Integration_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        SelectTIFButton            matlab.ui.control.Button
        tif_field                  matlab.ui.control.EditField
        RunButton                  matlab.ui.control.Button
        NucleiChannelSpinnerLabel  matlab.ui.control.Label
        NucleiChannelSpinner       matlab.ui.control.Spinner
        DCChannelSpinnerLabel      matlab.ui.control.Label
        DCChannelSpinner           matlab.ui.control.Spinner
        SelectCSVButton            matlab.ui.control.Button
        csv_field                  matlab.ui.control.EditField
        debug                      matlab.ui.control.CheckBox
        file_name                  matlab.ui.control.EditField
        ThresholdSliderLabel       matlab.ui.control.Label
        ThresholdSlider            matlab.ui.control.Slider
        CloseFiguresButton         matlab.ui.control.Button
        DisplayTIFChannelsButton   matlab.ui.control.Button
        DiskSizeSpinnerLabel       matlab.ui.control.Label
        DiskSizeSpinner            matlab.ui.control.Spinner
    end

    
    methods (Access = private)
        % reads tif into array
        function FinalImage = tif_reader(app,FileTif)
            InfoImage = imfinfo(FileTif);
            mImage = InfoImage(1).Width;
            nImage = InfoImage(1).Height;
            NumberImages = length(InfoImage);   
            %storage for image layers
            FinalImage1 = zeros(nImage, mImage, NumberImages, 'uint16');
            tifLink = Tiff(FileTif, 'r');
            for i = 1:NumberImages
                tifLink.setDirectory(i);
                FinalImage1(:,:,i) = tifLink.read();
            end
            tifLink.close();
            FinalImage = FinalImage1;
        end
        
        % reads ROI CSV and returns it as array
        function Selection_Image = selection_logical(app,filepath, debug)
            % Read .csv file from ImageJ interface
            Selection = csvread(filepath,1,0);
            X = Selection(:,1);
            Y = Selection(:,2);
            I = Selection(:,3);
            % Remove all selection outside image boundary          
            idx = (X <= 0) | (Y <= 0);
            X(idx) = [];
            Y(idx) = [];
            I(idx) = [];
            I = ones(size(I));
            % Create sparse representation
            Selection_Image = sparse(X,Y,I);    
            % Rotate to true position
            Selection_Image = fliplr(Selection_Image);
            Selection_Image = rot90(Selection_Image);
            % Convert to Logical
            Selection_Image = ~logical(Selection_Image);
            if debug == 1
                figure('Name','ROI Mask'); imagesc(Selection_Image); axis off;
            end
        end
        
        % main function
        function run(app)
            % Nindex - layer for Nuclei Information
            % DCindex - layer for Dendridic Information
            Nindex = app.NucleiChannelSpinner.Value; % nuclei channel
            DCindex = app.DCChannelSpinner.Value; % dendritic cell channel    
            SelectionPath = app.csv_field.Value; % ROI Selecion CSV path
            debug_mode = app.debug.Value; 
            TifPath = string(app.tif_field.Value); % 3-Channel image path
            File = app.file_name.Value; % Generated Name of output file, based on input TIF
            k_threshold = app.ThresholdSlider.Value;
            disk_size = app.DiskSizeSpinner.Value;% convolution disk radius
            % requires ROI to be saved as .csv from Fiji
            try
                Selection = selection_logical(app,SelectionPath, debug_mode);
%                 disp('Selection');
%                 disp(size(Selection));
            catch
                f = uifigure;
                uialert(f, 'Error in selection_logical.', 'ERROR')
                return
            end
            try
                FinalImage = tif_reader(app,TifPath);
%                 disp('FinalImage')
%                 disp(size(FinalImage));
            catch
                f = uifigure;
                uialert(f, 'Error in tif_reader.', 'ERROR')
                return
            end
            try
                % adjust size if ROI doesn't match image size
                Selection = Selection(1:size(FinalImage(:,:,Nindex),1),1:size(FinalImage(:,:,Nindex),2));
            catch
                %padd ones if size doesn't match and error is thrown
                Selection = padarray(Selection',[abs(size(FinalImage(:,:,Nindex),2)-size(Selection,2)) 2],1,'post')';
                Selection = padarray(Selection',[abs(size(FinalImage(:,:,Nindex),1)-size(Selection,1)) 1],1,'post')';
                
                Selection = Selection(1:size(FinalImage(:,:,Nindex),1),1:size(FinalImage(:,:,Nindex),2));
%                 figure;
%                 imagesc(Selection);
            end           
            
            %% Nuclei Centers
            
            % Segmentation of nuclei using DoG filter
            gaussian1 = fspecial('Gaussian',[25,1],5);
            gaussian2 = fspecial('Gaussian',[25,1],6);
            Nuclei = convn(convn((FinalImage(:,:,Nindex))',gaussian1,'same')',gaussian1,'same') - convn(convn((FinalImage(:,:,Nindex))',gaussian2,'same')',gaussian2,'same') ;
            
%             if debug_mode == 1
%                 %   generate a 'in process' image
%                 figure('Name','Magnified Convoluted Nuclie Image');
%                 imagesc(Nuclei(1:400,700:1100)); axis off;
%             end
            
            
            % Finding Centers of Nuclei
            
            Nuclei_Centers = imregionalmax(Nuclei);
            
            %   generate a 'in process' image
            se = strel('disk',3);
            NC_Display = imdilate(Nuclei_Centers,se);
            if debug_mode == 1
                figure('Name','imregionalmax Output'); imagesc(NC_Display); axis off;
            end
            
            % pull off items outside selection area
            Nuclei_Centers(~logical(full(Selection))) = 0;
            
            %   generate a 'in process' image
            se = strel('disk',3);
            NC_Display = imdilate(Nuclei_Centers,se);
            if debug_mode == 1
                figure('Name','imregionalmax Selection Removed'); imagesc(NC_Display); axis off;
            end
            
            % correction for imregionalmax
            f = FinalImage(:,:,Nindex);
            Nuclei_Centers(Nuclei < 0.05*median(f(:)) | f < 0.5*median(f(:))) = 0;
%             if debug_mode == 1
%                 figure('Name','imregionalmax Correction'); imagesc(Nuclei_Centers); axis off;
%             end
            
            %   Results images
            se = strel('disk',3);
            NC_Display = imdilate(Nuclei_Centers,se);
%             if debug_mode == 1
%                 figure('Name','imregionalmax Selection Removed'); imagesc(NC_Display); axis off;
%             end
                            
            if debug_mode == 1
                figure('Name','imregionalmax Corrected'); imagesc(NC_Display); axis off;
                %figure('Name','Nuclei Centers'); imagesc(NC_Display); axis off;
                %figure('Name','Nucei Image'); imagesc(FinalImage(:,:,Nindex)); title("Nuclei Image"); axis off;
                figure('Name','Nuclei Image After DoG Convolution'); imagesc(Nuclei); axis off;
                %figure('Name','Nuclei Center Locations'); imagesc(NC_Display); title("Nuclei Center locations"); axis off;
                %figure('Name','Nuclei Centers Over Nuclei Image'); imagesc(150*NC_Display + 0.1*double(FinalImage(:,:,Nindex))); title("Nuclei over Nuclei Image"); axis off;
                figure('Name','High Contrast Nuclei Centers Over Nuclei Image'); imshowpair(FinalImage(:,:,Nindex),NC_Display); axis off;
            end
            
            %% Filter dendridic cells
            
            %  gaussian3 = fspecial('Gaussian',15,2);
            disk_conv = fspecial('disk',disk_size);
            
            D_thresh = FinalImage(:,:,DCindex);
            D = reshape(D_thresh,1,size(D_thresh,1)*size(D_thresh,2));
            D_dub = double(D);
            D_std = std(D_dub);
            D_thresh(D_thresh < mean(D)+2*D_std) = 0;
            D_thresh(~Selection) = 0;
            
            % figure; imagesc(FinalImage(:,:,DCindex)); 
            % figure; imagesc(D_thresh);
            if debug_mode == 1
               figure('Name','Thresholded DC Channel'); imagesc(D_thresh); axis off;
            end
             
            Dendridic = convn(D_thresh, disk_conv,'same');
            Dendridic = Dendridic/max(max(Dendridic));
             
            if debug_mode == 1
               figure('Name','DC Channel After Disk Convolution'); imagesc(Dendridic); axis off;
            end
            
           % Threshold nuclei centers based on Dendridic Cells
                            
            Nuclei_Dendritic = Nuclei_Centers;
            %D = reshape(Dendridic,1,size(Dendridic,1)*size(Dendridic,2));
            D = Dendridic;
            D(Nuclei_Centers == 0) = 0;
            Nuclei_Dendritic(D<k_threshold) = 0;        
            
            se = strel('disk',3);
            ND_Display = imdilate(Nuclei_Dendritic,se);
            
            if debug_mode == 1
                figure('Name','Nuclei After Thresholding'); imagesc(ND_Display); axis off;
            end
            
            % generate final image
            figure('Name','3 Channel DC Centers'); 
            X = zeros(size(ND_Display,1),size(ND_Display,2),3);
            X(:,:,3) = double(FinalImage(:,:,Nindex))/4095;
            X(:,:,2) = double(FinalImage(:,:,DCindex))/4095;
            X(:,:,1) = double(ND_Display);
            X3 = X(:,:,3);
            X3(bwperim(full(~Selection))) = 1;
            X(:,:,3) = X3;
            X2 = X(:,:,2);
            X2(bwperim(full(~Selection))) = 1;
            X(:,:,2) = X2;
            X1 = X(:,:,1);
            X1(bwperim(full(~Selection))) = 1;
            X(:,:,1) = X1;

            imagesc(X); title(strcat('Nuclei, DC, and DC Centers Threshold=',num2str(k_threshold))); axis off;
            
            %save variables for Stats processing
            k_threshold_name = strrep(num2str(k_threshold,2),'.','_'); %replace decimals with underscore for file acess
            saved_file_name = strcat(File,'thresh',k_threshold_name);
            save(saved_file_name,'Selection','Nuclei_Centers','Nuclei_Dendritic'); %saves ROI, Nuclie centers, and DC nuclie centers as a 3D binary array
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Statistics Section
            
            % save ROI as a list of interior coordinates
            roi_data = Selection; 
            [i,j] = find(roi_data);
            XYROI = [i,j];
            ROIarea = size(XYROI,1);         
            if debug_mode == 1
                figure('Name', 'roi_data');
                imagesc(roi_data);
            end
            
            % compute the perimeter of the ROI; save as sparse array, and M x 2 list of coordinates
            border = sparse(bwperim(roi_data));
            [i,j] = find(border);
            XYborder = [i,j];           
            if debug_mode == 1
                figure('Name', 'border');
                imagesc(border);
            end

            % extract DC Nuclei centers; save as sparse array, and M x 2 list of coordinates
            nuclei_centers = sparse(Nuclei_Dendritic);
            [i,j] = find(nuclei_centers);
            XYnuclei_centers = [i,j];
            numNuclei = size(XYnuclei_centers,1);
            disp("Total number of cells initially: " + numNuclei);
            
            %% Calculate real nuclei statistics
            % compute distances between nuclei and nearest bdry point
            cell_ep_distances = pdist2( XYborder, XYnuclei_centers, 'euclidean', 'Smallest', 1 );
            cell_ep_dist_min = min(cell_ep_distances);
            cell_ep_dist_max = max(cell_ep_distances);
            cell_ep_dist_mean = mean(cell_ep_distances);
            
            % compute distances between nuclei and nuclei
            cell_cell_distances = pdist2( XYnuclei_centers, XYnuclei_centers, 'euclidean', 'Smallest', 2 );
            cell_cell_distances(1,:) = []; % remove self-distances
            cell_cell_dist_min = min(cell_cell_distances);
            cell_cell_dist_max = max(cell_cell_distances);
            cell_cell_dist_mean = mean(cell_cell_distances);
            
            %% Run simulations and get simulation data

            % border and nuclei-nuclei bootstrap
            numSimulations = 1000;
            bootstrap_cell_ep_distances = zeros(numSimulations, numNuclei);
            bootstrap_cell_cell_distances = zeros(numSimulations, numNuclei);
            
            for i = 1:numSimulations
                XYbootstrapcenters = XYROI( randperm(ROIarea,numNuclei), : ); % choose points w/o replacement
                % border stats
                bootstrap_cell_ep_distances(i,:) = pdist2( XYborder, XYbootstrapcenters, 'euclidean', 'Smallest', 1 );
                % nuclei-nuclei stats
                tmp = pdist2( XYbootstrapcenters, XYbootstrapcenters, 'euclidean', 'Smallest', 2 );
                bootstrap_cell_cell_distances(i,:) = tmp(2,:); % save only non-self distances
            end
            if debug_mode == 1
                % could show an artificial DC center image, here, too
            end            
            bootstrap_cell_ep_stats = [min(bootstrap_cell_ep_distances,[],2), max(bootstrap_cell_ep_distances,[],2), mean(bootstrap_cell_ep_distances,2)];
            bootstrap_cell_cell_stats = [min(bootstrap_cell_cell_distances,[],2), max(bootstrap_cell_cell_distances,[],2), mean(bootstrap_cell_cell_distances,2)];
            %% Plot the Real Data and Simulated data in relative count histograms

            % Epithelium distance
            figure('Name', 'DC--epithelium distances');
            histogram(bootstrap_cell_ep_distances(:), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 );  hold on;
            histogram(cell_ep_distances(:), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 );
            title('Distance b/w DC and epithelium')
            ylabel('Density')
            xlabel('Distance to Epithelium (pxl)')
            % set(gca, 'XLim', [0,120] );
            legend('Simulations', 'Real Data');     
            
            % Nearest neighbor distance
            figure('Name', 'DC--DC distances');
            histogram(bootstrap_cell_cell_distances(:), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 ); hold on;
            histogram(cell_cell_distances(:), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 ); 
            title('Distance b/w DC and nearest neighbor')
            ylabel('Density')
            xlabel('Distance to nearest neighbor (pxl)')
            % set(gca, 'XLim', [0,400] );
            legend('Simulations', 'Real Data');
            
            % Mean Epithelium distance
            figure('Name', 'Mean DC--epithelium distances');
            histogram(bootstrap_cell_ep_stats(:,3), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 ); hold on;
            xline(cell_ep_dist_mean); 
            title('Mean distance b/w DC and epithelium')
            ylabel('Density')
            xlabel('Distance to Epithelium (pxl)')
            % set(gca, 'XLim', [0,120] );
            legend('Simulations', 'Real Data');
            
            % Mean Nearest neighbor distance
            figure('Name', 'Mean DC--DC distances');
            histogram(bootstrap_cell_cell_stats(:,3), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 );  hold on;
            xline(cell_cell_dist_mean);
            title('Mean distance b/w DC and nearest neighbor')
            ylabel('Density')
            xlabel('Distance to nearest neighbor (pxl)')
            % set(gca, 'XLim', [0,400] );
            legend('Simulations', 'Real Data');
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diagnostice TXT file output
            DC_count = sum(Nuclei_Dendritic(:)); % Number of identified DC's
            nuc_count = sum(Nuclei_Centers(:)); % Number of identified Nuclei
            params = {'Nuclear Channel','Dendritic Channel', 'Threshold', 'Disk Size','Nuclei Count', 'DC Count', 'Real Max Cell Ep Dist', 'Real Min Cell Ep Dist', 'Real Mean Cell Ep Dist', 'Real Max Cell Cell Dist', 'Real Min Cell Cell Dist', 'Real Mean Cell Cell Dist'};
            param_values = [Nindex, DCindex, k_threshold, disk_size, nuc_count, DC_count, cell_ep_dist_max, cell_ep_dist_min, cell_ep_dist_mean, cell_cell_dist_max, cell_cell_dist_min, cell_cell_dist_mean];
            %disp(DC_count);
            txt_name = strcat(saved_file_name,'.txt');
            txt_id = fopen(txt_name,'w');
            for row = 1:length(params)
                fprintf(txt_id, '%s: %8.3f\n', params{row},param_values(row));
            end
            fclose(txt_id);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SelectTIFButton
        function select_tif(app, event)
            try
                [File, Path] = uigetfile({'*.tif'}, 'File Selector');
                FileTif = strcat(Path,File);
                app.tif_field.Value = FileTif;
                app.file_name.Value = erase(File,'.tif');
            catch
                f = uifigure;
                uialert(f, 'Could not load TIF.', 'ERROR')
                return
            end
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
           run(app)
        end

        % Button pushed function: SelectCSVButton
        function select_csv(app, event)
            [File, Path] = uigetfile({'*.csv'}, 'File Selector');
            Filecsv = strcat(Path,File);
            app.csv_field.Value = Filecsv;
        end

        % Callback function
        function close_figures(app, event)
            close all;
        end

        % Button pushed function: CloseFiguresButton
        function CloseFiguresButtonPushed(app, event)
            close all;
        end

        % Button down function: UIFigure
        function UIFigureButtonDown(app, event)
            
        end

        % Button pushed function: DisplayTIFChannelsButton
        function DisplayTIFChannelsButtonPushed(app, event)

            TifPath = string(app.tif_field.Value);
            Nindex = app.NucleiChannelSpinner.Value;
            DCindex = app.DCChannelSpinner.Value;
            Epindex = [];
            if Nindex == DCindex
                f = uifigure;
                uialert(f, 'Nindex cannot equal DCindex', 'ERROR')
                return
            end
            
            if Nindex + DCindex == 3
                Epindex = 3;
            end
            if Nindex + DCindex == 4
                Epindex = 2;
            end
            if Nindex + DCindex == 5
                Epindex = 1;
            end             
        
            try
                r_chnl = rescale(imread(TifPath,Epindex)); % Red channel
                g_chnl = rescale(imread(TifPath,DCindex)); % Green channel
                b_chnl = rescale(imread(TifPath,Nindex)); % Blue channel
                
                all_blk = zeros(size(r_chnl, 1), size(r_chnl, 2)); % all black channel
                just_red = cat(3, r_chnl, all_blk, all_blk);
                just_green = cat(3, all_blk, g_chnl, all_blk);
                just_blue = cat(3, all_blk, all_blk, b_chnl);
                
%                 figure('Name','TIF Channels');
%                 subplot(1,3,1)
%                 imshow(just_red)
%                 title('Epithelial Cells')
%                 subplot(1,3,2)
%                 imshow(just_green)
%                 title('Dendritic Cells')
%                 subplot(1,3,3)
%                 imshow(just_blue)
%                 title('Cell Nuclei')
                
                figure('Name','Epithelial Cells');
                imshow(just_red);
                figure('Name','Dendritic Cells');
                imshow(just_green);
                figure('Name','Cell Nuclei');
                imshow(just_blue);                        
    
            catch
                f = uifigure;
                uialert(f, 'Could not display TIF channels.', 'ERROR')
                return                
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 562 245];
            app.UIFigure.Name = 'UI Figure';
            app.UIFigure.ButtonDownFcn = createCallbackFcn(app, @UIFigureButtonDown, true);

            % Create SelectTIFButton
            app.SelectTIFButton = uibutton(app.UIFigure, 'push');
            app.SelectTIFButton.ButtonPushedFcn = createCallbackFcn(app, @select_tif, true);
            app.SelectTIFButton.Position = [32 193 99 23];
            app.SelectTIFButton.Text = 'Select TIF';

            % Create tif_field
            app.tif_field = uieditfield(app.UIFigure, 'text');
            app.tif_field.Position = [141 194 100 22];

            % Create RunButton
            app.RunButton = uibutton(app.UIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Position = [32 29 100 22];
            app.RunButton.Text = 'Run';

            % Create NucleiChannelSpinnerLabel
            app.NucleiChannelSpinnerLabel = uilabel(app.UIFigure);
            app.NucleiChannelSpinnerLabel.HorizontalAlignment = 'right';
            app.NucleiChannelSpinnerLabel.Position = [32 108 87 22];
            app.NucleiChannelSpinnerLabel.Text = 'Nuclei Channel';

            % Create NucleiChannelSpinner
            app.NucleiChannelSpinner = uispinner(app.UIFigure);
            app.NucleiChannelSpinner.Limits = [1 3];
            app.NucleiChannelSpinner.Position = [141 108 94 22];
            app.NucleiChannelSpinner.Value = 1;

            % Create DCChannelSpinnerLabel
            app.DCChannelSpinnerLabel = uilabel(app.UIFigure);
            app.DCChannelSpinnerLabel.HorizontalAlignment = 'right';
            app.DCChannelSpinnerLabel.Position = [32 70 71 22];
            app.DCChannelSpinnerLabel.Text = 'DC Channel';

            % Create DCChannelSpinner
            app.DCChannelSpinner = uispinner(app.UIFigure);
            app.DCChannelSpinner.Limits = [1 3];
            app.DCChannelSpinner.Position = [141 70 94 22];
            app.DCChannelSpinner.Value = 2;

            % Create SelectCSVButton
            app.SelectCSVButton = uibutton(app.UIFigure, 'push');
            app.SelectCSVButton.ButtonPushedFcn = createCallbackFcn(app, @select_csv, true);
            app.SelectCSVButton.Position = [32 149 99 23];
            app.SelectCSVButton.Text = 'Select CSV';

            % Create csv_field
            app.csv_field = uieditfield(app.UIFigure, 'text');
            app.csv_field.Position = [141 149 100 22];

            % Create debug
            app.debug = uicheckbox(app.UIFigure);
            app.debug.Text = 'debug';
            app.debug.Position = [303 149 55 22];

            % Create file_name
            app.file_name = uieditfield(app.UIFigure, 'text');
            app.file_name.Editable = 'off';
            app.file_name.Visible = 'off';
            app.file_name.Position = [141 214 100 22];

            % Create ThresholdSliderLabel
            app.ThresholdSliderLabel = uilabel(app.UIFigure);
            app.ThresholdSliderLabel.HorizontalAlignment = 'right';
            app.ThresholdSliderLabel.Position = [300 112 59 22];
            app.ThresholdSliderLabel.Text = 'Threshold';

            % Create ThresholdSlider
            app.ThresholdSlider = uislider(app.UIFigure);
            app.ThresholdSlider.Limits = [0 1];
            app.ThresholdSlider.MajorTicks = [0 0.2 0.4 0.6 0.8 1];
            app.ThresholdSlider.MinorTicks = [0.1 0.3 0.5 0.7 0.9];
            app.ThresholdSlider.Position = [380 121 150 3];
            app.ThresholdSlider.Value = 0.1;

            % Create CloseFiguresButton
            app.CloseFiguresButton = uibutton(app.UIFigure, 'push');
            app.CloseFiguresButton.ButtonPushedFcn = createCallbackFcn(app, @CloseFiguresButtonPushed, true);
            app.CloseFiguresButton.Position = [437 193 100 22];
            app.CloseFiguresButton.Text = 'Close Figures';

            % Create DisplayTIFChannelsButton
            app.DisplayTIFChannelsButton = uibutton(app.UIFigure, 'push');
            app.DisplayTIFChannelsButton.ButtonPushedFcn = createCallbackFcn(app, @DisplayTIFChannelsButtonPushed, true);
            app.DisplayTIFChannelsButton.Position = [300 194 131 22];
            app.DisplayTIFChannelsButton.Text = 'Display TIF Channels';

            % Create DiskSizeSpinnerLabel
            app.DiskSizeSpinnerLabel = uilabel(app.UIFigure);
            app.DiskSizeSpinnerLabel.HorizontalAlignment = 'right';
            app.DiskSizeSpinnerLabel.Position = [303 49 56 22];
            app.DiskSizeSpinnerLabel.Text = 'Disk Size';

            % Create DiskSizeSpinner
            app.DiskSizeSpinner = uispinner(app.UIFigure);
            app.DiskSizeSpinner.Limits = [1 10];
            app.DiskSizeSpinner.Position = [374 49 100 22];
            app.DiskSizeSpinner.Value = 5;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = dendritic_segmentation_Integration_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end