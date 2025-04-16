fileList = dir(fullfile('images', '*'));
fileList = fileList(~[fileList.isdir]);
files = {fileList.name};
qValues = [10, 25, 50, 75];
psValues = [4,8];

folderName = 'encodings';
if isfolder(folderName)
    rmdir(folderName, 's'); 
    fprintf('Deleted folder: %s\n', folderName);
end
mkdir(folderName);

folderName = 'decodings';
if isfolder(folderName)
    rmdir(folderName, 's'); 
    fprintf('Deleted folder: %s\n', folderName);
end
mkdir(folderName);

for j = 1:length(files)
    figure; % Create a single figure for all patch sizes
    hold on; % Ensure all plots are on the same figure
    grid on;
    
    [~, filename, ~] = fileparts(files{j});

    for m = 1:length(psValues)
        ps = psValues(m)
        rmse_values = zeros(size(qValues));
        bpp_values = zeros(size(qValues));
        
        for k = 1:length(qValues) 
            Q = qValues(k);
            inputfile = files{j};
            

            enc_filename =  sprintf('%s_Q%s_ps%s.%s', filename, num2str(Q),num2str(ps), 'myjpeg');
            disp(enc_filename)

            dec_filename = sprintf('%s_Q%s_ps%s.%s', filename, num2str(Q), num2str(ps),'png');
            disp(dec_filename)

            % call color encode with input q, file
            color_encode(qValues(k), ps, inputfile)
            % call color decode
            color_decode(enc_filename)

            encFile = fullfile('encodings', enc_filename)
            decFile = fullfile('decodings', dec_filename)
            inpFile = fullfile('images', inputfile)

            % RMSE
            orig_img = double(imread(inpFile));
            reconstr_img = double(imread(decFile));
                
            rmse_values(k) = rmse(reconstr_img, orig_img, 'all')
            % BPP
            info = dir(encFile)
            num_pixels = numel(imread(inpFile))         
            bpp_values(k) = (info.bytes * 8) / num_pixels  
            
        end

        plot(bpp_values, rmse_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'DisplayName', sprintf('ps = %d', ps));

        % Label each data point with the corresponding quality factor
        for i = 1:length(qValues)
            text(bpp_values(i), rmse_values(i), ...
                sprintf('Q = %d', qValues(i)), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                'FontSize', 15, 'Color', 'magenta'); 
        end

    end

    % Finalize the single plot
    xlabel('BPP', 'FontSize', 10);
    ylabel('RMSE', 'FontSize', 10);
    title(sprintf("Color JPEG: RMSE vs. BPP for %s", filename), 'FontSize', 12);
    legend('show', 'Location', 'Best');
    set(gca, 'FontSize', 10);

    % Save the figure
    plotfile = sprintf('%s_combined.jpg', filename);
    saveas(gcf, fullfile('color_plots', plotfile));
end
