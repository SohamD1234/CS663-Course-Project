function compare_plots()
    fileList = dir(fullfile('images', '*'));
    fileList = fileList(~[fileList.isdir]);
    files = {fileList.name};
    qValues = [10, 25, 50, 75];
    ps = 8;

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
            rmse_values = zeros(size(qValues));
            bpp_values = zeros(size(qValues));
            mat_rmse_values = zeros(size(qValues));
            mat_bpp_values = zeros(size(qValues));
            [~, filename, ~] = fileparts(files{j});
            for k = 1:length(qValues) 
                Q = qValues(k);
                inputfile = files{j};
    
                % matlab jpeg compression
                img = imread(fullfile('images', inputfile));
                compFilename = sprintf('mat_%s_Q%s.%s', filename, num2str(Q), 'jpeg')
                compFile = fullfile('matlab_compressed', compFilename)
                imwrite(img, compFile, 'jpeg', 'Quality', Q)
                % mat_rmse and mat_bpp
                mat_rmse_values(k) = rmse(double(imread(compFile)), double(img), 'all')
                totalPixels = numel(img);
                fileInfo = dir(compFile);
                mat_bpp_values(k) = fileInfo.bytes * 8 / totalPixels;
                % 
    
                enc_filename =  sprintf('%s_Q%s_ps%s.%s', filename, num2str(Q),num2str(ps), 'myjpeg');
                disp(enc_filename)
    
                dec_filename = sprintf('%s_Q%s_ps%s.%s', filename, num2str(Q), num2str(ps),'png');
                disp(dec_filename)
    
                % call encode with input q, file
                encode(qValues(k), ps, inputfile)
                % call decode
                decode(enc_filename)
    
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
    
            figure; 
            xtickformat('%.3f'); 
            plot(bpp_values, rmse_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
            hold on;
            plot(mat_bpp_values, mat_rmse_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
            grid on;
    
            xlabel('BPP', 'FontSize', 4); 
            ylabel('RMSE', 'FontSize', 4); 
            title(sprintf("RMSE vs. BPP ps=%s for %s ",num2str(ps), filename), 'FontSize', 4);
    
            % Label each data point with the corresponding quality factor
            for i = 1:length(qValues)
                text(bpp_values(i), rmse_values(i), ...
                    sprintf('Q = %d', qValues(i)), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                    'FontSize', 7, 'Color', 'magenta'); 
                text(mat_bpp_values(i), mat_rmse_values(i), ...
                    sprintf('Q = %d', qValues(i)), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                    'FontSize', 7, 'Color', 'blue'); 
            end
    
            set(gca, 'FontSize', 7); 
            legend('Grayscale JPEG Compression', 'MAT Compression', 'Location', 'Best');
            plotfile = sprintf('%s_ps%s.jpg',filename,num2str(ps));
            saveas(gcf, fullfile('compare_plots', plotfile)); 
    
        end
    
    
end