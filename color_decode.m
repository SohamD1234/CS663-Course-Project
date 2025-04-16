function color_decode(input_file)
% decode.m
% This script decodes an encoded JPEG file and reconstructs the image.

% Validate input arguments
    if nargin ~= 1
        error(['Incorrect number of input arguments.\n' ...
               'Usage: color_decode(input_file)\n' ...
               '  input_file: File name (string) containing the encoded data.\n']);
    end
    
    % Further validation for input types
    if ~ischar(input_file) && ~isstring(input_file)
        error(['Invalid argument type for input_file.\n' ...
               'Usage: color_decode(input_file)\n' ...
               '  input_file: File name (string) containing the encoded data.\n' ...
               'Example: color_decode("encoded_image.dat")']);
    end

clc;

% Prompt user for input encoded file
if ~isfile(['./encodings/', input_file])
    error('The specified encoded file does not exist.');
end

% Derive output image file name from input file
[~, output_image_file, ~] = fileparts(input_file);

% Read all data sets from disk
data = readFromDisk(['./encodings/', input_file]);
disp('Encoded data has been read from disk.');

% Validate the number of data sets
num_sets = length(data);
disp(['Number of data sets read: ', num2str(num_sets)]);

if num_sets ~= 3
    error('Expected exactly 3 data sets for Y, Cr, Cb channels.');
end
% Initialize a cell array to store the reconstructed channels
reconstructed_channels = cell(3, 1);
channel_names = {'Y', 'Cr', 'Cb'};

% Reconstruct each channel
for heyyy = 1:3
        reconstructed_channels{heyyy} = reconstructChannel(data(heyyy), channel_names{heyyy});
        fprintf('Channel %s reconstructed successfully.\n', channel_names{heyyy});

        figure;
        imshow(reconstructed_channels{heyyy}, []);
        title(['Reconstructed ', channel_names{heyyy}, ' Channel'])
end

% Resize reconstructed channels
if any(cellfun(@isempty, reconstructed_channels))
    error('One or more reconstructed channels are empty.');
end

resized_Y = reconstructed_channels{1};
resized_Cb = imresize(reconstructed_channels{2}, size(resized_Y));
resized_Cr = imresize(reconstructed_channels{3}, size(resized_Y));

% Combine channels into YCrCb image
YCrCb_image = cat(3, resized_Y, resized_Cb, resized_Cr);

disp('Combined Y, Cr, Cb channels into YCrCb image.');

% Convert YCrCb to RGB
try
    RGB_image = ycbcr2rgb(YCrCb_image);
    disp('Converted YCrCb image to RGB.');
catch ME
    error(['Error converting YCrCb to RGB: ', ME.message]);
end

figure;
imshow(RGB_image);
title('Reconstructed RGB Image from YCrCb');

% Save the reconstructed RGB image
try
    imwrite(RGB_image, ['./decodings/', output_image_file, '.png']);
    disp(['Decoded RGB image has been saved to: ', ['./decodings/', output_image_file, '.png']]);
catch ME
    error(['Error saving RGB image: ', ME.message]);
end

% Nested Functions

    function data = readFromDisk(input_file)
        % Helper function to convert uint8 array back to bitstream
        function bs = uint8ToBitstream(uint8_array, pad_count)
            bs = [];
            for i = 1:length(uint8_array)
                bs = [bs, de2bi(uint8_array(i), 8, 'left-msb')];
            end
            if pad_count > 0
                bs = bs(1:end - pad_count);
            end
        end
        
        % Open the file for reading
        fid = fopen(input_file, 'rb');
        if fid == -1
            error('Unable to open the file for reading.');
        end
        
        data = struct('varm3', {}, 'varm2', {}, 'varm1', {}, 'var0', {}, 'var1', {}, 'var2', {}, ...
            'bs1', {}, 'bs2', {}, 'bs3', {}, 'bs4', {});
        
        % Read until end of file
        while ~feof(fid)
            % Read the int64 variables
            temp_varm3 = fread(fid, 1, 'int64');
            if isempty(temp_varm3)
                break; % Reached EOF
            end
            temp_varm2 = fread(fid, 1, 'int64');
            temp_varm1 = fread(fid, 1, 'int64');
            temp_var0  = fread(fid, 1, 'int64');
            temp_var1  = fread(fid, 1, 'int64');
            temp_var2  = fread(fid, 1, 'int64');
            
            % Read the size and padding count for each bitstream
            size_bs1 = fread(fid, 1, 'int64');
            pc1      = fread(fid, 1, 'int64');
            
            size_bs2 = fread(fid, 1, 'int64');
            pc2      = fread(fid, 1, 'int64');
            
            size_bs3 = fread(fid, 1, 'int64');
            pc3      = fread(fid, 1, 'int64');
            
            size_bs4 = fread(fid, 1, 'int64');
            pc4      = fread(fid, 1, 'int64');
            
            % Calculate the number of uint8 bytes to read
            num_bytes_bs1 = floor(size_bs1 / 8)+1;
            num_bytes_bs2 = floor(size_bs2 / 8)+1;
            num_bytes_bs3 = floor(size_bs3 / 8)+1;
            num_bytes_bs4 = floor(size_bs4 / 8)+1;
            
            % Read the bitstreams as uint8 arrays
            bs1_uint8 = fread(fid, num_bytes_bs1, 'uint8');
            bs2_uint8 = fread(fid, num_bytes_bs2, 'uint8');
            bs3_uint8 = fread(fid, num_bytes_bs3, 'uint8');
            bs4_uint8 = fread(fid, num_bytes_bs4, 'uint8');
            
            % Handle cases where fread might return fewer bytes
            if length(bs1_uint8) < num_bytes_bs1 || ...
                    length(bs2_uint8) < num_bytes_bs2 || ...
                    length(bs3_uint8) < num_bytes_bs3 || ...
                    length(bs4_uint8) < num_bytes_bs4
                warning('Unexpected end of file encountered.');
                break;
            end
            
            % Convert uint8 arrays back to bitstreams
            temp_bs1 = uint8ToBitstream(bs1_uint8, pc1);
            temp_bs2 = uint8ToBitstream(bs2_uint8, pc2);
            temp_bs3 = uint8ToBitstream(bs3_uint8, pc3);
            temp_bs4 = uint8ToBitstream(bs4_uint8, pc4);
            
            % Create a struct for the current set
            current_set.varm3 = temp_varm3;
            current_set.varm2 = temp_varm2;
            current_set.varm1 = temp_varm1;
            current_set.var0  = temp_var0;
            current_set.var1  = temp_var1;
            current_set.var2  = temp_var2;
            current_set.bs1   = temp_bs1;
            current_set.bs2   = temp_bs2;
            current_set.bs3   = temp_bs3;
            current_set.bs4   = temp_bs4;
            
            % Append to the data struct array
            data(end+1) = current_set;
        end
        
        % Close the file
        fclose(fid);
    end

    function dict = bitstream2Dict(bitstream)
        % Convert a binary bitstream back to a Huffman dictionary
        % Input:
        %   - bitstream: A 1D array of bits (1s and 0s)
        % Output:
        %   - huff_dict: A struct where keys are integers, and values are binary arrays
        
        dict = {};
        idx = 1;
        while idx <= length(bitstream)
            % Extract the 12-bit key
            key_bits = bitstream(idx:idx+11);
            key = bi2de(key_bits(2:end), 'left-msb');
            if (key_bits(1) == 1)
                key = - key;
            end
            idx = idx + 12;
            
            % Extract the 8-bit code length
            length_bits = bitstream(idx:idx+7);
            code_length = bi2de(length_bits, 'left-msb');
            idx = idx + 8;
            
            % Extract the Huffman code of the given length
            huff_code = bitstream(idx:idx+code_length-1);
            idx = idx + code_length;
            
            % Add to the Huffman dictionary
            huff_code = double(huff_code);
            dict(end+1, :) = {[key], double(huff_code)};
        end
    end



    function ac_coeffs_vec = runLengthDecode(rle, huff_dict, num_patches, coeffs_per_patch)
        % Decode the run-length encoded AC coefficients using the Huffman dictionary
        ac_coeffs_vec = zeros(num_patches * coeffs_per_patch, 1);
        out_idx = 1;
        in_idx = 1;
        
        while in_idx <= length(rle)
            % Ensure enough bits for run_length and size_bits
            if in_idx + 7 > length(rle)
                warning('Insufficient bits for run_length and size_bits.');
                break;
            end
            
            % Extract run_length (4 bits) and size_bits (4 bits)
            run_length = bi2de(rle(in_idx:in_idx+3)', 'left-msb');
            size_bits = bi2de(rle(in_idx+4:in_idx+11)', 'left-msb');
            in_idx = in_idx + 12;
            
            if size_bits == 0
                if run_length == 15
                    out_idx = out_idx + 15;
                else
                    % Update out_idx to next multiple of coeffs_per_patch
                    new_idx = out_idx - mod(out_idx, coeffs_per_patch) + 1 - 2 * coeffs_per_patch;
                    while new_idx < out_idx
                        new_idx = new_idx + coeffs_per_patch;
                    end
                    out_idx = new_idx;
                end
                continue;
            end
            
            % Update out_idx based on run_length
            out_idx = out_idx + run_length;
            
            % Ensure enough bits for the Huffman code
            if in_idx + size_bits - 1 > length(rle)
                warning('Insufficient bits for Huffman code.');
                break;
            end
            
            % Extract the Huffman code
            huff_code = rle(in_idx:in_idx+size_bits-1);
            in_idx = in_idx + size_bits;
            
            % Decode the Huffman code
            try
                coeff = huffmandeco(huff_code, huff_dict);
            catch
                warning('Huffman decoding failed for a code.');
                coeff = 0; % Assign a default value or handle as needed
            end
            
            if isempty(coeff)
                warning('Decoded coefficient is empty. Assigning zero.');
                coeff = 0;
            end
            
            % Assign the coefficient
            if out_idx <= length(ac_coeffs_vec)
                ac_coeffs_vec(out_idx) = coeff;
                out_idx = out_idx + 1;
            else
                warning('AC coefficients vector exceeded expected length.');
                break;
            end
        end
        
        % Trim any excess zeros
        % ac_coeffs_vec = ac_coeffs_vec(1:out_idx-1);
    end

    function dequantized_dct_coeffs = dequantize(quantized_dct_coeffs, Q, ps)
        fifty_mtx = getFiftyMtx(ps);
        Q_mtx = fifty_mtx * (50 / Q);
        dequantized_dct_coeffs = quantized_dct_coeffs .* Q_mtx;
    end



    function reconstructed_img = DCT2im(dct_coeffs, ps, H, W, channel_name)
        % Perform inverse DCT on each block and reconstruct the image
        [padded_H, padded_W] = getPaddedDims(H, W, ps);
        reconstructed_img = zeros(padded_H, padded_W);
        
        % Reset patch index
        num_patches = size(dct_coeffs, 3);
        disp(['Number of patches in DCT coefficients: ', num2str(num_patches)]);
        
        for i = 1:num_patches
            % Inverse DCT for each block
            block = idct2(dct_coeffs(:,:,i));
            
            % Shift chrominance channels by 128
            block = block + 128;
            
            % Determine the block's position in the image
            row = ceil(i / ceil(padded_W / ps));
            col = mod(i-1, ceil(padded_W / ps)) + 1;
            
            % Place the reconstructed block in the image
            reconstructed_img((row-1)*ps + 1:row*ps, (col-1)*ps + 1:col*ps) = block;
        end
        
        % Crop the image to original dimensions
        reconstructed_img = reconstructed_img(1:H, 1:W);
    end

    function uzz = unzigzag(vec, num_patches, ps)
        % unzigzag the ac coefficients antidiagonal by antidiagonal and create
        % a temp array
        % put down dummy values at (1,1,:)
        % reshape the temp array and set the dc_idcs (starting) to []
        % return a 2d array of size (ps * ps - 1), num_patches
        vec = reshape(vec, [ps*ps-1, num_patches]);
        uzz = zeros([ps, ps, num_patches]);
        itr = 1;
        
        % Fill in the rest of the coefficients in antidiagonal order
        for sumIdx = 3:(2 * ps)
            if sumIdx <= ps + 1
                rowStart = 1;
                colStart = sumIdx-1;
            else
                colStart = ps;
                rowStart = sumIdx - ps;
            end
            
            while colStart >= 1 && rowStart <= ps
                uzz(rowStart, colStart, :) = vec(itr, :);
                colStart = colStart - 1;
                rowStart = rowStart + 1;
                itr = itr+1;
            end
        end
        uzz = reshape(uzz, [ps^2, num_patches]);
        uzz(1, :) = [];
    end

    function resized_fifty_mtx = getFiftyMtx(ps)
    % Define the standard 8x8 quantization matrix for quality factor = 50
        fifty_mtx = [16 11 10 16 24 40 51 61;
           12 12 14 19 26 58 60 55;
           14 13 16 24 40 57 69 56;
           14 17 22 29 51 87 80 62;
           18 22 37 56 68 109 103 77;
           24 35 55 64 81 104 113 92;
           49 64 78 87 103 121 120 101;
           72 92 95 98 112 100 103 99];
        resized_fifty_mtx = imresize(fifty_mtx, [ps, ps], 'bilinear');
    end

    function [padded_H, padded_W] = getPaddedDims(H, W, ps)
        pad_h = mod(-H, ps);
        pad_w = mod(-W, ps);
        padded_H = H + pad_h;
        padded_W = W + pad_w;
    end

    function reconstructed_img = reconstructChannel(current_data, channel_name)
        disp(['--- Reconstructing Channel: ', channel_name, ' ---']);
        
        % Extract variables for the current set
        ps = current_data.varm3;
        Q = current_data.varm2;
        H = current_data.varm1;
        W = current_data.var0;
        num_dc_coeffs = current_data.var1;
        first_dc_coeff_disk = current_data.var2;
        
        % Extract bitstreams for the current set
        huff_code_dc_disk = current_data.bs1';
        huff_code_ac_disk = current_data.bs2';
        dc_bs_disk = current_data.bs3;
        ac_bs_disk = current_data.bs4;
        
        % Convert bitstreams to strings
        % huff_code_dc_disk = num2str(huff_code_dc_disk');
        % huff_code_dc_disk = huff_code_dc_disk(~isspace(huff_code_dc_disk));
        
        % huff_code_ac_disk = num2str(huff_code_ac_disk');
        % huff_code_ac_disk = huff_code_ac_disk(~isspace(huff_code_ac_disk));
        
        % Reconstruct Huffman dictionaries
        huff_dict_ac_recovered = bitstream2Dict(ac_bs_disk);
        huff_dict_dc_recovered = bitstream2Dict(dc_bs_disk);
        
        % Reconstruct AC coefficients
        zz_dehuff_ac_coeffs = runLengthDecode(huff_code_ac_disk, huff_dict_ac_recovered, num_dc_coeffs, ps * ps - 1);
        dehuff_ac_coeffs = unzigzag(zz_dehuff_ac_coeffs, num_dc_coeffs, ps);
        disp(['Reconstructed AC coefficients for channel ', channel_name, ': ', mat2str(dehuff_ac_coeffs(:,1)')]);
        
        % Reconstruct DC coefficients
        try
            dehuff_dc_diffs = [huffmandeco(huff_code_dc_disk, huff_dict_dc_recovered)];
            dehuff_dc_coeffs = cumsum([first_dc_coeff_disk; dehuff_dc_diffs]);
            disp(['Reconstructed DC coefficients: ', mat2str(dehuff_dc_coeffs')]);
        catch ME
            error(['Huffman decoding failed for DC coefficients in channel ', channel_name, ': ', ME.message]);
        end
        
        
        
        % Interweave DC and AC coefficients
        try
            dehuff_dct_coeffs = [dehuff_dc_coeffs'; dehuff_ac_coeffs];
            dehuff_dct_coeffs = reshape(dehuff_dct_coeffs, [ps, ps, num_dc_coeffs]);
        catch
            error(['Error reshaping DCT coefficients for channel ', channel_name, '.']);
        end
        
        
        % Dequantize DCT coefficients
        dequantized_dct_coeffs = dequantize(dehuff_dct_coeffs, Q, ps);
        
        % Reconstruct the image from DCT coefficients
        try
            reconstructed_img = DCT2im(dequantized_dct_coeffs, ps, H, W,channel_name);
        catch ME
            error(['Error during inverse DCT for channel ', channel_name, ': ', ME.message]);
        end
        
        % Convert to uint8
        reconstructed_img = uint8(reconstructed_img);
    end
end