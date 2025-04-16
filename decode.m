function decode(input_file)
    % decode.m
    % This script decodes an encoded JPEG file and reconstructs the image.

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
    if ~isfile(['./encodings/', input_file])
        error('The specified encoded file does not exist.');
    end

    % Prompt user for output image file name
    % output_image_file = input('What to name the decoded image?: ', 's');
    [~, output_image_file, ~] = fileparts(input_file);

    [ps, Q, H, W, num_dc_coeffs, first_dc_coeff_disk, huff_code_dc_disk, huff_code_ac_disk, dc_bs_disk, ac_bs_disk] = readFromDisk(['./encodings/', input_file]);
    disp('Encoded data has been read from disk.');

    % Reconstruct Huffman dictionaries


    huff_code_dc_disk = huff_code_dc_disk';
    huff_code_ac_disk = huff_code_ac_disk';
    huff_dict_ac_recovered = bitstream2Dict(ac_bs_disk);
    huff_dict_dc_recovered = bitstream2Dict(dc_bs_disk);

    % Reconstruct DC coefficients
    dehuff_dc_diffs = [huffmandeco(huff_code_dc_disk, huff_dict_dc_recovered)];
    dehuff_dc_coeffs = cumsum([first_dc_coeff_disk; dehuff_dc_diffs]); % Vector size: num_patches

    % Reconstruct AC coefficients
    zz_dehuff_ac_coeffs = runLengthDecode(huff_code_ac_disk, huff_dict_ac_recovered, num_dc_coeffs, ps * ps - 1);
    dehuff_ac_coeffs = unzigzag(zz_dehuff_ac_coeffs, num_dc_coeffs);

    % Interweave DC and AC coefficients
    dehuff_dct_coeffs = [dehuff_dc_coeffs'; dehuff_ac_coeffs];
    dehuff_dct_coeffs = reshape(dehuff_dct_coeffs, [ps, ps, num_dc_coeffs]);

    dequantized_dct_coeffs = dequantize(dehuff_dct_coeffs, Q, ps);

    reconstructed_img = DCT2im(dequantized_dct_coeffs, ps, H, W);

    % Display the reconstructed image
    reconstructed_img = reconstructed_img(1:H, 1:W);
    figure;
    imshow(uint8(reconstructed_img));
    title('Reconstructed Image');

    % Save the reconstructed image
    imwrite(uint8(reconstructed_img), ['./decodings/', output_image_file, '.png']);
    disp(['Decoded image has been saved to: ', ['./decodings/', output_image_file, '.png']]);

    % Nested Functions

    function [varm3, varm2, varm1, var0, var1, var2, bs1, bs2, bs3, bs4] = readFromDisk(input_file)

    % Helper function to convert uint8 array back to bitstream
    function bs = uint8ToBitstream(uint8_array, pad_count)
        % Convert each uint8 back into a 8-bit binary array
        bs = [];
        for i = 1:length(uint8_array)
            % Convert each uint8 to an 8-bit binary representation
            bs = [bs, de2bi(uint8_array(i), 8, 'left-msb')];
        end
        
        % Remove padding
        if pad_count > 0
            bs = bs(1:end - pad_count);
        end
    end

    % Open the file for reading
    fid = fopen(input_file, 'rb');
    
    if fid == -1
        error('Unable to open the file for reading.');
    end
    
    % Read the int16 variable
    varm3 = fread(fid, 1, 'int64');
    varm2 = fread(fid, 1, 'int64');
    varm1 = fread(fid, 1, 'int64');

    var0 = fread(fid, 1, 'int64');
    var1 = fread(fid, 1, 'int64');
    var2 = fread(fid, 1, 'int64');
    
    % Read the size and padding count for each bitstream
    size_bs1 = fread(fid, 1, 'int64');
    pc1 = fread(fid, 1, 'int64');
    
    size_bs2 = fread(fid, 1, 'int64');
    pc2 = fread(fid, 1, 'int64');
    
    size_bs3 = fread(fid, 1, 'int64');
    pc3 = fread(fid, 1, 'int64');
    
    size_bs4 = fread(fid, 1, 'int64');
    pc4 = fread(fid, 1, 'int64');
    
    % Read the bitstreams as uint8 arrays
    bs1_uint8 = fread(fid, floor(size_bs1 / 8)  + 1, 'uint8');
    bs2_uint8 = fread(fid, floor(size_bs2 / 8)  + 1, 'uint8');
    bs3_uint8 = fread(fid, floor(size_bs3 / 8)  + 1, 'uint8');
    bs4_uint8 = fread(fid, floor(size_bs4 / 8)  + 1, 'uint8');
    
    % Convert uint8 arrays back to bitstreams
    bs1 = uint8ToBitstream(bs1_uint8, pc1);
    bs2 = uint8ToBitstream(bs2_uint8, pc2);
    bs3 = uint8ToBitstream(bs3_uint8, pc3);
    bs4 = uint8ToBitstream(bs4_uint8, pc4);
    
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
        % Extract the 8-bit key
        key_bits = bitstream(idx:idx+7);
        key = bi2de(key_bits(2:end), 'left-msb');
        if (key_bits(1) == 1)
            key = - key;
        end
        idx = idx + 8;

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

function dequantized_dct_coeffs = dequantize(quantized_dct_coeffs, Q, ps)
    fifty_mtx = getFiftyMtx(ps);
    Q_mtx = fifty_mtx * (double(50)/ Q);
    dequantized_dct_coeffs = quantized_dct_coeffs.*Q_mtx;
end
    function reconstructed_img = DCT2im(dct_coeffs, ps, H, W)
    % Initialize an empty matrix for the reconstructed image
    reconstructed_img = zeros(getPaddedDims(H, W, ps));

    [padded_H, padded_W] = getPaddedDims(H, W, ps);
    
    % Reset patch index
    i = 0;
    
    % Iterate over each patch position in the image and place reconstructed patch
    for j = 1:ps:padded_H - ps + 1
        for k = 1:ps:padded_W - ps + 1
            i = i + 1;
    
            % Inverse DCT on the coefficients to reconstruct the patch
            reconstructed_patch = idct2(dct_coeffs(:,:,i));
            
            % Place the reconstructed patch back into the image
            reconstructed_img(j:j+ps-1, k:k+ps-1) = reconstructed_patch+128;
        end
    end
end


function uzz = unzigzag(vec, num_patches)
    % unzigzag the ac coefficients antidiagonal by antidiagonal and create
    % a temp array
    % put down dummy values at (1,1,:)
    % reshape the temp array and set the dc_idcs (starting) to []
    % return a 2d array of size (ps * ps - 1), num_patches
    ps = sqrt(size(vec, 1) / num_patches + 1);
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

    function ac_coeffs_vec = runLengthDecode(rle, huff_dict, num_patches, coeffs_per_patch)
    ac_coeffs_vec = zeros(num_patches * coeffs_per_patch, 1);
    out_idx = 1;
    in_idx = 1;

    while in_idx <= length(rle)
        % Extract run_length (4 bits) and size_bits (4 bits)
        run_length = bi2de(rle(in_idx:in_idx+3)', 'left-msb');
        size_bits = bi2de(rle(in_idx+4:in_idx+11)', 'left-msb');
        in_idx = in_idx + 12;

        % Handle (15,0) case for 15 consecutive zeros
        if size_bits == 0
            if run_length == 15
                out_idx = out_idx + 15;
            else
                % Update out_idx to next multiple of coeffs_per_patch
                new_idx = out_idx - mod(out_idx, coeffs_per_patch) + 1 - 2*coeffs_per_patch;
                while new_idx < out_idx
                    new_idx = new_idx + coeffs_per_patch;
                end
                out_idx = new_idx;
            end
            continue;
        end
        out_idx = out_idx + run_length;
        % Extract the Huffman code based on size_bits
        huff_code = rle(in_idx:in_idx+size_bits-1);
        coeff = huffmandeco(huff_code, huff_dict);
        assert(~isempty(huff_code));
        ac_coeffs_vec(out_idx) = coeff;
        in_idx = in_idx + size_bits;
        out_idx = out_idx + 1;
    end
end

    
function [padded_H, padded_W] = getPaddedDims(H, W, ps)
    pad_h = mod(-H, ps);
    pad_w = mod(-W, ps);
    padded_H = H+pad_h;
    padded_W = W+pad_w;
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
end
