function encode(Q, ps, src_filename)
    % encode.m
    % This script encodes an image using JPEG compression and saves the encoded data.

    clc;

    if nargin ~= 3
        error(['Incorrect number of input arguments.\n' ...
               'Usage: color_encode(Q, ps, src_filename)\n' ...
               '  Q: Quantization parameter (e.g., integer value).\n' ...
               '  ps: Processing parameter (e.g., scalar or array).\n' ...
               '  src_filename: Source file name (string).\n']);
    end
    
    % Further validation for input types (optional)
    if ~isnumeric(Q) || ~isscalar(Q)
        error('Q must be a numeric scalar. Example: color_encode(10, 5, "image.jpg")');
    end
    if ~isnumeric(ps)
        error('ps must be a numeric value. Example: color_encode(10, [1 2 3], "image.jpg")');
    end
    if ~ischar(src_filename) && ~isstring(src_filename)
        error('src_filename must be a string. Example: color_encode(10, 5, "image.jpg")');
    end


    % Prompt user for input image file
    if ~isfile(['./images/', src_filename])
        error('The specified image file does not exist.');
    end

    im = imread(['./images/', src_filename]);
    figure;
    imshow(im);
    title('Original Image');

    if size(im, 3) == 3
        im = double(rgb2gray(im));
        figure;
        imshow(uint8(im));
        title('Converted to B/W Image');
    else
        im = double(im);
    end

    [H, W] = size(im); % Original image dimensions

    dct_coeffs = im2DCT(im, ps);

    % Set Quality factor
    %Q = 52; % Keep this at a value != 50 to ensure quantization effects

    quantized_dct_coeffs = quantize(dct_coeffs, Q, ps);

    % Huffman encode DC coefficients
    quantized_dc_coeffs = quantized_dct_coeffs(1, 1, :);
    quantized_dc_coeffs = quantized_dc_coeffs(:);
    quantized_dc_diffs = diff(quantized_dc_coeffs);
    huff_dict_dc = generateHuffDict(quantized_dc_diffs, false);
    huff_code_dc = huffmanenco(quantized_dc_diffs, huff_dict_dc); % Vector size: num_patches - 1
    num_dc_coeffs = size(quantized_dc_coeffs, 1);

    % Huffman encode AC coefficients
    zz_quantized_ac_coeffs = zigzag(quantized_dct_coeffs);
    huff_dict_ac = generateHuffDict(zz_quantized_ac_coeffs, true);
    huff_code_ac = runLengthEncode(zz_quantized_ac_coeffs, huff_dict_ac, num_dc_coeffs);

    first_dc_coeff = quantized_dc_coeffs(1);

    % Convert Huffman dictionaries to bitstreams
    ac_bs = dict2Bitstream(huff_dict_ac);
    dc_bs = dict2Bitstream(huff_dict_dc);

    [~, img_name, ~] = fileparts(src_filename);
    output_file_full = saveToDisk(ps, Q, im, num_dc_coeffs,first_dc_coeff, huff_code_dc', huff_code_ac', dc_bs, ac_bs,src_filename, ['./encodings/', img_name, '_Q', num2str(Q), '_ps', num2str(ps), '.myjpeg']);

    disp('Encoding completed successfully.');

    % Nested Functions

    
    function [padded_H, padded_W] = getPaddedDims(im, ps)
        [H,W] = size(im);
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
    
    
    function dct_coeffs = im2DCT(im, ps)
        [H,W] = size(im);
        i = 0;
        % patch size 8 x 8
        numPatches = ceil(H/ps) * ceil(W/ps);
        patches = zeros(ps,ps,numPatches); % 3D array to store patches - size of array is 8 x 8 x N
        dct_coeffs = zeros(ps, ps, numPatches);
        
        [padded_H, padded_W] = getPaddedDims(im, ps);
        % pad the image with 0s (to the end of each dimension)
        % in case that dimension of given image is not divisible by 8
        padded_img = padarray(im, [padded_H - H, padded_W - W], 0, 'post');
        
        for j=1:ps:padded_H - ps + 1
            for k=1:ps:padded_W -ps+1
                patch = padded_img(j:j+ps-1,k:k+ps-1)-128;
                i = i+1;
                %  DCT coefficients
                dct_coeffs(:, :, i) = dct2(patch);
                % store patches if needed else remove it
                patches(:,:,i) = patch;
            end
        end
    end
    
    function quantized_dct_coeffs = quantize(dct_coeffs, Q, ps)
        fifty_mtx = getFiftyMtx(ps);
        Q_mtx = fifty_mtx * (double(50)/ Q);
        quantized_dct_coeffs = round(dct_coeffs./Q_mtx);
    end
    
    function huff_dict = generateHuffDict(vec, remove_zeroes)
        % if remove_zeroes, dont include zeroes in the dict
        if remove_zeroes
            vec = vec(vec ~= 0);
        end
        [unique_els, ~, idx] = unique(vec);       % vec = unique_els(idx)
        counts = accumarray(idx, 1);
        probs = counts / numel(vec);
        [huff_dict,~] = huffmandict(unique_els,probs);
        % unique_els(i) occurs with a probability probs(i)
    end
    
    function zz = zigzag(array)
        % ignore all (1,1,:)
        % zigzag the ac coefficients antidiagonal by antidiagonal
        % return a column vector. this vec has size (ps * ps - 1) * num_patches
    
        % Get the size of the input array
        [ps, ~, num_patches] = size(array);
        
        % Preallocate the result as a column vector
        zz = [];        % target shape: (ps*ps-1) * num_patches
        
        % Loop through each possible sum of indices (i + j)
        for sumIdx = 3:2 * ps
            % Determine the range for the current antidiagonal
            if sumIdx - 1 <= ps
                colStart = sumIdx - 1;
                rowStart = 1;
            else
                colStart = ps;
                rowStart = sumIdx - ps;
            end
            
            % Collect elements along the current antidiagonal
            while colStart >= 1 && rowStart <= ps
                colVec = array(rowStart, colStart, :); colVec = colVec(:);
                zz = [zz; colVec'];
                colStart = colStart - 1;
                rowStart = rowStart + 1;
            end
        end
    
        % shape: [(ps*ps-1), num_patches]
        zz = zz(:);
        % target shape: (ps*ps-1) * num_patches
    end
    
        
    function rle = runLengthEncode(ac_coeffs_vec, huff_dict, num_patches)
        function [size_bits, huff_code] = getHuffmanCode(coeff, huff_dict)
            huff_code = huffmanenco([coeff;], huff_dict);
            size_bits = length(huff_code);
        end
    
        coeffs_per_patch = size(ac_coeffs_vec, 1) / num_patches;
        rle = [];  % Initialize as numeric array
    
        for patch_idx = 1:num_patches
            run_length = 0;
            for i = 1 + (patch_idx-1) * coeffs_per_patch : patch_idx * coeffs_per_patch
                coeff = ac_coeffs_vec(i);
                
                if coeff == 0
                    run_length = run_length + 1;
                    if run_length == 16
                        % Append (15, 0) as 8 bits: [1111, 0000]
                        rle = [rle, [1, 1, 1, 1], zeros(1, 8)];
                        run_length = 1;
                    end
                else
                    [size_bits, huff_code] = getHuffmanCode(coeff, huff_dict);
    
                    % Convert run_length and size_bits to 4-bit binary arrays
                    run_length_bin = dec2bin(run_length, 4) - '0';
                    size_bits_bin = dec2bin(size_bits, 8) - '0';
    
                    % Append run_length, size_bits, and Huffman code
                    rle = [rle, run_length_bin, size_bits_bin, huff_code'];
                    assert(isequal(size(run_length_bin), [1, 4]));
                    assert(isequal(size(size_bits_bin), [1, 8]));
                    assert(isequal(size(huff_code'), [1, size_bits]));
                    assert(bi2de(run_length_bin, 'left-msb') == run_length);
                    assert(bi2de(size_bits_bin, 'left-msb') == size_bits);
                    run_length = 0;
                end
            end
            % Append end-of-block marker [00000000]
            rle = [rle, zeros(1, 12)];
        end
        rle = rle';
    end
        
    
    function bitstream = dict2Bitstream(dict)
        % Convert a Huffman dictionary to a binary bitstream
        % Input:
        %   - huff_dict: A struct where keys are integers, and values are binary arrays
        % Output:
        %   - bitstream: A 1D array of bits (1s and 0s)
        
        bitstream = [];
        % fields = fieldnames(dict);
    
        for i = 1:size(dict, 1)
            key = dict{i, 1};
            huff_code = dict{i, 2};
            assert(iscell(dict));
    
            % Convert the key to binary representation
            if key >= 0
                key_bits = dec2bin(key, 8) - '0';  % Use 8 bits to store the key
            else
                key_bits = [1, dec2bin(-key, 7) - '0'];
            end
    
            % Store the length of the Huffman code as a 8-bit binary number
            code_length = length(huff_code);
            % assert(code_length <= 15 && code_length > 0);
            length_bits = dec2bin(int16(code_length), 8) - '0';
    
            assert(bi2de(length_bits, 'left-msb') == code_length, 'Code length is actually %d, but encoded as %d i.e %s', code_length, bi2de(length_bits, 'left-msb'), mat2str(length_bits));
            assert(bi2de(key_bits(2:end), 'left-msb') == abs(key));
            
            % Append the key, length of the code, and the Huffman code to the bitstream
            bitstream = [bitstream, key_bits, length_bits, huff_code];
        end
    
    end
    
    
    function output_file = saveToDisk(varm3, varm2, im, var1, var2, bs1, bs2, bs3, bs4, src_filename, output_path)
    
        % Helper function to convert bitstream (1s and 0s) to uint8 representation
        function [uint8_array, pad_count] = bitstreamToUint8(bs)
            % Ensure that the bitstream length is a multiple of 8
            len = size(bs, 2);
         
            pad_count = 8 - mod(len, 8);
            bs = [bs, zeros(1, 8 - mod(len, 8))];
            
            % Reshape the bitstream into groups of 8 bits
            bs = reshape(bs, 8, [])';
            
            % Convert each group of 8 bits into a uint8 number
            uint8_array = zeros(1, size(bs, 1), 'uint8');
            for i = 1:size(bs, 1)
                uint8_array(i) = bi2de(bs(i, :), 'left-msb');
            end
        end
        
        % Convert the input bitstreams to a uint8 representation (8 bits)
        [bs1_uint8, pc1] = bitstreamToUint8(bs1);
        [bs2_uint8, pc2] = bitstreamToUint8(bs2);
        [bs3_uint8, pc3] = bitstreamToUint8(bs3);
        [bs4_uint8, pc4] = bitstreamToUint8(bs4);
    
        % Open a file to save the data (the file name will be generated based on current time)
        output_file = output_path;
        fid = fopen(output_file, 'wb');
    
        if fid == -1
            error('Unable to create file for writing.');
        end
        
        % Write the int16 variable
        [varm1, var0] = size(im);
        fwrite(fid, varm3, 'int64');
        fwrite(fid, varm2, 'int64');
        fwrite(fid, varm1, 'int64');%hight
        fwrite(fid, var0, 'int64');    %width
        fwrite(fid, var1, 'int64');
        fwrite(fid, var2, 'int64');
    
        fwrite(fid, size(bs1, 2), 'int64');
        fwrite(fid, pc1, 'int64');
    
        fwrite(fid, size(bs2, 2), 'int64');
        fwrite(fid, pc2, 'int64');
    
        fwrite(fid, size(bs3, 2), 'int64');
        fwrite(fid, pc3, 'int64');
    
        fwrite(fid, size(bs4, 2), 'int64');
        fwrite(fid, pc4, 'int64');
    
        % Write the bitstreams as uint8 arrays
        fwrite(fid, bs1_uint8, 'uint8');
        fwrite(fid, bs2_uint8, 'uint8');
        fwrite(fid, bs3_uint8, 'uint8');
        fwrite(fid, bs4_uint8, 'uint8');
        
        % Close the file
        fclose(fid);
    
        % Display the name of the saved file
        disp(['Data has been saved to: ', output_file]);
        disp(sprintf('The size of the file should be around %d Bytes', size(bs4, 2) / 8 + size(bs1, 2) / 8 + size(bs2, 2) / 8 + size(bs3, 2) / 8));
    
    end
end
