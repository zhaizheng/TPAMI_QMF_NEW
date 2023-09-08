function [clean, noisy] = cryo_em_data(show_or_not)

    filein='img_clean.bin';           
    fid=fopen(filein,'rb');            
    data=fread(fid,'double');   
    fclose(fid);   
    clean = reshape(data,[2000, 64, 64]);
    
    filein2='img_noisy.bin';           
    fid2=fopen(filein2,'rb');            
    data2=fread(fid2,'double');   
    fclose(fid2);   
    noisy = reshape(data2,[2000, 64, 64]);

    data = {clean, noisy};
    if show_or_not == 1
        m = 10;n = 8;
        im = zeros(64*m/2,64*n);
        im_w = [];
        for k = 1:2
            for i = 1:m/2
                for j = 1:n
                    im((i-1)*64+1:i*64,(j-1)*64+1:j*64) = data{k}((i-1)*n+j,:,:);
                end
            end
            im_w = [im_w; im];
        end
        image(im_w*200);
    end
end

