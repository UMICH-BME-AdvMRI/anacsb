function x = ifftdim(x,dims)
    for i = dims
        x = fftshift(ifft(ifftshift(x,i),[],i),i);
    end
end