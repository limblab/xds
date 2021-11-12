function [time_frame,curs_p,curs_v,curs_a] = resample_kin(xds, bin_size)

n = bin_size/xds.bin_width;
L = int32(floor(size(xds.EMG, 1)/n));
idx = [];
for i =1:L-1
    idx = [idx, int32(floor(i*n))+1];
end

time_frame = xds.time_frame(idx);
curs_p = xds.kin_p(idx, :);
curs_v = xds.kin_v(idx, :);
curs_a = xds.kin_a(idx, :);

end

