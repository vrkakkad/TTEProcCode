function disp = correctPhaseWrap(disp,par,options)

orig_disp = disp;
crop_disp = disp(:,:,[1:options.dispEst.dims(1),sum(options.dispEst.dims(1:3))+1:end]);
threshold = 1000*par.lambda/4;
for i = 1:size(crop_disp,1)
    for j = 1:size(crop_disp,2)
        jump = find(abs(diff(crop_disp(i,j,:)))>threshold);
        jump(end+1) = size(crop_disp,3);
        for k=1:length(jump)-1
            try
                offset = crop_disp(i,j,jump(k)+1) - crop_disp(i,j,jump(k)) - mean([crop_disp(i,j,jump(k))-crop_disp(i,j,jump(k)-1) crop_disp(i,j,jump(k)+2)-crop_disp(i,j,jump(k)+1)]);
            catch
                if jump(k)==size(crop_disp,3)-1
                    offset = crop_disp(i,j,jump(k)+1) - crop_disp(i,j,jump(k)) - crop_disp(i,j,jump(k))-crop_disp(i,j,jump(k)-1);
                elseif jump(k)==1
                    offset = crop_disp(i,j,jump(k)+1) - crop_disp(i,j,jump(k)) - crop_disp(i,j,jump(k)+2)-crop_disp(i,j,jump(k)+1);
                end
            end
            crop_disp(i,j,jump(k)+1:jump(k+1)) = crop_disp(i,j,jump(k)+1:jump(k+1)) - offset;
        end
        crop_disp(i,j,:) = crop_disp(i,j,:) - crop_disp(i,j,options.dispEst.ref_idx);
    end
end

disp(:,:,1:options.dispEst.dims(1)) = crop_disp(:,:,1:options.dispEst.dims(1));
disp(:,:,sum(options.dispEst.dims(1:3))+1:end) = crop_disp(:,:,options.dispEst.dims(1)+1:end);
