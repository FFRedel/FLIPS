function x_output = patch2image(x_input,im,store_inp)
% Reconstructing image from the extracted/reconstructed (sliding) image patches

[imdim1,imdim2]=size(im);
indices_store_patches_inp= reshape(1:imdim1*imdim2,[imdim1 imdim2]);

rebuild_inp_rec     =  zeros(imdim1,imdim2);
for i=1:imdim1
    for j=1:imdim2
        idx_inp = store_inp==indices_store_patches_inp(i,j);
        rebuild_inp_rec(i,j)=mean(x_input(idx_inp),'all');
    end
end
x_output = rebuild_inp_rec;

end