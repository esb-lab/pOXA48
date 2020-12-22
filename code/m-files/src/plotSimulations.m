function plotSimulations(simsPath, which_sims, which_sigmas, which_Ms, which_Cs)




rows=max([length(which_sims) length(which_sigmas) length(which_Ms) length(which_Cs)]);
cols=1;

figure(); clf('reset');
set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');
set(gcf, 'Units','normalized','Position',[0. 0.2 0.6 .3*rows]);


cellno=1;
for i=1:length(which_sims)
    
    if length(which_sims)>1
        subaxis(rows,cols,cellno);
    end
    
    for m=1:length(which_Ms)
        if length(which_Ms)>1
            subaxis(rows,cols,cellno);
        end
        
        for c=1:length(which_Cs)
            if length(which_Cs)>1
                subaxis(rows,cols,cellno);
            end
            
            for s=1:length(which_sigmas)
                
                if length(which_sigmas)>1
                    subaxis(rows,cols,cellno);
                end
                
                imagePath=[simsPath, 'PNG/sim',num2str(sprintf('%0*d',3,which_sims(i))),'_sigma',num2str(round(which_sigmas(s)*100)),'e-2_M',num2str(which_Ms(m)),'_conjRate',num2str(which_Cs(c),'%1.0e'),'.png'];
                disp(imagePath);
                if exist(imagePath,'file')
                    im = imread(imagePath);
                    imshow(im)
                end
                
                cellno=cellno+1;
            end
        end
    end
    
end