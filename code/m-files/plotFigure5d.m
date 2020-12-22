% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Figure 5d in the article
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% PARAMETERS

conj_rates=logspace(-11, -8, 101); %Conjugation rates used to determine optimal 

mus=linspace(0., 0.5, 501);   %DFE mean (cost=1-mu)
sigmas=sqrt([0 0.0075 0.015]);  %DFE standard deviation 

expeDir='../../data/runs/figure5d/';
figDir='../../figures/';

%% Find critical conjugation rates that maintain plasmids (for different mus)

conj_star=zeros(1,length(mus));
for si=1:length(sigmas)
    for mi=1:length(mus)

        matFile=[expeDir,'conjRate_assay_mu',num2str(round(mus(mi)*1e3)),'e-3_sigma',num2str(round(sigmas(si)*1e3)),'e-3.mat'];
        if exist(matFile, 'file')
            disp(['Loading ',matFile]);
            load(matFile);
            if length(conj_rates)>1
                Vq = interp1(expe_pFraction(:,1),conj_rates,.5);
                conj_star(si, mi)=Vq;
            else
                conj_star(si, mi)=Vq;
            end
        else
            conj_star(si, mi)=nan;
        end
    end
end

%% PLOT COST VS CONJUGATION RATE (FOR DIFFERENT SIGMAS)

cmap_sigmas=cbrewer('seq','Reds',length(sigmas)+2);
cmap_sigmas=cmap_sigmas(2:end,:);

figure( );
clf('reset'); 
set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);

leg_sigma={};
pmu=[];
for si =1:length(sigmas)
    
    istar=find(conj_star(si, :)>=0);
    xs=mus(istar);
    ys_star=smooth(conj_star(si,istar));
    if ~isempty(xs)
        %pf = polyfit(xs,ys_star,4);
        %ys=polyval(pf,xs);
        ys=ys_star';
        
        if si==1
            patch([ [0 xs]  fliplr([0 xs]) ], [(([1e-11 ys]))'; 1e-8*ones(length(xs)+1,1)],[1 1 0])  ; hold on;
            leg_sigma{si}=['\sigma^2=',sprintf('%.0f', sigmas(si)^2)];
        else
            leg_sigma{si}=['\sigma^2=',sprintf('%.3f', sigmas(si)^2)];

        end


        patch([xs  fliplr(xs) ], [1e-11*ones(length(xs),1);  fliplr((ys))'],cmap_sigmas(si, :))  
        alpha(0.2)   
        set(gca, 'YScale', 'log')

        plot([xs(1) xs(1)], [1e-11 ys(1)],'Color',cmap_sigmas(si, :),'LineWidth',3); hold on;
        pmu(si)=semilogy(xs,ys, '-','Color',cmap_sigmas(si, :), 'LineWidth',3); hold all;  
    end
    
end
text(.25, 5e-10,'Unstable','FontSize',20)
text(.05, 5e-9,'Stable','FontSize',20)
legend(pmu, leg_sigma,'Location','SouthEast')
xlim([0 max(mus)]);
set(gca,'FontSize',20);
xlabel('Plasmid cost (1-\mu)','FontSize',24);
ylabel('Conjugation rate (\gamma)','FontSize',24);

eval(['export_fig ',figDir,'Figure5d.png']);
