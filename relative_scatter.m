% Extra data evaluation for the response to the Referee #1 comment 10
% related to the scatter of results for wind-parallel vs wind-perpendicular
% flight segments


load([myprojectpath,filesep,'results.mat']);


vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
        'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
        'slp_psd_UX','slp_psd_VY','slp_psd_W'};

for i_p = 1:numel(planes)
    plane = planes{i_p};
    
    if ismember(plane,{'ATR-EUREC4A','TO-POST'})
        fprintf('%s\n',plane)
        
        MOM = MOM_vec{i_p};
        
        G_std = groupsummary(MOM,{'level','dir2'},{'std'},vars);
        
        G_std = addvars(G_std, groupsummary(MOM,{'level','dir2'},{'mean'},{'alt'}).mean_alt,...
            'After',3,'NewVariableNames','alt');
                
        G_std_sqrtN = G_std;
        G_std_sqrtN{:,5:end} = G_std{:,5:end}./repmat(sqrt(G_std.GroupCount),1,size(G_std,2)-4);
        
        print_table(G_std_sqrtN,{'level','dir2'},["GroupCount","std_"+vars],0,0,3)
        
    end
end