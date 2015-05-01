function  TrtEffvsDurPlot(params,plot_handle)

    %Initiate symbolics variables
    syms dTc dVc dTh dVh Tcstar Vcstar Thstar Vhstar 

    %Write equations as strings

    dTc = 'sc + r1*Tc*(1-(Tc +Ic)/Tcmax) - dc*Tc -beta_c*Tc*Vc';
    dIc = 'beta_c*Tc*Vc + r2*Ic*(1-(Tc +Ic)/Tcmax)-delta_c*(1+alph*Th)*Ic';
    dVc = 'p*Ic -c*Vc';
    dTh = 'sh*(1+gamm*Ic) - dh*Th - beta_h*VL*Th';

    %Solve for equilibria
    [Tcstar, Icstar,  Vcstar, Thstar] = solve(dTc, dIc, dVc, dTh,  'Tc','Ic','Vc','Th');

    % Assign equilibrium index
    equil_index = 1;
    
    %calculate adjusted delta
    params.delta_c = adjust_delta(params.delta_c_orig,Thstar, params, equil_index);

    % Set up HIV viral loads to use
    HIV_VLs = [0, 8e3, 2.5e4, 1e5,  1e6];
    
    % Set up Treatment efficacies to use
    eff = [.6 .65 .7 .7725 .85 .95 .97 .99 1];
    
    % Set up treatment durations to use
    dur = 7*[12:6:60];
    
    % pre allocate matrix for indication of SVR
    cleared_vec = zeros(length(dur),length(HIV_VLs));
    
    % cycle through viral loads, durations, and efficacies to pull out
    % which combinations result in SVR
    for k = 1:length(HIV_VLs)
        params.VL = HIV_VLs(k);
        clearedmat = zeros(length(eff),length(dur));
        for j = 1:length(dur)
            for i = 1:length(eff)
                while i <length(eff)+1
                    % Run function to test if patient clears with efficacy
                    % i for duration j
                    clearedmat(i,j) = clear_test(eff(i),dur(j), Tcstar, Icstar, Vcstar, Thstar,equil_index,params);
                    % if clearance is achieved, move to next duration
                    if clearedmat(i,j) ==1
                        i = length(eff)+1;
                    else
                        i= i+1;
                    end
                end
                %find index of minimum efficacy that leads to clearance
                eff_index = min(find(clearedmat(:,j)==1));
                if ~isempty(eff_index)
                    cleared_vec(j,k) = eff(eff_index);
                else
                    cleared_vec(j,k) = 1.1;
                end
            end
        end    
    end

    %Calculate CD4s for legend
    CD4s = zeros(size(HIV_VLs));
    paramnames = fieldnames(params);

    for i = 1:numel(paramnames)
        eval(char(strcat(paramnames(i),'= params.',paramnames{i},';')))
    end
    for j = 1:length(HIV_VLs)
        VL = HIV_VLs(j);
        CD4s(j) = eval(char(Thstar(equil_index)));
        CD4s(j) = round(CD4s(j));
    end
    
    %Plot combinations
    %set up y axis
    yaxis_vec = [0:.1:1];

    hold on
    set(gca,'ColorOrder', bone(7))
    %plot duration vs efficacy on input axes
    plot(plot_handle,dur/7,cleared_vec,'LineWidth',3)
    % make y axis into percentages
    pc_vec = repmat('%',length(yaxis_vec),1);
    yticklabel_pc = [num2str(yaxis_vec'*100) pc_vec];
    set(gca, 'YLim',[0 1],'XLim',[min(dur/7) max(dur/7)])
    set(gca,'XTick',dur/7,'YTick', yaxis_vec, 'YTickLabel', yticklabel_pc,'FontSize',20)
    ylabel('Treatment Efficacy Needed for SVR', 'FontSize',20)
    xlabel('Treatment Duration (weeks)','FontSize',20)
    legend(['CD4 = ', num2str(CD4s(1))],['CD4 = ', num2str(CD4s(2))],['CD4 = ', num2str(CD4s(3))],['CD4 = ', num2str(CD4s(4))],['CD4 = ', num2str(CD4s(5))],'Location','SouthEast')
    text(18,.5, 'No SVR','FontSize',20,'FontWeight', 'Normal','Parent',plot_handle)
    text(55,1.01, 'SVR','FontSize',20,'FontWeight', 'Normal','Parent',plot_handle)

