function  genSecrecyResultsNew()
clc;close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        params.M = 1;                                 
        params.P = 1; 
        params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.la_u =   600e-6 ;               % users density (users/m2)
        params.Ps = 100e-3;
        params.BW = 20e6;
        params.No = 10^(-17.4) * params.BW;
         params.channel = 'Rice';
         params.rice.Km = 32;    
         params.rice.Ke = 0;

sim = 'Noise' ;        
%sim = 'LeakageAnalytical';
%sim = 'UserDensity';
%sim = 'EvesDensity';
%sim = 'CellDensity';
%sim = 'test'
switch (sim)
    case 'Noise'
         params.la_s =  [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1] ; 
         params.channel = 'Rayleigh';
         params.la_e =   500e-6 ;
         params.alpha = 4;
         [Rm_math_noise , Re_math_noise  , Rs_math_noise] = com_secrecy_results_with_smallcell_density_noise(params);
         [Rm_math , Re_math  , Rs_math] = com_secrecy_results_with_smallcell_density(params);
         
         f1 = semilogx(params.la_s ,Rs_math_noise,'k-' ,params.la_s,Rs_math,'ko');
            set(f1,'MarkerSize',15);
            set(f1,'LineWidth',4);
            legend({'Exact' , 'Approximate'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
    case 'test'
         params.la_e =   60e-6 ;
         params.la_s =   1e-1 ;
         params.alpha = 4;
         tic
         [a ,b  , c]  =  com_secrecy_results_with_eves_density(params);
         toc

    case 'LeakageAnalytical'

        params.la_u = 600e-6;
%         params.channel = 'Rice';
%         params.rice.Km = 32;    
%         params.rice.Ke = 32;
        params.alpha = 4;
        params.la_e = [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1 0.2:0.1:1 2:1:1e1] .* params.la_u;
        
        la_s = [1e-3 1e-1];
        for i = 1:numel(la_s)
            params.la_s = la_s(i);
            [Rm_math_4(i,:) , Re_math_4(i,:)  , Rs_math_4(i,:)]  =  com_secrecy_results_with_eves_density(params);
        end
        
        subplot(1,2,2);
        f1 = semilogx(params.la_e ./ params.la_u ,Re_math_4,'k-');
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1(1:2),{'$\lambda_s = 10^{-3} cells/m^2$' , '$\lambda_s = 10^{-1} cells/m^2$' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
        title('$\alpha = 4$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        
         params.alpha = 3;
         
         
         for i = 1:numel(la_s)
            params.la_s = la_s(i);
            [Rm_math_3(i,:) , Re_math_3(i,:)  , Rs_math_3(i,:)]  =  com_secrecy_results_with_eves_density(params);
         end
        
        subplot(1,2,1);
        f2 = semilogx(params.la_e ./ params.la_u ,Re_math_3,'k-');
        set(f2,'MarkerSize',15); 
        set(f2,'LineWidth',4);
        legend(f2(1:2),{'$\lambda_s = 10^{-3} cells/m^2$' , '$\lambda_s = 10^{-1} cells/m^2$' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
        title('$\alpha = 3$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
         
         
         
    case 'CellDensity'
        
        params.la_s =  [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1] ; 
        
        

%         params.channel = 'Rice';
%         params.rice.Km = 32;    
%         params.rice.Ke = 32;
        params.alpha = 4;
        la_e = [0.1 0.5 1] .* params.la_u;
        
        for i = 1:numel(la_e)
            params.la_e = la_e(i);
            [Rm_math_4(i,:) , Re_math_4(i,:)  , Rs_math_4(i,:)]  =  com_secrecy_results_with_smallcell_density(params);
            [Rm_simul_4(i,:), Re_simul_4(i,:), Rs_simul_4(i,:)]  =  gen_secrecy_results_with_smallcell_density_multiple_assoc(params);
            %[a, b, c]  =  gen_secrecy_results_with_smallcell_density_multiple_assoc(params);
        end    

            subplot(1,2,2);
            f1 = semilogx(params.la_s ,Re_simul_4,'ko' ,params.la_s,Re_math_4,'k-');
            set(f1,'MarkerSize',15);
            set(f1,'LineWidth',4);
            legend(f1(4:6),{'$\lambda_e / \lambda_u = 0.1$' , '$\lambda_e / \lambda_u = 0.5$' , '$\lambda_e / \lambda_u = 1$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');

        params.alpha = 3;

        for i = 1:numel(la_e)
            params.la_e = la_e(i);
            [Rm_math_3(i,:), Re_math_3(i,:), Rs_math_3(i,:)]  = com_secrecy_results_with_smallcell_density(params);
            [Rm_simul_3(i,:), Re_simul_3(i,:), Rs_simul_3(i,:)]  =  gen_secrecy_results_with_smallcell_density_multiple_assoc(params);
        end    

            subplot(1,2,1);
            f2 = semilogx(params.la_s ,Re_simul_3,'ko' ,params.la_s,Re_math_3,'k-');
            set(f2,'MarkerSize',15);
            set(f2,'LineWidth',4);
            legend(f2(4:6),{'$\lambda_e / \lambda_u = 0.1$' , '$\lambda_e / \lambda_u = 0.5$' , '$\lambda_e / \lambda_u = 1$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 3$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            figure(2);
            subplot(1,2,2);
            f3 = semilogx(params.la_s ,Rs_simul_4,'ko' ,params.la_s,Rs_math_4,'k-');
            set(f3,'MarkerSize',15);
            set(f3,'LineWidth',4);
            legend(f3(4:6),{'$\lambda_e / \lambda_u = 0.1$' , '$\lambda_e / \lambda_u = 0.5$' , '$\lambda_e / \lambda_u = 1$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            
            subplot(1,2,1);
            f4 = semilogx(params.la_s ,Rs_simul_3,'ko' ,params.la_s,Rs_math_3,'k-');
            set(f4,'MarkerSize',15);
            set(f4,'LineWidth',4);
            legend(f4(4:6),{'$\lambda_e / \lambda_u = 0.1$' , '$\lambda_e / \lambda_u = 0.5$' , '$\lambda_e / \lambda_u = 1$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 3$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
           
            figure(3);
            subplot(1,2,2);
            f5 = semilogx(params.la_s ,Rm_simul_4,'ko' ,params.la_s,Rm_math_4,'k-');
            set(f5,'MarkerSize',15);
            set(f5,'LineWidth',4);
            legend(f5(4:6),{'$\lambda_e / \lambda_u = 0.1$' , '$\lambda_e / \lambda_u = 0.5$' , '$\lambda_e / \lambda_u = 1$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average main link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            
            subplot(1,2,1);
            f6 = semilogx(params.la_s ,Rm_simul_3,'ko' ,params.la_s,Rm_math_3,'k-');
            set(f6,'MarkerSize',15);
            set(f6,'LineWidth',4);
            legend(f6(4:6),{'$\lambda_e / \lambda_u = 0.1$' , '$\lambda_e / \lambda_u = 0.5$' , '$\lambda_e / \lambda_u = 1$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
            ylabel('Average main link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 3$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
      
    case 'EvesDensity'
            la_s =  [1e-3 5e-3 1e-2];
            params.la_e =   (0.1:0.1:1) * params.la_u ;              
%             params.channel = 'Rice';
%             params.rice.K = 32;
             params.alpha = 4;
            
            for i = 1:numel(la_s)
                params.la_s = la_s(i);
                [Rm_math_4(i,:) , Re_math_4(i,:)  , Rs_math_4(i,:)]  =  com_secrecy_results_with_eves_density(params);
                [Rm_simul_4(i,:), Re_simul_4(i,:), Rs_simul_4(i,:)]  =  gen_secrecy_results_with_eves_density_multiple_assoc(params);
                %[a, b, c]  =  gen_secrecy_results_with_smallcell_density_multiple_assoc(params);
            end    

            subplot(1,2,1);
            f1 = semilogx(params.la_e ./ params.la_u ,Re_simul_4,'ko' ,params.la_e ./ params.la_u,Re_math_4,'k-');
            set(f1,'MarkerSize',15);
            set(f1,'LineWidth',4);
            legend(f1(4:6),{'$\lambda_s = 10^{-3} cells/m^2$' , '$\lambda_s = 5\times10^{-3} cells/m^2$' , '$\lambda_s = \times10^{-2} cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('$Eves relative density \lambda_e / \lambda_s (Eves/User)$  ','Interpreter','LaTex');
            ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');



        params.alpha = 5;

        for i = 1:numel(la_s)
            params.la_s = la_s(i);
            [Rm_math_5(i,:), Re_math_5(i,:), Rs_math_5(i,:)]  = com_secrecy_results_with_eves_density(params);
            [Rm_simul_5(i,:), Re_simul_5(i,:), Rs_simul_5(i,:)]  =  gen_secrecy_results_with_eves_density_multiple_assoc(params);
        end   

            subplot(1,2,2);
            f2 = semilogx(params.la_e ./ params.la_u ,Re_simul_5,'ko' ,params.la_e ./ params.la_u,Re_math_5,'k-');
            set(f2,'MarkerSize',15);
            set(f2,'LineWidth',4);
            legend(f2(4:6),{'$\lambda_s = 10^{-3} cells/m^2$' , '$\lambda_s = 5\times10^{-3} cells/m^2$' , '$\lambda_s = \times10^{-2} cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('$Eves relative density \lambda_e / \lambda_s (Eves/User)$  ','Interpreter','LaTex');
            ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 5$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
              
            figure(2);
            
            subplot(1,2,1);
            f3 = semilogx(params.la_e ./ params.la_u ,Rm_simul_4,'ko' ,params.la_e ./ params.la_u,Rm_math_4,'k-');
            set(f3,'MarkerSize',15);
            set(f3,'LineWidth',4);
            legend(f3(4:6),{'$\lambda_s = 10^{-3} cells/m^2$' , '$\lambda_s = 5\times10^{-3} cells/m^2$' , '$\lambda_s = \times10^{-2} cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('$Eves relative density \lambda_e / \lambda_s (Eves/User)$  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            subplot(1,2,2);
            f4 = semilogx(params.la_e ./ params.la_u ,Rm_simul_5,'ko' ,params.la_e ./ params.la_u,Rm_math_5,'k-');
            set(f4,'MarkerSize',15);
            set(f4,'LineWidth',4);
            legend(f4(4:6),{'$\lambda_s = 10^{-3} cells/m^2$' , '$\lambda_s = 5\times10^{-3} cells/m^2$' , '$\lambda_s = \times10^{-2} cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('$Eves relative density \lambda_e / \lambda_s (Eves/User)$  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
    
    case 'UserDensity'
        params.la_u = (100:100:1000) .* 1e-6; 
        params.la_e = 0.5 .* params.la_u;
%         params.channel = 'Rice';
%         params.rice.Km = 32;    
%         params.rice.Ke = 0;
        params.alpha = 4;
        la_s = [1e-3 1e-2 1e-1] ;
        
        for i = 1:numel(la_s)
            params.la_s = la_s(i);
            [Rm_math_4(i,:) , Re_math_4(i,:)  , Rs_math_4(i,:)]  =  com_secrecy_results_with_users_density(params);
            [Rm_simul_4(i,:), Re_simul_4(i,:), Rs_simul_4(i,:)]  =  gen_secrecy_results_with_users_density_multiple_assoc(params);
            %[a, b, c]  =  gen_secrecy_results_with_smallcell_density_multiple_assoc(params);
        end    

            subplot(1,2,2);
            f1 = plot(params.la_u * 1e6 ,Re_simul_4,'ko' ,params.la_u* 1e6,Re_math_4,'k-');
            set(f1,'MarkerSize',15);
            set(f1,'LineWidth',4);
            legend(f1(4:6),{'$\lambda_s = 10^{-3}~cells/m^2$' , '$\lambda_s = 10^{-2}~cells/m^2$' , '$\lambda_s = 10^{-1}~cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Users density (users/$km^2$)  ','Interpreter','LaTex');
            ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');

        params.alpha = 3;

        for i = 1:numel(la_s)
            params.la_s = la_s(i);
            [Rm_math_3(i,:), Re_math_3(i,:), Rs_math_3(i,:)]  = com_secrecy_results_with_users_density(params);
            [Rm_simul_3(i,:), Re_simul_3(i,:), Rs_simul_3(i,:)]  =  gen_secrecy_results_with_users_density_multiple_assoc(params);
        end    

            subplot(1,2,1);
            f2 = plot(params.la_u * 1e6 ,Re_simul_3,'ko' ,params.la_u * 1e6,Re_math_3,'k-');
            set(f2,'MarkerSize',15);
            set(f2,'LineWidth',4);
            legend(f2(4:6),{'$\lambda_s = 10^{-3}~cells/m^2$' , '$\lambda_s = 10^{-2}~cells/m^2$' , '$\lambda_s = 10^{-1}~cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Users density (users/$km^2$)  ','Interpreter','LaTex');
            ylabel('Average leakage link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 3$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            figure(2);
            subplot(1,2,2);
            f3 = plot(params.la_u * 1e6 ,Rs_simul_4,'ko' ,params.la_u * 1e6,Rs_math_4,'k-');
            set(f3,'MarkerSize',15);
            set(f3,'LineWidth',4);
            legend(f3(4:6),{'$\lambda_s = 10^{-3}~cells/m^2$' , '$\lambda_s = 10^{-2}~cells/m^2$' , '$\lambda_s = 10^{-1}~cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Users density (users/$km^2$)  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            
            subplot(1,2,1);
            f4 = plot(params.la_u * 1e6 ,Rs_simul_3,'ko' ,params.la_u * 1e6,Rs_math_3,'k-');
            set(f4,'MarkerSize',15);
            set(f4,'LineWidth',4);
            legend(f4(4:6),{'$\lambda_s = 10^{-3}~cells/m^2$' , '$\lambda_s = 10^{-2}~cells/m^2$' , '$\lambda_s = 10^{-1}~cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Users density (users/$km^2$)  ','Interpreter','LaTex');
            ylabel('Average secrecy rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 3$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
           
            figure(3);
            subplot(1,2,2);
            f5 = plot(params.la_u * 1e6 ,Rm_simul_4,'ko' ,params.la_u * 1e6,Rm_math_4,'k-');
            set(f5,'MarkerSize',15);
            set(f5,'LineWidth',4);
            legend(f5(4:6),{'$\lambda_s = 10^{-3}~cells/m^2$' , '$\lambda_s = 10^{-2}~cells/m^2$' , '$\lambda_s = 10^{-1}~cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Users density (users/$km^2$)  ','Interpreter','LaTex');
            ylabel('Average main link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 4$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
            
            
            subplot(1,2,1);
            f6 = plot(params.la_u * 1e6 ,Rm_simul_3,'ko' ,params.la_u * 1e6 ,Rm_math_3,'k-');
            set(f6,'MarkerSize',15);
            set(f6,'LineWidth',4);
            legend(f6(4:6),{'$\lambda_s = 10^{-3}~cells/m^2$' , '$\lambda_s = 10^{-2}~cells/m^2$' , '$\lambda_s = 10^{-1}~cells/m^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
            xlabel('Users density (users/$km^2$)  ','Interpreter','LaTex');
            ylabel('Average main link rate (nats/s/Hz)','Interpreter','LaTex'); 
            title('$\alpha = 3$','Interpreter','LaTex')
            grid on;
            set(gca, 'FontSize', 30);
            set(gca, 'FontWeight', 'Bold');
    case 'Interference'
        params.alpha = 4;
        %params.la_s =  [1e-4 3e-4 5e-4 8e-4 1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_s =  [1e-4 1e-3 1e-2 1e-1] ;
        params.la_e = 1 .* params.la_u;
         [Im,Ie] = gen_secrecy_results_with_smallcell_density_Interference(params);
         f2 = plot(params.la_s ,Im,'ko' ,params.la_s ,Ie,'k-');
%checkAssociation(params)
        
end


%     subplot(1,2,1);
%     f1 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-',params.la_s,R_math_ref,'r-');
%     set(f1,'MarkerSize',15);
%     set(f1,'LineWidth',4);
%     legend(f1([1 params.M+1]),{'Simulation' , 'Theorem (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
%     xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
%     ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex'); 
%     title('$$\alpha = 4$$','Interpreter','LaTex')
%     grid on;
%     set(gca, 'FontSize', 30);
%     set(gca, 'FontWeight', 'Bold');
%     
%     params.alpha = 5; 
%     R_math  = com_downlink_rate_with_smallcell_density(params);
%     R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
%     R_simul = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params);
%         
%     subplot(1,2,2);
%     f2 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-',params.la_s,R_math_ref,'r-');
%     title('$$\alpha = 5$$','Interpreter','LaTex')
%     set(f2,'MarkerSize',15);
%     set(f2,'LineWidth',4);
%     xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
%     ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex'); 
%     set(gca, 'FontSize', 30);
%     set(gca, 'FontWeight', 'Bold');
%     legend(f2([1 params.M+1]),{'Simulation' , 'Theorem (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
%     grid on;

end

% function  genSecrecyResults()
% clc;close all; clear;
% %format long e;
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %      Simulation Parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.channel = 'Rayleigh';
% params.la_s =  1e-2 ; 
% params.la_u =   (100:100:1000) * 1e-6 ;               % users density (users/m2)
% params.alpha = 4;                      % pass loss exponent
% params.M = 1;                                 
% params.P = 1; 
% 
% params.simulation_area_side = [-100 100];    % simulation area side to side
% params.space_realizations = 100;
% params.time_slots = 10;
% 
% 
% R_math  = com_downlink_rate_with_user_density(params);
% R_simul = gen_downlink_rate_with_user_density_ray_multiple_assoc(params);
% 
% 
%  
%     subplot(1,2,1);
%     f1 = plot(params.la_u*1e6 ,R_simul,'ko' ,params.la_u*1e6 ,R_math,'k-');
%     set(f1,'MarkerSize',15);
%     set(f1,'LineWidth',4);
%     legend(f1([1 params.M+1]),{'Simulation' , 'Theorem (2)'},'FontSize',25,'FontWeight','bold','Location','northeast');
%     xlabel('User density (users/$Km^2$)  ','Interpreter','LaTex');
%     ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex'); 
%     title('$$\alpha = 4$$','Interpreter','LaTex')
%     grid on;
%     set(gca, 'FontSize', 30);
%     set(gca, 'FontWeight', 'Bold');
%     xlim([100 1000]);
%     
%     params.alpha = 5; 
%     R_math  = com_downlink_rate_with_user_density(params);
%     R_simul = gen_downlink_rate_with_user_density_ray_multiple_assoc(params);
%         
%     subplot(1,2,2);
%     f2 = plot(params.la_u*1e6 ,R_simul,'ko' ,params.la_u*1e6,R_math,'k-');
%     title('$$\alpha = 5$$','Interpreter','LaTex')
%     set(f2,'MarkerSize',15);
%     set(f2,'LineWidth',4);
%     xlabel('User density (cells/$Km^2$)  ','Interpreter','LaTex');
%     ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex'); 
%     set(gca, 'FontSize', 30);
%     set(gca, 'FontWeight', 'Bold');
%     legend(f2([1 params.M+1]),{'Simulation' , 'Theorem (2)'},'FontSize',25,'FontWeight','bold','Location','northeast');
%     xlim([100 1000]);
%     grid on;
%     
%     disp(sprintf ( '\n') )
% end

%%
function checkAssociation(params)
        
        simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;

        mu_s = params.la_s(10) * simulation_area;
        mu_u = params.la_u * simulation_area;
        mu_e = params.la_e * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        N_eaves = poissrnd(mu_e);
        cmap = hsv(N_users);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        eaves_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_eaves, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        r_se = pdist2(cells_pos,eaves_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        MDE = zeros(N_users,1);  % association matrix

        % Association  [Long_Term Association i.e. the channel is averaged to 1] 
         [r_u , Server_u] = min(r_su);
         
        for i = 1:N_users;
            [~ , MDE(i)] = min(r_se(Server_u(i),:));
            
        end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %      Network Plot
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 s = plot(cells_pos(:,1),cells_pos(:,2),'sb');
 for i = 1:N_cells
    text(cells_pos(i,1) - 2,cells_pos(i,2),num2str(i),'Color','blue','FontSize',14);
 end
 set(s,'MarkerSize',20);
 xlabel('x (meters)');
 ylabel('y (meters)');
 set(gca, 'FontName', 'Arial');
 set(gca, 'FontSize', 20);
 set(gca, 'FontWeight', 'Bold');
 
 hold on;
 
 u = plot(users_pos(:,1),users_pos(:,2),'or');
 set(u,'MarkerSize',15);
 for i = 1:N_users
    text(users_pos(i,1) - 2,users_pos(i,2),num2str(i),'Color','red','FontSize',14);
 end

 
 e = plot(eaves_pos(:,1),eaves_pos(:,2),'pk');
 set(e,'MarkerSize',30);
 for i = 1:N_eaves
    text(eaves_pos(i,1) - 2,eaves_pos(i,2),num2str(i),'Color','black','FontSize',14);
 end
 
         for i = 1:N_users;
             x = [users_pos(i,1) , cells_pos(Server_u(i),1) ];
             y = [users_pos(i,2) , cells_pos(Server_u(i),2) ];
             xe = [eaves_pos(MDE(i),1) , cells_pos(Server_u(i),1) ];
             ye = [eaves_pos(MDE(i),2) , cells_pos(Server_u(i),2) ];
            line(x ,y,'Color','red','LineStyle','-','LineWidth',3)
            line(xe ,ye,'Color',cmap(i,:),'LineStyle',':','LineWidth',2)
            
        end
end

function [Rm,Re,Rs] = com_secrecy_results_with_eves_density(params)

points = numel(params.la_e);

R_math = zeros(points,params.M);

% Rice Channel Parameters
    
for p = 1:points
    k = params.la_s / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
    if strcmp(params.channel ,'Rayleigh')
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
    elseif strcmp(params.channel ,'Rice')
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  (1 - (2/(mo*a+2))*(mo ./(omega.*z)).^mo  .*  hyp2f1(mo,mo + 2/a,mo + 1 + 2/a,-mo./(omega.*z))) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
    end
    
     kes = params.la_s / params.la_e(p);
     Fnl = @(z,a,mo,omega,mi,po) (1 - (1 + z.*omega./mo).^-mo)./(z .*(1 - (1-po).*(kes) + (1-po).*(kes) .* ((1 + z./mi).^-mi + ...
                                      mi.*(mi).^mi.*(1 - 2/a).^-1 .* z .* (z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,z.*(z+mi).^-1)))) ;
   Re(p) = integral(@(z)Fnl(z,params.alpha,m_e,omega_e,1,po),0,inf);
   Rm(p) = R_math(p)
   Rs(p) = Rm(p) - Re(p);
end

end

function [Rm,Re,Rs] = com_secrecy_results_with_smallcell_density(params)

points = numel(params.la_s);

Rm = zeros(points,params.M);
Rs = zeros(points,params.M);
Re = zeros(points,1);


    
for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
      if strcmp(params.channel ,'Rayleigh')
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
      elseif strcmp(params.channel ,'Rice')
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  (1 - (2/(mo*a+2))*(mo ./(omega.*z)).^mo  .*  hyp2f1(mo,mo + 2/a,mo + 1 + 2/a,-mo./(omega.*z))) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
    end
    
    kes = params.la_s(p) / params.la_e;
    Fnl = @(z,a,mo,omega,mi,po) (1 - (1 + z.*omega./mo).^-mo)./(z .*(1 - (1-po).*(kes) + (1-po).*(kes) .* ((1 + z./mi).^-mi + ...
                                      mi.*(mi).^mi.*(1 - 2/a).^-1 .* z .* (z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,z.*(z+mi).^-1)))) ;
    Re(p) = integral(@(z)Fnl(z,params.alpha,1,1,1,po),0,inf);
    Rm(p) = R_math(p);
    Rs(p) = Rm(p) - Re(p);
    if(Rs(p) < 0)
       Rs(p)= 0;
    end
end

end

function [Rm,Re,Rs] = com_secrecy_results_with_smallcell_density_noise(params)

points = numel(params.la_s);

Rm = zeros(points,params.M);
Rs = zeros(points,params.M);
Re = zeros(points,1);


    
for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
      if strcmp(params.channel ,'Rayleigh')
        switch(params.M)
            case 1
                
                snr = params.Ps / params.No ; 
                F = @(y,a,mo,omega,mi,po,snr,n) (1 - (1 + y.*snr.*omega./mo).^-mo)./(y .*(po + (1-po) .* ((1 + params.M.*y.*snr./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* y.*snr .* (params.M.*y.*snr+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*y.*snr.*(params.M.*y.*snr+mi).^-1))).^n) ;
                Iy(p) = integral(@(y)F(y,params.alpha,1,1,1,po,snr,params.M),0,inf); 
                
                G = @(y,r,a,mo,omega,mi,po,snr,n) a * r.^(a-1) .* exp(- y .* r.^a) .* exp(-pi.*params.la_s(p).*r.^2.* (po + (1-po) .* ((1 + params.M.*y.*snr./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* y.*snr .* (params.M.*y.*snr+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*y.*snr.*(params.M.*y.*snr+mi).^-1))).^n) .* (1 - (1 + y.*snr.*omega./mo).^-mo)./(y.*snr .*(po + (1-po) .* ((1 + params.M.*y.*snr./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* y.*snr .* (params.M.*y.*snr+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*y.*snr.*(params.M.*y.*snr+mi).^-1))).^n) ;
                                   
                Iyr(p) = integral2(@(y,r) G(y,r,params.alpha,1,1,1,po,snr,params.M),0,inf,0,inf);
                R_math(p) = Iy(p) - Iyr(p);
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
      elseif strcmp(params.channel ,'Rice')
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  (1 - (2/(mo*a+2))*(mo ./(omega.*z)).^mo  .*  hyp2f1(mo,mo + 2/a,mo + 1 + 2/a,-mo./(omega.*z))) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
    end
    
    kes = params.la_s(p) / params.la_e;
    Fnl = @(y,a,mo,omega,mi,po) (1 - (1 + y.*snr.*omega./mo).^-mo)./(y .*(1 - (1-po).*(kes) + (1-po).*(kes) .* ((1 + y.*snr./mi).^-mi + ...
                                      mi.*(mi).^mi.*(1 - 2/a).^-1 .* y.*snr .* (y.*snr+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,y.*snr.*(y.*snr+mi).^-1)))) ;
    Iy(p) = integral(@(y)Fnl(y,params.alpha,1,1,1,po),0,inf);
    
    G = @(y,r,a,mo,omega,mi,po,snr) a * r.^(a-1) .* exp(- y .* r.^a) .* exp(-pi.*params.la_s(p).*r.^2.* (1 - (1-po).*(kes) + (1-po).*(kes) .* ((1 + params.M.*y.*snr./mi).^-mi + ...
        mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* y.*snr .* (params.M.*y.*snr+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*y.*snr.*(params.M.*y.*snr+mi).^-1)))) .* (1 - (1 + y.*snr.*omega./mo).^-mo)./(y.*snr .*(1 - (1-po).*(kes) + (1-po).*(kes) .* ((1 + params.M.*y.*snr./mi).^-mi + ...
        mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* y.*snr .* (params.M.*y.*snr+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*y.*snr.*(params.M.*y.*snr+mi).^-1)))) ;
    
    Iyr(p) = integral2(@(y,r) G(y,r,params.alpha,1,1,1,po,snr),0,inf,0,inf);
    
    Re(p) = Iy(p) - Iyr(p)
    Rm(p) = R_math(p);
    Rs(p) = Rm(p) - Re(p);
    if(Rs(p) < 0)
       Rs(p)= 0;
    end
end

end

function [Rm,Re,Rs] = com_secrecy_results_with_users_density(params)

points = numel(params.la_u);

Rm = zeros(points,params.M);
Rs = zeros(points,params.M);
Re = zeros(points,1);


    
for p = 1:points
    k = params.la_s / params.la_u(p);   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
      if strcmp(params.channel ,'Rayleigh')
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
      elseif strcmp(params.channel ,'Rice')
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
        switch(params.M)
            case 1
                F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p) = integral(@(z)F(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 2
                F21 = @(z,a,mo,omega,mi,po,n)  (1 - (2/(mo*a+2))*(mo ./(omega.*z)).^mo  .*  hyp2f1(mo,mo + 2/a,mo + 1 + 2/a,-mo./(omega.*z))) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F21(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F22(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            case 3
                F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F31(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F32(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F33(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

            case 4
                F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;

                R_math(p,1) = integral(@(z)F41(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,2) = integral(@(z)F42(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,3) = integral(@(z)F43(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 

                F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                              mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
                R_math(p,4) = integral(@(z)F44(z,params.alpha,m_m,omega_m,1,po,params.M),0,inf); 
            otherwise
                disp('Invalid Multicell size (M)')
        end
    end
    
    kes = params.la_s / params.la_e(p);
    Fnl = @(z,a,mo,omega,mi,po) (1 - (1 + z.*omega./mo).^-mo)./(z .*(1 - (1-po).*(kes) + (1-po).*(kes) .* ((1 + z./mi).^-mi + ...
                                      mi.*(mi).^mi.*(1 - 2/a).^-1 .* z .* (z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,z.*(z+mi).^-1)))) ;
    Re(p) = integral(@(z)Fnl(z,params.alpha,m_e,omega_e,1,po),0,inf);
    Rm(p) = R_math(p);
    Rs(p) = Rm(p) - Re(p);
    if(Rs(p) < 0)
       Rs(p)= 0;
    end
end

end

function [R_math_ref] = com_downlink_rate_with_smallcell_densityRef(params)

points = numel(params.la_s);

R_math_ref = zeros(points,1);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K; 
    mo = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end
    
for p = 1:points
   
   
    F = @(z,a,mo,omega,mi,n) (1 - (1 + z.*omega./mo).^-mo)./(z .* ((1 + z./mi).^-mi + ...
                                  mi.*(mi).^mi.*(1 - 2/a).^-1 .*  z .* (z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1)).^n) ;
    R_math_ref(p) = integral(@(z)F(z,params.alpha,1,1,1,1),0,inf); 
        
    
end

end

function [Rm,Re,Rs] = gen_secrecy_results_with_eves_density_multiple_assoc(params)
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_e);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
Rate_S = zeros(points ,params.M, params.space_realizations , params.time_slots );
Rate_E = zeros(points ,params.M, params.space_realizations , params.time_slots );
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Eves Density: ' , num2str(params.la_e(p))]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u * simulation_area;
        mu_e = params.la_e(p) * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        N_eaves = poissrnd(mu_e);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        eaves_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_eaves, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        r_se = pdist2(cells_pos,eaves_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        MDE = zeros(N_users,1);  
        S_u = zeros(N_users,params.M);
        S_e = zeros(N_users,1);  
        I_e = zeros(N_users,1); 
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1] 
        %[r_u_ , Server_u_] = min(r_su);
         for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
            [~ , MDE(i)] = min(r_se(Server_u(i),:));
         end
         
         
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            He = exprnd(1,N_cells,N_eaves);
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            R_e =  params.P *  He.* r_se.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                    S_e(i,s) = R_e(Server_u(i,s),MDE(i));
                end

            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = Int_all - Int_servers;
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                    
                end
                Hum = gamrnd(mo,omega/mo) ;
                S_u(i) = params.P *  Hum.* r_u(i).^-params.alpha;
                I_e(i) = sum(Ac.*R_e(:,MDE(i))) - S_e(i);
            end
            Rate_u = log(1 + SIR_u);
            SIR_e = S_e./I_e;
            Rate_e = log(1 + SIR_e);
            Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_s = Rate_u - Rate_e;
            Rate_s = Rate_s(Rate_s > 0);
            Rate_S(p,:,m,t) = sum(Rate_s) / N_users;
            Rate_E(p,:,m,t) = sum(Rate_e) / N_users;
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Rate_U(isnan(Rate_U)) = 0;
Rate_E(isnan(Rate_E)) = 0;
Rate_S(isnan(Rate_S)) = 0;
Rate_U(isinf(Rate_U)) = 0;
Rate_E(isinf(Rate_E)) = 0;
Rate_S(isinf(Rate_S)) = 0;

Rate_U(Rate_U == 0) = max(Rate_U(:));
Rate_E(Rate_E == 0) = max(Rate_E(:));
Rate_S(Rate_S == 0) = max(Rate_S(:));

Rm = sum(sum(Rate_U,4),3) / normfact;
Re = sum(sum(Rate_E,4),3) / normfact;
Rs = sum(sum(Rate_S,4),3) / normfact;
end

function [Rm,Re,Rs] = gen_secrecy_results_with_smallcell_density_multiple_assoc(params)
H = exprnd(1,1e5,1e3) ;
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
    
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
Rate_S = zeros(points ,params.M, params.space_realizations , params.time_slots );
Rate_E = zeros(points ,params.M, params.space_realizations , params.time_slots );
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp(['Eves Density: ' , num2str(params.la_e)]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        mu_e = params.la_e * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        N_eaves = poissrnd(mu_e);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        eaves_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_eaves, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        r_se = pdist2(cells_pos,eaves_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        MDE = zeros(N_users,1);  
        S_u = zeros(N_users,params.M);
        S_e = zeros(N_users,1);  
        I_e = zeros(N_users,1); 
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1] 
        %[r_u_ , Server_u_] = min(r_su);
         for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
            [~ , MDE(i)] = min(r_se(Server_u(i),:));
         end
         
         
        for t = 1:params.time_slots
            toss = round(rand(1)*100)+1;
            %Hui = exprnd(1,N_cells,N_users) ;
            Hui = H(t:N_cells+t-1,t:N_users+t-1) ;
            He = H(toss:N_cells+toss-1,toss:N_eaves+toss-1) ;
            %He = exprnd(1,N_cells,N_eaves);
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            R_e =  params.P *  He.* r_se.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                    S_e(i,s) = R_e(Server_u(i,s),MDE(i));
                end

            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = Int_all - Int_servers;
                    dominter = getInterfernceOfDominantInterferers(Server_u,MDE,N_eaves,R_e);
            for i = 1:N_users;
                I_e(i) = sum(Ac.*R_e(:,MDE(i))) - dominter (MDE(i)) ;%- S_e(i);
                Hum = gamrnd(m_m,omega_m/m_m) ;
                Hee = gamrnd(m_e,omega_e/m_e) ;
                S_u(i) = params.P *  Hum.* r_u(i).^-params.alpha;
                S_e(i) = params.P *  Hee.* r_se(Server_u(i,1),MDE(i)).^-params.alpha;
                for s = 1:params.M
                    SIR_u(i) = S_u(i) / I_u(i);
                    
                end 
            end
            Rate_u = log(1 + SIR_u);
            SIR_e = S_e./I_e;
            Rate_e = log(1 + SIR_e);
            Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_s = Rate_u - Rate_e;
            Rate_s = Rate_s(Rate_s > 0);
            Rate_S(p,:,m,t) = sum(Rate_s) / N_users;
            Rate_E(p,:,m,t) = sum(Rate_e) / N_users;
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Rate_U(isnan(Rate_U)) = 0;
Rate_E(isnan(Rate_E)) = 0;
Rate_S(isnan(Rate_S)) = 0;
Rate_U(isinf(Rate_U)) = 0;
Rate_E(isinf(Rate_E)) = 0;
Rate_S(isinf(Rate_S)) = 0;

Rate_U(Rate_U == 0) = max(Rate_U(:));
Rate_E(Rate_E == 0) = max(Rate_E(:));
Rate_S(Rate_S == 0) = max(Rate_S(:));

Rm = sum(sum(Rate_U,4),3) / normfact;
Re = sum(sum(Rate_E,4),3) / normfact;
Rs = sum(sum(Rate_S,4),3) / normfact;
end

function [Rm,Re,Rs] = gen_secrecy_results_with_users_density_multiple_assoc(params)
H = exprnd(1,1e5,1e3) ;
    Km = params.rice.Km; 
    m_m = (Km+1)^2/(2*Km+1);
    sigma_m = 1;
    omega_m = (2*Km+1) * sigma_m^2;
    
    Ke = params.rice.Ke; 
    m_e = (Ke+1)^2/(2*Ke+1);
    sigma_e = 1;
    omega_e = (2*Ke+1) * sigma_e^2;  
    
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_u);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
Rate_S = zeros(points ,params.M, params.space_realizations , params.time_slots );
Rate_E = zeros(points ,params.M, params.space_realizations , params.time_slots );
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Smallcell Density: ' , num2str(params.la_s)]);
    disp(['Eves Density: ' , num2str(params.la_e(p))]);
    disp(['Users Density: ' , num2str(params.la_u(p))]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u(p) * simulation_area;
        mu_e = params.la_e(p) * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        N_eaves = poissrnd(mu_e);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        eaves_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_eaves, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        r_se = pdist2(cells_pos,eaves_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        MDE = zeros(N_users,1);  
        S_u = zeros(N_users,params.M);
        S_e = zeros(N_users,1);  
        I_e = zeros(N_users,1); 
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1] 
        %[r_u_ , Server_u_] = min(r_su);
         for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
            [~ , MDE(i)] = min(r_se(Server_u(i),:));
         end
         
         
        for t = 1:params.time_slots
            toss = round(rand(1)*100)+1;
            %Hui = exprnd(1,N_cells,N_users) ;
            Hui = H(t:N_cells+t-1,t:N_users+t-1) ;
            He = H(toss:N_cells+toss-1,toss:N_eaves+toss-1) ;
            %He = exprnd(1,N_cells,N_eaves);
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            R_e =  params.P *  He.* r_se.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                    S_e(i,s) = R_e(Server_u(i,s),MDE(i));
                end

            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = Int_all - Int_servers;
                    dominter = getInterfernceOfDominantInterferers(Server_u,MDE,N_eaves,R_e);
            for i = 1:N_users;
                I_e(i) = sum(Ac.*R_e(:,MDE(i))) - dominter (MDE(i)) ;%- S_e(i);
                Hum = gamrnd(m_m,omega_m/m_m) ;
                Hee = gamrnd(m_e,omega_e/m_e) ;
                S_u(i) = params.P *  Hum.* r_u(i).^-params.alpha;
                S_e(i) = params.P *  Hee.* r_se(Server_u(i,1),MDE(i)).^-params.alpha;
                for s = 1:params.M
                    SIR_u(i) = S_u(i) / I_u(i);
                    
                end 
            end
            Rate_u = log(1 + SIR_u);
            SIR_e = S_e./I_e;
            Rate_e = log(1 + SIR_e);
            Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_s = Rate_u - Rate_e;
            Rate_s = Rate_s(Rate_s > 0);
            Rate_S(p,:,m,t) = sum(Rate_s) / N_users;
            Rate_E(p,:,m,t) = sum(Rate_e) / N_users;
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Rate_U(isnan(Rate_U)) = 0;
Rate_E(isnan(Rate_E)) = 0;
Rate_S(isnan(Rate_S)) = 0;
Rate_U(isinf(Rate_U)) = 0;
Rate_E(isinf(Rate_E)) = 0;
Rate_S(isinf(Rate_S)) = 0;

Rate_U(Rate_U == 0) = max(Rate_U(:));
Rate_E(Rate_E == 0) = max(Rate_E(:));
Rate_S(Rate_S == 0) = max(Rate_S(:));

Rm = sum(sum(Rate_U,4),3) / normfact;
Re = sum(sum(Rate_E,4),3) / normfact;
Rs = sum(sum(Rate_S,4),3) / normfact;
end

function [Im,Ie] = gen_secrecy_results_with_smallcell_density_Interference(params)
 
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Interference_U = zeros(points, params.space_realizations , params.time_slots);
Interference_E = zeros(points, params.space_realizations , params.time_slots );
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp(['Eves Density: ' , num2str(params.la_e)]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        mu_e = params.la_e * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        N_eaves = poissrnd(mu_e);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        eaves_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_eaves, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        r_se = pdist2(cells_pos,eaves_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        MDE = zeros(N_users,1);  
        S_u = zeros(N_users,params.M);
        S_e = zeros(N_users,1);  
        I_e = zeros(N_users,1);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1] 
        %[r_u_ , Server_u_] = min(r_su);
         for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
            [~ , MDE(i)] = min(r_se(Server_u(i),:));
         end
         
         
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            He = exprnd(1,N_cells,N_eaves);
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            R_e =  params.P *  He.* r_se.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                    S_e(i,s) = R_e(Server_u(i,s),MDE(i));
                end

            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = Int_all - Int_servers;
            
    
            for i = 1:N_users;
                I_e(i) = sum(Ac.*R_e(:,MDE(i))) - S_e(i);
            end
            Interference_u = I_u;
            Interference_e = I_e;
            
            %Interference_e(Interference_e >= 100 * min(Interference_e)) = min(Interference_e);
            Interference_U(p,m,t) = sum(Interference_u) / N_users;
            Interference_E(p,m,t) = sum(Interference_e) / N_users;
            
            
        end
    end
    
end

normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Im = sum(sum(Interference_U,3),2) / normfact;
Ie = sum(sum(Interference_E,3),2) / normfact;


end

function [servers,r_servers] = getServers(user,distances, N)
    servers = zeros(N,1);
    r_servers = zeros(N,1);
    for i = 1:N
        [r_servers(i) , servers(i)] = min(distances(:,user));
        distances(servers(i),user) = inf;
    end
end


function inter = getInterfernceOfDominantInterferers(servers,MDE,No_Eves,R_signal)
    inter = zeros(No_Eves,1);
    for mde = 1 : No_Eves
        dominterferes = unique(servers(MDE == mde));
        if(isempty(dominterferes))
            inter(mde)= 0;
            continue;
        end
        inter(mde) = sum(R_signal(dominterferes,mde));
    end
end

