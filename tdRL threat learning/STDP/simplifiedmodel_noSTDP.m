function V = simplifiedmodel_noSTDP

    V_rest = -60;
    V_th = -40;
    V_react = 20;
    V_reset = -70;
    tau = 10;
    g_inh = 0.015;
    g_ext = 0.075;
    
    totaltime = 1000;
    neuron_num = 100;
    
    interneuron_num = 40;
    groupA_num = 30;
    groupB_num = 30;
    
    W = zeros (neuron_num, neuron_num);
    for i = 1: neuron_num
        for j = 1: neuron_num
            if i<=interneuron_num && interneuron_num<j && j<=interneuron_num+groupA_num
                s = rand;
                if s>0.8
                    W(i,j)=-1;
                end
            elseif  i<=interneuron_num && j > interneuron_num+groupA_num
                s = rand;
                if s > 0.5
                    W(i,j) = -1;
                end
            elseif interneuron_num<i && i<=neuron_num && interneuron_num<j && j<=neuron_num
                s = rand;
                if s > 0.9
                    W(i,j) = 1;
                end
            end
        end
    end
    
    
    V = zeros(totaltime, neuron_num) + V_rest;
    N = poissrnd(10, totaltime, neuron_num)-10;
    firing = zeros(totaltime,neuron_num);
    window = zeros(totaltime, neuron_num);
    max_LTP = 0;
    max_LTD = 0;
    updated = zeros(totaltime, neuron_num);
    
    % 五次学习
    I = zeros(totaltime,neuron_num);
    I (50:59, interneuron_num+1:neuron_num) = 10; 
    I  (100: 109,  interneuron_num+1:neuron_num) = 10;
    I (200: 209,  interneuron_num+1:neuron_num) = 10;
    I (300:309,  interneuron_num+1:neuron_num) = 10;
    I (400: 409,  interneuron_num+1:neuron_num) = 10;
    
    I (50:53, 1: interneuron_num) = 10; 
    I  (100: 103,  1: interneuron_num) = 15;
    I (200: 203,  1: interneuron_num) = 20;
    I (300:303,  1: interneuron_num) = 25;
    I (400: 403,  1: interneuron_num) = 30;
    

    W_update = W;
    for i = 2:totaltime
        for j = 1: neuron_num
            dV = -(V(i-1, j)-V_rest)/tau + I(i-1,j) + N(i-1,j);
            interactionV = V(i-1,:).*(V(i-1, :)>V_th);
            dS  = sum(W_update(1:interneuron_num, j).*interactionV(1:interneuron_num)'.* g_inh) + sum(W_update(interneuron_num+1:neuron_num, j).*interactionV(interneuron_num+1:neuron_num)'.* g_ext);
            dV = dV + dS;
            V(i,j) = V(i-1, j)+dV;
    
             if firing(i-1,j)==1 && firing(i,j)==0
                 V(i,j) = V_reset;
             end
    
             if V(i-1, j) + dV >V_th && firing(i-1, j)==0
                 upper_limit = min([i+5, totaltime]);
                 firing(i:upper_limit, j) = 1;
                  V(i,j) = V_react;
             end        
        end
    
      for j = interneuron_num+1:neuron_num
          if firing(i-1,j)==0 && firing(i,j)==1
              lower_limit = max([1, i-5]);
              for m = 1: i - lower_limit
                  window(i-m, j) = 1+max_LTP*(1-1/5*(m-1));
              end
          end
          if V(i,j)==V_reset
              upper_limit = min([totaltime, i+5]);
              for m = 1:upper_limit-i
                 window(i+m, j) = 1-max_LTD*(1-1/5*(m-1));
              end
          end
      end
      
      for j = interneuron_num+1:neuron_num
          temp = window(1:i, j);
          temp2 = updated(1:i,j);
          update_time = find(temp>0 & temp2==0);
          input_cell_idx = find(W(:, j)>0);
          for m = 1:length(input_cell_idx)
              temp3 = firing(1:i, input_cell_idx(m));
              firing_time = find(temp3>0);
              if ~isempty(intersect(update_time, firing_time))
                  update_idx = intersect(update_time, firing_time);
                  if mean(temp(update_idx))>1
                      W_update(input_cell_idx(m),j) = W_update(input_cell_idx(m),j)*max(temp(update_idx));
                  end
                  if mean(temp(update_idx))<1
                      W_update(input_cell_idx(m),j) = W_update(input_cell_idx(m),j)*min(temp(update_idx));
                  end
                  updated(update_idx, j) = 1;
              end
          end
      end



%         for j = 1:neuron_num
%             if i > 2 &&  interneuron_num<j
%                 temp_act = V(i-1, j);
%                 if temp_act>1
%                     potentiation_input = find(max(V(i-2:i-1, :), [], 1)>V_th);
%                     potentiation_input(potentiation_input<=interneuron_num)=[];
%                     depression_input = find(V(i, :)>V_th);
%                     depression_input(depression_input<=interneuron_num)=[];
%                     W_update (potentiation_input, j) = W_update (potentiation_input, j).*1.1;
%                     W_update (depression_input, j) = W_update (depression_input, j).*0.9;
%                     potentiation_output = find(max(V(i-1:i, :), [], 1)>V_th);
%                     depression_output = find(V(i-2, :)>V_th);
%                     potentiation_output(potentiation_output<=interneuron_num)=[];
%                     depression_output(depression_output<=interneuron_num)=[];
%                     W_update (j, potentiation_output) = W_update (j, potentiation_output).*1.1;
%                     W_update (j, depression_output) = W_update (j, depression_output).*0.9;
%                 end
%             end
%         end
    end
    
    W_change = W_update-W;
    group1_up= sum(W_change(41:70, 41:70)>0, 'all')/900;
    group1_dn= sum(W_change(41:70, 41:70)<0, 'all')/900;
    group2_up= sum(W_change(71:100, 71:100)>0, "all")/900;
    group2_dn= sum(W_change(71:100, 71:100)<0, "all")/900;
    
    inter_up= sum([W_change(71:100, 41:70); W_change(41:70, 71:100)]>0, "all")/900;
    inter_dn= sum([W_change(71:100, 41:70); W_change(41:70, 71:100)]<0, "all")/900;
    
    inter = [W_change(71:100, 41:70); W_change(41:70, 71:100)];
    inter_change = reshape(inter, [],1);
    inter_nonzero = inter_change(~inter_change==0);
%  figure 
%  imagesc(V')
%  figure
%  plot(mean(V(:, 41:70), 2));
%  hold on
%  plot(mean(V(:, 71:100), 2))
%  hold off
% 
% figure
% imagesc(W_change)    
end



    
