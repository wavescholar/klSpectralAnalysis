function [best_adj, best_score] = max_card_match(adj, best_score, Nedgs, lim_score)
% [best_adj, best_score] = max_card_match((adj > 0), 0, Nedgs, floor(Nnds/2))
% Finds Maximum Cardinality matching in un-directed non  bipartite graph
% represented by (symmetric) ADJacency matrix  "adj" with 1's
% OUTPUT: ADJacency matrix of the maximal matching and the number of
%         edges in that matching. 
%  Relies on symmetric adj mat = undirected graph, best_score = 0 initially
% This lame function works by recursive depth-first pruning and 
% has exponential complexity
% There are polinomial (linear) algorithms  

best_adj = adj;
if sum(sum(adj) > 1) > 0    % If there are nodes with >1 edge  ==> not finished pruning
    if (Nedgs - 1 > best_score)      % potentially better matching 
        [i, j] = find(triu(adj));          % Make arrays of indices of edges  
        bkp_adj = adj;             % Save a backup copy of ADJ 
        for ndx = 1:Nedgs              
          if (best_score >= Nedgs-1) 
              % fprintf('Best score %d, NOT calling on %d edged graph, passing up \n', best_score,Nedgs-1);
              return
          end
          adj = bkp_adj;                                   % restore the original ADJ 
          adj(i(ndx),j(ndx)) = 0; adj(j(ndx),i(ndx)) = 0;         % erase one edge
          % fprintf('Best score %d, calling on %d edged graph \n', best_score,Nedgs-1);
          [adj, score] = max_card_match(adj, best_score, Nedgs-1, lim_score);
          if (score > best_score)
              best_score = score; 
              best_adj = adj;
              if (best_score == lim_score) % found SOME matching of best possible score
                  break
              end
          end
        end
    else
        % fprintf('Best score %d, NOT calling on %d edged graph \n', best_score,Nedgs-1);
        best_score = 0;
        best_adj = 0;  
    end
else
    best_score = Nedgs; 
    fprintf('Found B-matching on %d edges\n', Nedgs);
end