% SIGMA_WRITE_DOT
% Ross M. Richardson
% Fan Chung Graham's Research Group, UCSD
% Begun 7/13/05

function sigma_write_dot(G,W,name, cutoff,grey)

% This file is part of a private research program to compute
% the sigma function and associated weight vector of a graph G
% (see Spectral Graph Theory, by Dr. Fan Chung Graham, pg. 104).
%
% Writes a DOT file to filename containing the graph information
% plus a weight for each edge spring.
if (nargin == 3)
  cutoff = 0;	
end
                
filename = [name, '.dot'];
fid = fopen(filename, 'w');

if (fid == -1)
  error('Error opening file');
end

fprintf(fid,'graph %s {\n', name);

[m,n] = size(G);

% Normalize W
W = W / max(vec(W));
for i = 1:n
  for j = i:n
    if (G(i,j) == 1 & i ~= j)
      if (W(i,j) <= cutoff)
        if (grey == 1)
          fprintf(fid,'%d -- %d [weight="%e", color="0.0 0.0 %e"];\n', i, ...
                  j, W(i,j), 0.75 - W(i,j)*0.75); % DO GREYSCALE FOR PUBLICATION 
        else
          fprintf(fid,'%d -- %d [weight="%e", color="%e 1.0 1.0"];\n', i, ...
                  j, W(i,j), W(i,j)/2 + .5);
        end
      end
    end
  end
end

% Normalize W
W = W / max(vec(W));
for i = 1:n
  for j = i:n
    if (G(i,j) == 1 & i ~= j)
      if (W(i,j) > cutoff)
        if (grey == 1)
          fprintf(fid,'%d -- %d [weight="%e", color="0.0 0.0 %e" style="bold"];\n', i, ...
                  j, W(i,j), 0.75 - W(i,j)*0.75); % DO GREYSCALE FOR PUBLICATION 
        else
          fprintf(fid,'%d -- %d [weight="%e", color="%e 1.0 1.0" style="bold"];\n', i, ...
                  j, W(i,j), W(i,j)/2 + .5);
        end
      end
    end
  end
end



sums = sum(W);
maxdeg =  max(sum(W));
for i = 1:n
  %fprintf(fid,'%d [color="%e 1.0 1.0"];\n', i, sum(G(i,:))/ ...
  %        (2*maxdeg) + .5);
  if (grey == 1)
    fprintf(fid,'%d [color="0.0 0.0 %e"];\n', i, 0.75 - ... 
	(sum(W(i,:))/ maxdeg) * 0.75); 
  else
    fprintf(fid,'%d [color="%e 1.0 1.0"];\n', i, sum(W(i,:))/ ...
            (2*maxdeg) + .5);
  end
end

fprintf(fid,'}\n');

fclose(fid);
