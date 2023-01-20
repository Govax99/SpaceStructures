function [K,F] = fem(data_file)
%FEM apply fem method to a structure of beam and trusses
%
% PROTOTYPE:
%     [K,F] = fem(data_file)
%
% INPUT:
%    data_file[dim]       file containing all the data, follows example
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE OF DATA_FILE BEGINS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEMENT
% model: truss %you must put model or shape functions
% type: axial
% EA: 1e5
% length: 1000
% shape functions: [ (1-x/l) x/l] %not needed if you put model
% global index: [4 5]
% inclination: 30
% 
% ELEMENT
% model: beam %you must put model or shape functions
% type: bending
% EJ: 1e7
% length: 1000
% global index: [1 2 0 0]
% FORCE
% type: concentrated
% index: 4 %index or element -> depends on if distributed or concentrated
% value: 7
% 
% FORCE
% type: distributed
% element: 1
% shape: 3*(x/1000)
% 
% SPRING
% index: 1
% k: 3000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE OF DATA_FILE ENDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT:
%    K[dim]           global stiffness matrix [-]
%    F[dim]           global force vector [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    [elements,forces,springs] = read_data_file(data_file);

    [elements, forces] = build_local(elements, forces);

    [K,F] = build_global(elements,forces,springs);
end

function [elements,forces,springs] = read_data_file(data_file)
%READ_DATA_FILE 
% subfunction 1): extract information from text file, store
% it in structures.
%
    elements = cell(0);
    forces = cell(0);
    springs = cell(0);
    fid = fopen(data_file,'r');
    i = 0;
    j = 0;
    k = 0;
    while true
      thisline = fgetl(fid);
      if ~ischar(thisline); break; end  %end of file

      if isempty(thisline); continue; end

      % check if it is an element
      if strcmp(thisline,"ELEMENT")
          last = "ELEMENT";
          i = i+1;
          continue;
      end

      if strcmp(thisline,"FORCE")
          last = "FORCE";
          j = j+1;
          continue;
      end

      if strcmp(thisline,"SPRING")
          last = "SPRING";
          k = k+1;
          continue;
      end

      fields = split(thisline,':');
      command = fields{1};
      value = fields{2};

      if strcmp(last,"ELEMENT")
          switch command
              case "model"
                  elements{i}.model = strtrim(value);
              case "type"
                  elements{i}.type = strtrim(value);
              case "EA"
                  elements{i}.EA = str2func(strcat("@(x)",value));
              case "EJ"
                  elements{i}.EJ = str2func(strcat("@(x)",value));
              case "GJ"
                  elements{i}.GJ = str2func(strcat("@(x)",value));
              case "length"
                  elements{i}.length = str2double(value);
              case "shape functions"
                  elements{i}.shapes = str2func(strcat("@(x)",value));
              case "global index"
                  value = erase(value,{'[',']'});
                  value = strtrim(value);
                  scan = textscan(value, '%f','Delimiter',{' ',','});
%                   scan = scan{1}(~isnan(scan{1}));
                  elements{i}.index = scan{1};
              case "inclination"
                  elements{i}.alfa = str2double(value);
          end
      end

      if strcmp(last,"FORCE")
          switch command
              case "type"
                  forces{j}.type = strtrim(value);
              case "index"
                  forces{j}.index = str2double(value);
              case "element"
                  forces{j}.element = str2double(value);
              case "q"
                  forces{j}.q = str2func(strcat("@(x)",value));
              case "value"
                  forces{j}.f = str2double(value);
              
          end
      end

      if strcmp(last,"SPRING")
          switch command
              case "index"
                  springs{k}.index = str2double(value);
              case "k"
                  springs{k}.k = str2double(value);         
          end
      end

    end
    fclose(fid);

end

function [elements, forces] = build_local(elements, forces)
%BUILD_LOCAL 
% subfunction 2): use the structures to build local stiffness, local forces
% and local rotation matrices to global.
%    


    %for the elements which do not have shapes N build them
    for i = 1:length(elements)
        if isfield(elements{i},'shapes')
            continue;
        else
            l = elements{i}.length;
            switch elements{i}.model
                case "truss"
                    elements{i}.shapes = @(x) [(1-x/l); x/l];
                case "beam"
                    elements{i}.shapes = @(x) [1 - 3*(x/l)^2 + 2*(x/l)^3;
                                               l*(-x/l + 2*(x/l)^2 - (x/l)^3);
                                               3*(x/l)^2 - 2*(x/l)^3;
                                               l*((x/l)^2 - (x/l)^3)];
            end
        end
    end

    %after we have built all the shapes we generate K matrices

    for i = 1:length(elements)
        l = elements{i}.length;
        switch elements{i}.type
            case "axial"
                EA = elements{i}.EA;
                N = elements{i}.shapes;
                Nx = matlabFunction(diff(sym(N)));
                if nargin(Nx) ~= 0
                    K = integral(@(x) Nx(x)*EA(x)*Nx(x).',0,l,'ArrayValued',true);
                else
                    K = integral(@(x) Nx()*EA(x)*Nx().',0,l,'ArrayValued',true);
                end
                elements{i}.K = K;
            case "bending"
                EJ = elements{i}.EJ;
                N = elements{i}.shapes;
                Nxx = matlabFunction(diff(sym(N),2));
                if nargin(Nxx) ~= 0
                    K = integral(@(x) Nxx(x)*EJ(x)*Nxx(x).',0,l,'ArrayValued',true);
                else
                    K = integral(@(x) Nxx()*EJ(x)*Nxx().',0,l,'ArrayValued',true);
                end
                elements{i}.K = K;
            case "torsional"
                GJ = elements{i}.GJ;
                N = elements{i}.shapes;
                Nx = matlabFunction(diff(sym(N)));
                if nargin(Nx) ~= 0
                    K = integral(@(x) Nx(x)*GJ(x)*Nx(x).',0,l,'ArrayValued',true);
                else
                    K = integral(@(x) Nx()*GJ(x)*Nx().',0,l,'ArrayValued',true);
                end
                elements{i}.K = K;
        end

        if isfield(elements{i},'alfa')
            alfa = elements{i}.alfa;
            switch elements{i}.model
                case "truss"
                    elements{i}.M = [cosd(alfa), 0; sind(alfa), 0; 0 cosd(alfa); 0 sind(alfa)];
                case "beam"
                    elements{i}.M = [-sind(alfa) 0 0 0;
                                     cosd(alfa)  0 0 0;
                                              0  1 0 0;
                                     0 0 -sind(alfa) 0;
                                     0 0 -cosd(alfa) 0;
                                     0 0           0 1];
            end
        end
    end

    for j = 1:length(forces)
        if strcmp(forces{j}.type,"distributed")
            n = forces{j}.element;
            q = forces{j}.q;
            l = elements{n}.length;
            forces{j}.index = elements{n}.index;
            N = elements{n}.shapes;
            forces{j}.f = integral(@(x) N(x)*q(x),0,l,'ArrayValued',true);
        end
    end

end

function [K,F] = build_global(elements,forces,springs)
%BUILD_GLOBAL 
% subfunction 3): assemble all local informations into the global stiffness
% matrix and the global force vector.
%
    num_dof = 0;
    for i = 1:length(elements)
        n = max(elements{i}.index);
        if n > num_dof
            num_dof = n;
        end
    end

    K = zeros(num_dof);
    F = zeros(num_dof,1);

    for i = 1:length(elements)
        index_loc = 1:length(elements{i}.index);
        index_glb = elements{i}.index(elements{i}.index ~= 0);
        index_loc = index_loc(elements{i}.index ~= 0);

        if isfield(elements{i},'M')
            M = elements{i}.M;
            Ktr = M*elements{i}.K*M';
        else
            Ktr = elements{i}.K;
        end

        

        K(index_glb,index_glb) = K(index_glb,index_glb) + Ktr(index_loc,index_loc);
    end

    for j = 1:length(forces)
        

        if isfield(forces{j},'element')
            n_el = forces{j}.element;
            index_glb = elements{n_el}.index;
            
            if isfield(elements{n_el},'M')
                M = elements{n_el}.M;
                Ftr = M*forces{j}.f;
            else
                Ftr = forces{j}.f;
            end
        else
            index_glb = forces{j}.index;
            Ftr = forces{j}.f;
        end

        index_loc = 1:length(index_glb);
        index_glb = index_glb(index_glb ~= 0);
        index_loc = index_loc(index_glb ~= 0);

        F(index_glb,1) = F(index_glb,1) + Ftr(index_loc,1);

    end

    for k = 1:length(springs)
        index = springs{k}.index;
        K(index,index) = K(index,index) + springs{k}.k;
    end

end