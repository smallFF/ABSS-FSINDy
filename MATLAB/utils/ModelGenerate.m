classdef ModelGenerate < handle
    properties
        Charset;
        Xi;
        func_name;
        ode_model;
    end

    properties (Hidden)
        filename;
        fw;
        func_file;
    end
    
    methods (Static)
        function r = format_transform(item)
            if item == '1'
                r = item;
            else
                r = [];
                str_splited = strsplit(item, 'x');
                for i = 2:length(str_splited)
                    if i == length(str_splited)
                        r = [r, 'X(',num2str(str_splited{i}), ')'];
                    else
                        r = [r, 'X(',num2str(str_splited{i}), ')*'];
                    end
                end
            end
        end
    end
    methods
        function obj = ModelGenerate(Charset, Xi, func_name)
            obj.Charset = Charset;
            obj.Xi = Xi;
            obj.func_name = func_name;
            obj.filename = [func_name, '.m'];
            obj.initial_function_file();
            obj.write_to_file();
            obj.fw.delete();
        end
        
        function obj = initial_function_file(obj)
            obj.fw = FileWriter(obj.filename);
            output_params = 'dX';
            input_params = 't, X';
            obj.func_file = WriteFunctionFile(obj.fw, obj.func_name, output_params, input_params);
        end
        
        function obj = write_to_file(obj)
            % contents
            dim = size(obj.Xi, 2);
            
            obj.func_file.write(['dX=zeros(',num2str(dim),',1);']);
            tmp_model = cell(dim, 1);
            nzind = obj.Xi~=0;
            for i = 1:dim
                expr = ['dX(',num2str(i),')='];
                nz_count = sum(nzind(:, i));
                Charset_i = obj.Charset(nzind(:,i));
                Xi_i    = obj.Xi(nzind(:,i),i);
                for j = 1:nz_count
                   if Xi_i(j) > 0
                       flag = '+';
                   else
                       flag = '-';
                   end
                   
                   formated_Charset_item = ModelGenerate.format_transform(Charset_i{j});
                   % disp('formated_Charset_item')
                   % disp(formated_Charset_item);
                   
                   item = [num2str(abs(Xi_i(j))), '*', formated_Charset_item]; 
                   % disp('item')
                   % disp(item);
                   if j == 1
                        if flag == '-'
                            expr = [expr, flag, item];
                        else
                            expr = [expr, item];
                        end
                   else
                       expr = [expr, flag, item];
                   end
                end
                tmp_model{i} = expr;
                
                obj.func_file.write([expr, ';']);
            %     disp(model{i});
            end
            
            % end
            obj.ode_model = tmp_model;
            
        end
    end
end