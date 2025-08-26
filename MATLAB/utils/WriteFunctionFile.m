classdef WriteFunctionFile < handle
   properties
       fw;
       func_name;
       output_params;
       input_params;
   end
   
   methods
       function obj = WriteFunctionFile(FileWriter_obj, func_name, output_params, input_params)
          obj.fw = FileWriter_obj;
          obj.func_name = func_name;
          obj.output_params = output_params;
          obj.input_params = input_params;
          obj.write_head_line();
       end
       
       function write_head_line(obj)
          head_line = sprintf('function [%s]=%s(%s)', obj.output_params, obj.func_name, obj.input_params);
          obj.write(head_line); 
       end
       
       function write(obj, text_str)
           obj.fw.write_to_file(text_str);
       end
   end
end