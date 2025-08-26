classdef FileWriter < handle
   properties (Access = private)
      file_ID; 
   end
   
   methods
       function obj = FileWriter(filename)
          obj.file_ID = fopen(filename, 'w'); 
       end
       
       function write_to_file(obj, text_str)
          fprintf(obj.file_ID, '%s\n', text_str); 
       end
       
       function delete(obj)
          fclose(obj.file_ID); 
       end
       
   end
end