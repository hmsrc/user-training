function [count,index]=searchit(mystring,mytarget)
  
 % this function returns the total count of mytarget
 index=strfind(mystring,mytarget);  

 number=length(index);
 count=number;
end