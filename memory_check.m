function [mem_flag,free_mem]=memory_check(num_numbers,perc_remain)
% function [mem_flag]=memory_check(num_numbers)
% 
% Checks memory to make sure there is still at least 25% free after
% num_numbers will be created
%
% INPUT
%   num_numbers = number of numbers that will be created using nan(num_numbers) or zeros(num_numbers)
%   perc_remain = percentage of free memory that you still want remaining after reating num_numbers .
%
% OUTPUT
%   mem_flag = (0) insufficient memory (1) sufficient memory
% 
if nargin<2; perc_remain=10; end

[~,sys_mem]=memory;  free_mem=sys_mem.PhysicalMemory.Available; free_mem=free_mem*(1-(perc_remain/100)); % available memory in bytes minus 500 MBytes
if free_mem<0
%     fprintf('ERROR! Not enough free memory. Please free up memory at least %.f MegaBytes. \nProgram Terminated.\n',free_mem/1e6);
    mem_flag=0;
end
needed_mem=num_numbers*8;  % 8bytse per number
mem_diff=(needed_mem-free_mem);
if mem_diff>=0
    allow_cont=mem_diff/(num_numbers*8);
%     fprintf('ERROR! Not enough free memory.\nPlease free up memory at least %.f MegaBytes. Reduce number of numbers to %.f\nProgram Terminated.\n',mem_diff/1e6,allow_cont);
    mem_flag=0;
else
%     fprintf('Sufficient Memory.  Available Memory(%.f MB) - Needed Memory(%.f MB)  =  Remaining(%.f MB)\n',sys_mem.PhysicalMemory.Available/1e6,  needed_mem/1e6,(sys_mem.PhysicalMemory.Available-needed_mem)/1e6);
    mem_flag=1;
end
