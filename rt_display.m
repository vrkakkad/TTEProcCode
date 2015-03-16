target_dir = 'C:\Users\vrk4\Tools\Transthoracic_Clinical\Data\Temp';
cd(target_dir)

orig = size(dir('arfi_par*'),1);

while 1
    
    new = size(dir('arfi_par*'),1);
    
    if new-orig>0
        pause(5)
        fprintf(1,'Files Detected...')
        runTTE_raw(pwd,-1)
        orig = size(dir('arfi_par*'),1);
    end
end