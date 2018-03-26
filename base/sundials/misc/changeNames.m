function [] = changeNames

ex_p={
    'pvanx'  'cvsadjnonx_p'
    'pvakx'  'cvsadjkryx_p'
    'pvfkx'  'cvsfwdkryx_p'
    'pvfnx'  'cvsfwdnonx_p'
    'pvkxb'  'cvskryx_bbd_p'  % this one MUST come before pvkx
    'pvkx'   'cvskryx_p'
    'pvnx'   'cvsnonx_p'
      };

ex_s={
    'cvdxe'   'cvsdenx_uw'    % this one MUST come before cvdx!!!
    'cvdx'    'cvsdenx'
    'cvbx'    'cvsbanx'
    'cvkxb'   'cvskryx_bp'         % this one MUST come before cvkx
    'cvkxt'   'cvskrydem_lin'  % this one MUST come before cvkx
    'cvkx'    'cvskryx'
    'cvdemd'  'cvsdirectdem'
    'cvdemk'  'cvskrydem_pre'
    'cvfnx'   'cvsfwdnonx'
    'cvfdx'   'cvsfwddenx'
    'cvfkx'   'cvsfwdkryx'
    'cvadx'   'cvsadjdenx'
    'cvabx'   'cvsadjbanx'
    'cvakxb'  'cvsadjkryx_pnt'     % this one MUST come before cvakx
    'cvakx'   'cvsadjkryx_int'
      };


% Process examples_ser

fprintf('PROCESS EXAMPLES_SER\n');
[rm_ser,add_ser] = process_files('examples_ser',ex_s);

% Process examples_par

fprintf('PROCESS EXAMPLES_PAR\n');
[rm_par, add_par] = process_files('examples_par',ex_p);

% Generate CVS commit log file

fid = fopen('cvs_msg','w');
fprintf(fid,'Changed example program names. Updated documentation files.\n');
fprintf(fid,'Serial examples\n');
for i = 1:size(ex_s,1)
  fprintf(fid,'%20s  ->  %20s\n',ex_s{i,1},ex_s{i,2});
end
fprintf(fid,'Parallel examples\n');
for i = 1:size(ex_p,1)
  fprintf(fid,'%20s  ->  %20s\n',ex_p{i,1},ex_p{i,2});
end
fclose(fid);

% Execute CVS commit

r = input('Commit examples directories to CVS? [y/n]  ','s');

if r ~= 'y'
  fprintf('\nOK, I will not run any CVS commands.\n');
  fprintf('I created a CVS log message file (''cvs_msg'') that you can use.\n');
  fprintf('Here are the commands you should run.\n');
end

cmd1 = sprintf('examples_ser/%s ',rm_ser{:});
cmd2 = sprintf('examples_par/%s ',rm_par{:});
cmd = sprintf('cvs rm %s %s',cmd1,cmd2);
if r == 'y'
  %system(cmd);
else
  fprintf('\n%s\n',cmd);
end

cmd1 = sprintf('examples_ser/%s ',add_ser{:});
cmd2 = sprintf('examples_par/%s ',add_par{:});
cmd = sprintf('cvs add %s',cmd1,cmd2);
if r == 'y'
  %system(cmd);
else
  fprintf('\n%s\n',cmd);
end
  
cmd = sprintf('cvs ci -F cvs_msg examples_ser examples_par');
if r == 'y'
  %system(cmd);
else
  fprintf('\n%s\n',cmd);
end
  



%
% ============================================
%

function [rm_files,add_files] = process_files(d,list)

nf = 0;
for i = 1:size(list,1)

  old = list{i,1};
  fprintf('  process %s -> ',old);
  new = list{i,2};
  fprintf('%s -> ',new);

% Create LaTeX name
  iu = strfind(new,'_');
  nu = length(iu); 
  if nu == 0
    newL = new;
  else
    chunk = new(1:iu(1)-1);
    newL = sprintf('%s_',chunk);
    for j = 2:length(iu)
      chunk = new(iu(j-1)+1:iu(j)-1);
      newL = sprintf('%s%s_',newL,chunk);
    end
    chunk = new(iu(nu)+1:end);
    newL = sprintf('%s%s',newL,chunk);
  end
  fprintf('%s\n',newL);

% Replace C file names
  cmd = sprintf('mv %s/%s.c %s/%s.c',d,old,d,new);
  %system(cmd);

  nf = nf+1;
  rm_files{nf} = sprintf('%s.c',old);
  add_files{nf} = sprintf('%s.c',new);
  
% Replace output file names
  outfiles = sprintf('%s/%s.out*',d,old);
  lfiles = dir(outfiles);
  for i = 1:length(lfiles)
    outfile = lfiles(i).name;
    ii = strfind(outfile,'.');
    ext = outfile(ii+1:end);
    cmd = sprintf('mv %s/%s %s/%s.%s',d,outfile,d,new,ext);
    %system(cmd);
    nf = nf+1;
    rm_files{nf} = lfiles(i).name;
    add_files{nf} = sprintf('%s.%s',new,ext);
  end
  
% Replace occurences in example files
  cmd = sprintf('sed ''s#%s#%s#g'' %s/%s.c > tmp',old,new,d,new);
  %system(cmd);
  cmd = sprintf('mv tmp %s/%s.c',d,new);
  %system(cmd);
  
% Replace occurences in Makefile.in
  cmd = sprintf('sed ''s#%s#%s#g'' %s/Makefile.in > tmp',old,new,d);
  %system(cmd);
  cmd = sprintf('mv tmp %s/Makefile.in',d);
  %system(cmd);
  
% Replace occurences in LaTeX files
  lfiles = dir('doc/*.tex');
  for i = 1:length(lfiles)
    cmd = sprintf('sed ''s#%s#%s#g'' doc/%s > tmp',old,newL,lfiles(i).name);
    %system(cmd);
    cmd = sprintf('mv tmp doc/%s',lfiles(i).name);
    %system(cmd);
  end
  
end

