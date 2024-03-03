function Answer=dispInput(Menu,Menu_GUI,titl,Default_Answer)
disp(['----- '  titl ' ----- '])
Rows=size(char(Menu),1);
Answer=cell(Rows,1); % столбец
if Menu_GUI
   numlines=1; 
   options.Resize='on'; % 'on' 'off'
   options.WindowStyle='modal' ; % 'normal' 'modal'
   options.Interpreter='tex';
   Answer=inputdlg(Menu,titl,numlines,Default_Answer,options);  
   if isempty(Answer)
%          Answer=cell2mat(Default_Answer);
   return
   end
   for i=1:Rows
        disp( [Menu{i} ' =  ' strvcat(Answer{i})]);  
   end
   Answer=str2num(strvcat(Answer));  
else   
     for i=1:Rows
        Answer{i}=input([Menu{i} ' = ']);         
     end
     Answer=cell2mat(Answer);
%      Answer=str2mat(Answer);
end
% if isempty( Answer)
%     Answer=cell2mat(Default_Answer);
% end 
disp(' ')