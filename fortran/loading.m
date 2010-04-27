for i= 1:41
% strk = ['sk',int2str(i),'.mat'];
stra = ['sa',int2str(i),'.mat'];
prefixa = 'aem/alex_program/';
% prefixk = 'aem/matlab/point_source/';
fa = ['sa',int2str(i)];
% fk = ['sk',int2str(i)];
commanda = ['load ',prefixa,stra,';'];
% commandk = ['load ',prefixk,strk,';'];
% eval(commandk);
% eval([fk,'= coeff;']);
eval(commanda);
eval([fa,'= out;']);
end
