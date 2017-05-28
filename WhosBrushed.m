function Names = WhosBrushed()
% Names = WhosBrushed
%
% goes through all children of current plot, finds which objects are
% brushed, and returns names stored in UserData

ch = get(gca, 'Children');
Names = {};
for i=1:length(ch)
    if ~(isprop(ch(i), 'brushdata')  & isprop(ch(i), 'UserData') )
        continue;
    end
    brushed = get(ch(i), 'brushdata');
    MyNames = get(ch(i), 'UserData');
    BrushedNames = MyNames(find(brushed));
    Names = vertcat(Names, BrushedNames(:));
end