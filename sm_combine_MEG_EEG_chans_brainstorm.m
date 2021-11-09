% This is a work around function to combine EEG and MEG channels to make default anatomy for SimMEEG
% these is probably a much simpler way of doing this in Brainstorm, but at the time I wasn't familiar enough with Brainstorm to figure it out.


%% Need to export channel information for MEG  as 'grad') and EEG as 'elec') from Brainstorm interface under subjects "common files" 

chans = grad; 

chans.Comment = '186 MEG (CTF) + 128 EEG (Easy Cap)';
chans.TransfEeg = elec.TransfEeg;
chans.TransfEegLabels = elec.TransfEegLabels;

for vv=1:length(elec.Channel)
    v = length(grad.Channel)+vv;
    chans.Channel(v).Loc = elec.Channel(vv).Loc;
    chans.Channel(v).Orient = elec.Channel(vv).Orient;
    chans.Channel(v).Comment = elec.Channel(vv).Comment;
    chans.Channel(v).Weight = elec.Channel(vv).Weight;
    chans.Channel(v).Type = elec.Channel(vv).Type;
    chans.Channel(v).Name = elec.Channel(vv).Name;
    
end

% then "import" in Brainstorm under subjects "common files" 


