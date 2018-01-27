function DrawPts(xp,nPts,texture)
% draw particles on screen
tot_n = 0;
col = 'yellow';
for n = 1:length(nPts)
    switch texture(n,:)
        case 'Jello'
            col = 'red';
        case 'Snow'
            col = 'green';
        case 'Water'
            col = 'blue';
    end
    scatter(xp(1,tot_n+1:tot_n+nPts(n)),xp(2,tot_n+1:tot_n+nPts(n)),5.,col);
    hold on;
    tot_n = tot_n + nPts(n);
end
