%Sector ETF Momentum Strategy
clear;

ETF = {'XLB','XLE','XLF','XLI','XLK','XLP','XLU','XLV','XLY','SPY'};
for i=1:length(ETF)
    data = fetch(yahoo, ETF(i,1),{'Open','Close'}, '01/01/2000', today);
    Shares.Name(1,i) = ETF(i,1);
    Shares.Open(1:size(data,1),i) = data(:,2);
    Shares.Close(1:size(data,1),i) = data(:,3);
    Shares.Date(1:size(data,1),i) = data(:,1);
end
close(y);
Shares.Open = flipud(Shares.Open);
Shares.Date = flipud(Shares.Date);
Shares.Close = flipud(Shares.Close);
save ETF Shares
