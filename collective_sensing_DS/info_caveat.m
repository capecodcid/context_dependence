%first example

p=zeros(2,2);

 

for act = 1:49

    for bct = 1:49

        

p(1,1)=.01*act; %anything up to .5

p(1,2)=.5 - p(1,1);

p(2,1)=.01*bct; %anything up to .5

p(2,2)=1-p(1,1)-p(1,2)-p(2,1);
% index one is location of goat
% index 2 is guess
 

qa = p(1,:) + p(2,:);
%marginalize over location of goat -> just guesses

qa=qa/sum(qa);

qb = p(:,1) + p(:,2);
%marginalize over guess -> just goat

qb=qb/sum(qb);

 

info=0;

for ict=1:2

    for jct=1:2

        info = info + p(ict,jct)*log2(p(ict,jct)/(qa(jct)*qb(ict)));
        % KL divergence
    end

end

 

infoSAVE(act,bct) = real(info);

 

p_correct(act,bct) = p(1,1) + p(2,2);

 

    end

end

 

 

%% second example

 

X = rand(2,2);X = X./sum(X(:))  % we create a random probability matrix

ps = linspace(0,1,1000);

pTotal = sum(X(2,:));

for i=1:1000

Xtemp = X;

Xtemp(2,:) = pTotal.*[ps(i) 1-ps(i)];

infos(i)=sum(sum(Xtemp.*log(Xtemp./(repmat(sum(Xtemp),2,1).*repmat(sum(Xtemp,2),1,2)))));

end

plot(ps.*pTotal,infos,'o-')




