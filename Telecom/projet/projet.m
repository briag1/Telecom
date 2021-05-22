%2
%information binaire a transmettre
    Fe=24000;
    Te=1/Fe;
    Rb=3000;
    %longueur du message en bits
    nb_bits=10000;
    %generation de l'information binaire a transmettre
    bits=randi([0 1],1,nb_bits);
 %Modulateur/demodulateur1
    %nombre de symbole possibles pour le modulateur1 
    M=2;
    %debit symbole
    Rs=Rb/log2(M);
    %Duree symbole en nombre d'echantillons (Ts=NsTe)
    Ns=Fe/Rs;
    %reponse_impulsionnellle_globale 
    %modulateur1
    h1=ones(1,Ns);
    %Mapping binaire a? moyenne nulle : 0->-1, 1->1 
    Symboles=2*bits-1;
    %Surechantillionnage modulateur1
    suite_diracs=kron(Symboles,[1 , zeros(1,Ns-1)]);
    %filtrage
    X1=filter(h1,1, suite_diracs);
    %signal recu
    X1r=filter(h1,1,X1);
    %echantillonnage
    n0=Ns;
    X1e=X1r(n0:Ns:length(X1r));
    %decisions /demaping
    X1e(X1e>0)=1;
    X1e(X1e<0)=0;
    %taux d'erreur binnaire
    erreur=mean(abs(X1e - bits))
%% 3
    %sequence binaire
    bits=[1 1 1 0 0 1 0];
    %Modulateur/demodulateur1
    %nombre de symbole possibles pour le modulateur1 
    M=2;
    %debit symbole
    Rs=Rb/log2(M);
    %Duree symbole en nombre d'echantillons (Ts=NsTe)
    Ns=Fe/Rs;
    %reponse_impulsionnellle_globale 
    %modulateur1
    h1=ones(1,Ns);
    %Mapping binaire a? moyenne nulle : 0->-1, 1->1 
    Symboles=2*bits-1;
    %Surechantillionnage modulateur1
    suite_diracs=kron(Symboles,[1 , zeros(1,Ns-1)]);
    g1=conv(h1,h1);
    g2=cat(2, zeros(1,Ns),1/2*g1);
    g1l=cat(2,g1,zeros(1,Ns));
    g=g1l+g2;
    figure;
    plot(1:length(g1l),g1l,1:length(g1l),g2);
    figure;
    plot(g);
    
    Xr=filter(g,1,suite_diracs);
    figure;
    plot(Xr);
    
    figure;
    plot(reshape(Xr,Ns,length(Xr)/Ns));
    
    %echantillionnage
    Xe=Xr(Ns:Ns:end);
    
    %decisions/demapping 
    Xe(Xe>0)=1;
    Xe(Xe<0)=0;
    TEB=mean(abs(Xe-bits))
    %% 4
    %sequence binaire
   %longueur du message en bits
    nb_bits=10000;
    %generation de l'information binaire a transmettre
    bits=randi([0 1],1,nb_bits);
    %Modulateur/demodulateur1
    %nombre de symbole possibles pour le modulateur1 
    M=2;
    %debit symbole
    Rs=Rb/log2(M);
    %Duree symbole en nombre d'echantillons (Ts=NsTe)
    Ns=Fe/Rs;
    %reponse_impulsionnellle_globale 
    %modulateur1
    h1=ones(1,Ns);
    %Mapping binaire a? moyenne nulle : 0->-1, 1->1 
    Symboles=2*bits-1;
    
    %Surechantillionnage modulateur1
    suite_diracs=kron(Symboles,[1 , zeros(1,Ns-1)]);
    X=conv(suite_diracs,h1,'same');
    g1=conv(h1,h1);
    g2=cat(2, zeros(1,Ns),1/2*g1);
    g1l=cat(2,g1,zeros(1,Ns));
    
    
    Xr=filter(g,1,suite_diracs);
    Px=mean(abs(X).^2);
    rapportEN=linspace(1,6,100)';
    sigma=sqrt(Px*Ns./(2*log2(M).*rapportEN));
    bruitR=sigma*randn(1,length(X));
    bruitfiltre=zeros(size(bruitR));
   
    bruitfiltre(1,:)=conv(bruitR(1,:),h1,'same');
    bruitfiltre(2,:)=conv(bruitR(2,:),h1,'same');
    for i=1:length(rapportEN)
        bruitfiltre(i,:)=conv(bruitR(i,:),h1,'same');
    end
    Xr=repmat(Xr,length(rapportEN),1);
    Xbr=Xr+bruitfiltre;
    %echantillionnage
    Xe=Xbr(:,Ns:Ns:end);
    
    %decisions/demapping 
    Xe(Xe>0)=1;
    Xe(Xe<0)=0;
    TEB=mean(abs(Xe-bits),2);
    figure;
    semilogy(10*log10(rapportEN),TEB);
    
    
    
    %% 3.2
    %sequence binaire
   %longueur du message en bits
    nb_bits=10000;
    %generation de l'information binaire a transmettre
    bits=randi([0 1],1,nb_bits);
    bits=zeros(1,length(bits));
   bits(1)=1;
    %Modulateur/demodulateur1
    %nombre de symbole possibles pour le modulateur1 
    M=2;
    %debit symbole
    Rs=Rb/log2(M);
    %Duree symbole en nombre d'echantillons (Ts=NsTe)
    Ns=Fe/Rs;
    %reponse_impulsionnellle_globale 
    %modulateur1
    h1=ones(1,Ns);
    %Mapping binaire a? moyenne nulle : 0->-1, 1->1 
    Symboles=2*bits-1;
    %Surechantillionnage modulateur1
    suite_diracs=kron(Symboles,[1 , zeros(1,Ns-1)]);
    g1=conv(h1,h1);
    g2=cat(2, zeros(1,Ns),1/2*g1);
    g1l=cat(2,g1,zeros(1,Ns));
    g=g1l+g2;
    Xr=filter(g,1,suite_diracs);
    %echantillionnage
    Xe=Xr(Ns:Ns:end);
    %decisions/demapping 
    Xe(Xe>0)=1;
    Xe(Xe<0)=0;
    TEB=mean(abs(Xe-bits))
    Z=zeros(length(bits),2);
    no=Ns;
    for i=1:length(bits)
        for j=1:min(i,2)
            
            Z(i,j)=Xr(no+Ns*(i-j));
        end
    end
    Y=zeros(length(bits),1);
    Y(1)=1;
    C=(Z'*Z)\(Z'*Y);
    