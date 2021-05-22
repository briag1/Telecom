clear all;

Fe=24000;
Rb=6000;
%longueur du message en bits
nb_bits=1000;

%generation de l'information binaire a transmettre
bits=randi([0 1],1,nb_bits);
  

%Modulateur1
    %nombre de symbole possibles
    M=2;
    %debit symbole
    Rs=Rb/log2(M); % ici Rs=Rb
    %Dure?e symbole en nombre d'echantillons (Ts=NsTe)
    Ns=Fe/Rs;
    h1=ones(1,Ns);
    %Mapping binaire a? moyenne nulle : 0->-1, 1->1 
    Symboles=2*bits-1;
    %Surechantillionnage modulateur1
    suite_diracs=kron(Symboles,[1 , zeros(1,Ns-1)]);
    %filtrage
    X1=filter(h1,1, suite_diracs);
 
%Modulateur2
    %nombre de symbole possibles
    M2=4;
    %debit symbole
    Rs2=Rb/log2(M2);
    %Dure?e symbole en nombre d'echantillons (Ts=NsTe)
    Ns2=Fe/Rs2;
    h2=ones(1,Ns2);
    %surechantillonnage
    Symboles2=bits;
    for i=1:2:nb_bits
        if Symboles2(i)==0 && Symboles2(i+1)==0
            Symboles2(i)=-3;
            Symboles2(i+1)=-3;
        elseif Symboles2(i)==0 && Symboles2(i+1)==1
            Symboles2(i)=-1;
            Symboles2(i+1)=-1;
        elseif Symboles2(i)==1 && Symboles2(i+1)==1
            Symboles2(i)=1;
            Symboles2(i+1)=1;
        else
            Symboles2(i)=3;
            Symboles2(i+1)=3;
        end
    end
    suite_diracs2=kron(Symboles2,[1 , zeros(1,Ns2-1)]);
    X2=filter(h2, 1, suite_diracs2);   
    
%Modulateur3
    %nombre de symbole possibles
    M3=2;
    %debit symbole
    Rs3=Rb/log2(M3);
    %Dure?e symbole en nombre d'echantillons (Ts=NsTe)
    Ns3=Fe/Rs3;
    %h3=cat(2,ones(1,Ns/2),-ones(1,Ns/2),2);
    h3=[1,1,-1,-1];
     %Mapping binaire a? moyenne nulle : 0->-1, 1->1 
    Symboles3=2*bits-1;
    %Surechantillionnage modulateur1
    suite_diracs3=kron(Symboles3,[1 , zeros(1,Ns3-1)]);
    %filtrage
    X3=filter(h3,1, suite_diracs3);
    
    
    
%Modulateur4
    %coeff de roaloff 
    alpha=0.5;
    %nombre de symbole possibles
    M4=2;
    %debit symbole
    Rs4=Rb/log2(M4);
    %Dure?e symbole en nombre d'echantillons (Ts=NsTe)
    Ns4=Fe/Rs4;
    %h3=cat(2,ones(1,Ns/2),-ones(1,Ns/2),2);
    h4=rcosdesign(alpha,8,Ns4);
    %surechantillionnage
    Symboles4=2*bits-1;
    suite_diracs4=kron(Symboles4,[1 , zeros(1,Ns4-1)]);
    %filtrage
    X4=filter(h4,1,suite_diracs4);
 
%trace des siganaux obtenues
    % signa1 X1
        Te=1/Fe;
        fs1=figure;
        t1=linspace(0,Ns*nb_bits*Te,length(X1));
        plot(t1,X1);
        xlabel('s');
        ylabel('X1(t)');
        title('Signal X4');
    %signal X2
        fs2=figure;
        t2=linspace(0,Ns*nb_bits*Te,length(X2));
        plot(t2,X2);
        xlabel('s')
        xlabel('X2(t)');
    %signal X3
        fs3=figure;
        t3=linspace(0,Ns*nb_bits*Te,length(X3));
        plot(t3,X3);
        xlabel('s')
        xlabel('X3(t)');
        title('Signal X3');
    %signal x4
        fs4=figure;
        t4=linspace(0,Ns*nb_bits*Te,length(X4));
        plot(t4,X4);
        xlabel('s')
        xlabel('X4(t)');
        title('Signal X4');
%densité spectrale depuissance des  signaux
  
    %signal x1
        %
        PSD=pwelch(X1,[],[],[],Fe,'twosided')';
        f1=linspace(0,Fe,length(PSD));
        fig1=figure;
        %figure("psd signal 1");
        %plot(f,PSD);
        % densite spectrale de puissance theorique
        Te=1/Fe;
        Ts=Ns*Te;  
    
        psdth=Ts*sinc(f1*Ts).^2;
        plot(f1,PSD,f1,psdth);
       title('densite spectrale de puissance du signal X1');
        xlabel('f')
        legend('psd','psdth');
    %signal x2
        PSD2= pwelch(X2,[],[],[],Fe,'twosided')';
        f2=linspace(0,Fe,length(PSD2));
        fig2=figure;
    
        Ts2=Ns2*Te;
        psdth2=5*Ts2*((sin(pi*f2*Ts2))./(pi*f2*(Ts2))).^2;
        plot(f2,PSD2,f2,psdth2);
        xlabel('f');
        title('densite spectrale de puissance du signal X2');
        legend('PSD2','psdth2')
    %signal x3
        Ts3=Ns3*Te;
        PSD3=pwelch(X3,[],[],[],Fe,'twosided')';
        f3=linspace(0,Fe,length(PSD3));
        psd3th=Ts3*sin(pi*f3*Ts/2).^4./(pi*f3*Ts/2).^2;
        fig3=figure;
        plot(f3,PSD3,f3, psd3th);
        title('densite spectral de puissance du signal X3');
        xlabel('f');
        legend('PSD3','psdth3');
        
        
    
    
    %signal x4
        fig4=figure;
        alpha=0.5;
        PSD4=pwelch(X4,[],[],[],Fe,'twosided')';
        f4=linspace(0,Fe,length(PSD4));
        psd4th=zeros(1,length(PSD4));
        i=find(abs(f4)<(1-alpha)/(2*Ts));
        j=find((abs(f4)>((1-alpha)/(2*Ts))) & (abs(f4)<((1+alpha)/(2*Ts))));
        psd4th(i)=10^(-4);
        psd4th(j)=10^(-4)/2*(1+cos(pi*((Ts/alpha)*(abs(f4(j))-(1-alpha)/(2*Ts)))));
        plot(f4,PSD4, f4,psd4th);
        xlabel('f');
        title('densite spectral de puissance du signal X4');
        legend('PSD4','psdth4');