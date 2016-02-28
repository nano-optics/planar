%%%%                testing 1D grating against APL paper
%%%%                                    
%%%% Le réseau peut être 1D (périodique dans une direction) ou 2D (périodique dans 2 directions)
%%%% Le réseau est supposé sub-longueur d'onde (période < longueur d'onde), seul l'ordre de diffraction (0,0) 
%%%% 
%%%% Grandeurs calculées en fonction de la longueur d'onde : 
%%%% - Coefficient de réflexion dans l'ordre 0
%%%% - Coefficient de transmission dans l'ordre 0 
%%%% 
%%%% On peut également calculer pour une longueur d'onde
%%%% - Une coupe du champ électromagnétique (au choix parmi les 3 possibilités x=x0, y=y0 ou z=z0)
%%%%
%%%% On fait le calcul soit avec RCWA standard (op_granet=0) soit avec RCWA modifiée (op_granet=1) pour améliorer la convergence
%%%% Cette option peut être utile pour les réseaux métalliques. A priori pas nécessaire pour les réseaux diélectriques
%%%%
%%%% Unité de longueur : µm  (la longueur d'onde et toutes les longueurs sont en µm)
%%%%
%
%
% Vue de dessus du réseau
% -----------------------------
%                              periodicity_x
%                             <------------->                                                                             
%         *******       *******       *******       *******       ******* ^
%         *******       *******       *******       *******       ******* |                                                                  
%         *******       *******       *******       *******       ******* | periodicity_y                                                             
%                                                                         |                                                                                                                                         
%                                                                         |                                
%         *******       *******       *******       *******       ******* v                                                                   
%         *******       **ngm**       *******       *******       *******                                                                    
%         *******       *******       *******       *******       *******
%                                ng                                                                                                                                                                                  
%                                                                                                         
%         *******       *******       *******       *******       ******* ^                                                                   
% y ^     *******       *******       *******       *******       ******* | diameter_y                                                                  
%   |     *******       *******       *******       *******       ******* v                                                                  
%   |                                 <----->
%   --->  x                           diameter_x
%                                                                                                                                                                                         
%
% Sorties du code (vecteurs de la même taille que wavelength, sauf les champs)
% ------------------------------------------------------------------------------------------------------------------------------------
% On calcule les coefficients de réflexion dans l'ordre (0,0) uniquement (réseau supposé sub-longueur d'onde, ordres supérieurs évanescents)
%
% ref_TE_TE : coeff de réflexion en amplitude pour une onde incidente TE et une onde réfléchie TE
% ref_TE_TM : coeff de réflexion en amplitude pour une onde incidente TE et une onde réfléchie TM
% ref_TM_TM : coeff de réflexion en amplitude pour une onde incidente TM et une onde réfléchie TM
% ref_TM_TE : coeff de réflexion en amplitude pour une onde incidente TM et une onde réfléchie TE
% R0_TE_TE : coeff de réflexion en intensité pour TE --> TE 
% R0_TE_TM : coeff de réflexion en intensité pour TE --> TM 
% R0_TM_TM : coeff de réflexion en intensité pour TM --> TM
% R0_TM_TE : coeff de réflexion en intensité pour TM --> TE
%
% Note : Il n'y a que quand theta(1)~=0 ET theta(2)~=0 que les 4 coefficients de réflexion sont non nuls. 
% Sinon, dès que l'un des angles theta est égal à zéro, le programme utilise la symmétrie du problème et un seul coefficient est non nul.
% La polarisation de l'onde incidente est dans ce cas fixée avec la variable pol (cf paramètres numériques)
%
%
% Si trace_champ=1  (dans ce cas, une seule longueur d'onde)
% ------------------------------------------------------------
% Ex,Ey,Ez : coupe du champ E aux points zz,xx,yy 
% Hx,Hy,Hz : coupe du champ H aux points zz,xx,yy
% On choisit une coupe parmi x=x0, y=y0 ou z=z0

close all;
clear;retio;
[prv,vmax]=retio([],inf*i);             % PAS écriture sur fichiers (inutile sur des versions récentes de Matlab)

%%%%%% Longueur d'onde et angle d'incidence
wavelength=linspace(0.8,1.1,200); wavelength = 0.95;%wavelength = 0.7424;
theta=[0,0];                            % angle d'incidence en degrés
lambda0 = 0.9;                          % DBR central wvl

%%%%%% Indices de réfracion
nh=1;                                   % indice du milieu incident
ng=nh;                                   % indice entre les plots du réseau
ngm=retindice(wavelength,2.7);%ngm=nh;            % indice des plots du réseau (retindice 2.7=Au JC)
nm=3.5;             % indice de la couche entre le Bragg et le réseau
n1=3.5;            % indice de la couche 1 du miroir de Bragg 
n2=3;           % indice de la couche 2 du miroir de Bragg
nsub=n1;                               % indice du substrat 

%%%%%% Paramètres géométriques
periodicity_x=0.300;                     % période en x
periodicity_y=periodicity_x;            % période en y
diameter_x=periodicity_x;           % taille des plots en x 
diameter_y=periodicity_y;               % taille des plots en y (diameter_y=periodicity_y si réseau 1D)
hg=0.05;hg=0.02;                                % épaisseur du réseau
hm=0.1; hm=0.00;                                % épaisseur de la couche entre le Bragg et le réseau
h1=lambda0/4/n1;%h1=0;                        % épaisseur de la couche 1 du miroir de Bragg
h2=lambda0/4/n2;%h2=0;                        % épaisseur de la couche 2 du miroir de Bragg
nrep=15;                                 % nombre de périodes dans le miroir de Bragg (2xnrep couches)

qx = nh*sin(theta(1)*pi/180);
qy = nh*sin(theta(2)*pi/180);

%%%%%% Paramètres numériques

sens=2;                              % sens=1 calcul du champ incident du haut
                                     % sens=2 calcul du champ incident du
                                     % substrat                    
                                     
pol=2;                               % polarisation incidente, TM pol=2  TE pol=0. TE (TM) = champ E (champ H) Transverse au plan d'incidence
                                     % En incidence normale, TM veut dire H//y et TE est E//y
Mx=1;                               % nombre de termes de Fourier en x
My=0;                                % nombre de termes de Fourier en y (My=0 si réseau 1D)
op_granet=0;                         % si 1, RCWA modifiée pour améliorer la convergence (transformée de coordonnées réelles sur les discontinuités) 
if(length(wavelength) == 1)
trace_champ=1;                       % si 1, calcul d'une coupe du champ
                                     % ATTENTION: dans ce cas seulement une
                                     % longueur d'onde     
else trace_champ=0;
end

if sens == 2                         % in-plane wavevector corresponding to incidence from the substrate 
    qx = qx * nsub /nh;
    qy = qy * nsub /nh;
end

x0=[];                                % coupe en x=x0 si trace_champ=1 ([] si on ne fait pas de coupe en x0)
y0=0;                                % coupe en y=y0 si trace_champ=1 ([] si on ne fait pas de coupe en y0)                       
z0=[];                                % coupe en z=z0 si trace_champ=1 ([] si on ne fait pas de coupe en z0)
                                     % z=0 est situé tout en bas de la structure, à une hauteur h_sub dans le substrat
h_air=0.05;%h_air=2*wavelength;                          % épaisseur au-dessus pour tracé du champ (trace_champ=1)
h_sub=0.05;%h_sub=2*wavelength;                          % épaisseur dans le substrat pour tracé du champ (trace_champ=1)
op_objet=0;                          % si 1, tracé de la géométrie pour vérifier la structure calculée

if My==0
    if op_granet==1;Bx=500;Ax=0.02/Bx;By=[];Ay=[];xdisc=[-diameter_x/2,diameter_x/2];ydisc=[];end;
else
    if op_granet==1;Bx=500;Ax=0.02/Bx;By=Bx;Ay=Ax;xdisc=[-diameter_x/2,diameter_x/2];ydisc=[-diameter_y/2,diameter_y/2];end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_TE_TE=zeros(1,length(wavelength));ref_TE_TM=ref_TE_TE;ref_TM_TM=ref_TE_TE;ref_TM_TE=ref_TE_TE;
refsub_TE_TE=zeros(1,length(wavelength));refsub_TE_TM=refsub_TE_TE;refsub_TM_TM=refsub_TE_TE;refsub_TM_TE=refsub_TE_TE;
trans_TE_TE=zeros(1,length(wavelength));trans_TE_TM=trans_TE_TE;trans_TM_TM=trans_TE_TE;trans_TM_TE=trans_TE_TE;
R0_TE_TE=zeros(1,length(wavelength));R0_TE_TM=R0_TE_TE;R0_TM_TM=R0_TE_TE;R0_TM_TE=R0_TE_TE;
R0sub_TE_TE=zeros(1,length(wavelength));R0sub_TE_TM=R0sub_TE_TE;R0sub_TM_TM=R0sub_TE_TE;R0sub_TM_TE=R0sub_TE_TE;
T0_TE_TE=zeros(1,length(wavelength));T0_TE_TM=T0_TE_TE;T0_TM_TM=T0_TE_TE;T0_TM_TE=T0_TE_TE;
Ex=[];Ey=[];Ez=[];Hx=[];Hy=[];Hz=[];Einc=[];Hinc=[];
x=[];y=[];z=[];indice=[];
Ntre=1;
if trace_champ==1;op_retcouche=1;else op_retcouche=0;end;
if trace_champ==1&&isempty(x0)==1&&isempty(y0)==1&&isempty(z0)==1;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
if trace_champ==1&&isempty(x0)==0&&isempty(y0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
if trace_champ==1&&isempty(x0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
if trace_champ==1&&isempty(y0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;

for zou=1:length(wavelength)
    lwa=wavelength(zou);
    k0=2*pi/lwa;
    tic
    period=[periodicity_x,periodicity_y];

    if length(ng)==1;nng=ng;else nng=ng(zou);end;
    if length(ngm)==1;nngm=ngm;else nngm=ngm(zou);end;
    if length(nm)==1;nnm=nm;else nnm=nm(zou);end;
    if length(n1)==1;nn1=n1;else nn1=n1(zou);end;
    if length(n2)==1;nn2=n2;else nn2=n2(zou);end;
    if length(nsub)==1;ns=nsub;else ns=nsub(zou);end;

    if theta(1)==0 && theta(2)==0;sym=[pol-1,pol-1,0,0];end; 
    if theta(1)~=0 && theta(2)==0;sym=[0,pol-1,0,0];end; 
    if theta(1)==0 && theta(2)~=0;sym=[1-pol,0,0,0];end;
    if theta(1)~=0 && theta(2)~=0;sym=[];end;

    kx=k0*qx;
    ky=k0*qy;
    beta=[kx,ky];
    
    uh=retu(period,{nh,k0});                                    
    ug=retu(period,{nng,[0,0,diameter_x,diameter_y,nngm,Ntre],k0});
    um=retu(period,{nnm,k0});
    u1=retu(period,{nn1,k0});
    u2=retu(period,{nn2,k0});
    ub=retu(period,{ns,k0});
    
    if op_granet==1;
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym,{[],{xdisc,Ax,Bx,ydisc,Ay,By}});ah=retcouche(init,uh,1);
    else
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym);ah=retcouche(init,uh,op_retcouche);
    end
    ag=retcouche(init,ug,op_retcouche);
    am=retcouche(init,um,op_retcouche);
    a1=retcouche(init,u1,op_retcouche);
    a2=retcouche(init,u2,op_retcouche);
    ab=retcouche(init,ub,op_retcouche);
        
    if op_objet==1
        struc_test=cell(1,2*nrep+4);
        struct_test{1}={0.1,uh};
        struct_test{2}={hg,ug};
        struct_test{3}={hm,um};
        for az=0:1:nrep-1
           struct_test{(2*az+1)+3}={h1,u1};
           struct_test{(2*az+1)+4}={h2,u2};
        end
        struct_test{2*nrep+4}={0.05,ub};
        rettestobjet(period,struct_test,[0,-1.5,1],[2,2,-period/2]);
    end

    if op_granet==1
        init_pasgranet=retinit(period,[-Mx,Mx,-My,My],beta,sym);
        ah_pasgranet=retcouche(init_pasgranet,uh,1);
        ab_pasgranet=retcouche(init_pasgranet,ub,1);
        sh_passage=retpassage(init,ah_pasgranet,ah,1);
        sb_passage=retpassage(init,ab_pasgranet,ab,-1);
        [sh,haut,beth,ch,anglh]=retb(init_pasgranet,ah_pasgranet,1e-4);
        [sb,bas,betb,cb,anglb]=retb(init_pasgranet,ab_pasgranet,-1e-4);
        sh=retss(sh,sh_passage,retc(ah,2*lwa));
        sb=retss(retc(ab,2*lwa),sb_passage,sb);
    else
        [sh,haut,beth,ch,anglh]=retb(init,ah,1e-4);
        [sb,bas,betb,cb,anglb]=retb(init,ab,-1e-4);
    end
    ordre=[0,0];
    K=[kx,ky]+ordre*2*pi./period;
    inch=find(((K(1)-beth(:,1)).^2+(K(2)-beth(:,2)).^2)<1e-8);
    difh=inch;
    incb=find(((K(1)-betb(:,1)).^2+(K(2)-betb(:,2)).^2)<1e-8);
    difb=incb;
    sh=rettronc(sh,haut(inch,1),haut(difh,1),1);
    sb=rettronc(sb,bas(incb,1),bas(difb,1),-1);
    
    sg=retc(ag,hg);
    sm=retc(am,hm);
    s1=retc(a1,h1);
    s2=retc(a2,h2);
    stemp=retss(sg,sm);
    for az=1:nrep
        stemp=retss(stemp,s1,s2);
    end
    stot=retss(sh,stemp,sb);
    ef=retreseau(init,stot,betb,cb,anglb,incb,difb,beth,ch,anglh,inch,difh);
    if isempty(ef.amplitude(ef.dif.TEh,ef.inc.TEh))==0;ref_TE_TE(zou)=ef.amplitude(ef.dif.TEh,ef.inc.TEh);end;
    if isempty(ef.amplitude(ef.dif.TMh,ef.inc.TEh))==0;ref_TE_TM(zou)=ef.amplitude(ef.dif.TMh,ef.inc.TEh);end;
    if isempty(ef.amplitude(ef.dif.TMh,ef.inc.TMh))==0;ref_TM_TM(zou)=ef.amplitude(ef.dif.TMh,ef.inc.TMh);end;
    if isempty(ef.amplitude(ef.dif.TEh,ef.inc.TMh))==0;ref_TM_TE(zou)=ef.amplitude(ef.dif.TEh,ef.inc.TMh);end;
    R0_TE_TE(zou)=abs(ref_TE_TE(zou))^2;R0_TE_TM(zou)=abs(ref_TE_TM(zou))^2;R0_TM_TM(zou)=abs(ref_TM_TM(zou))^2;R0_TM_TE(zou)=abs(ref_TM_TE(zou))^2;
    if isempty(ef.amplitude(ef.dif.TEb,ef.inc.TEb))==0;refsub_TE_TE(zou)=ef.amplitude(ef.dif.TEb,ef.inc.TEb);end;
    if isempty(ef.amplitude(ef.dif.TMb,ef.inc.TEb))==0;refsub_TE_TM(zou)=ef.amplitude(ef.dif.TMb,ef.inc.TEb);end;
    if isempty(ef.amplitude(ef.dif.TMb,ef.inc.TMb))==0;refsub_TM_TM(zou)=ef.amplitude(ef.dif.TMb,ef.inc.TMb);end;
    if isempty(ef.amplitude(ef.dif.TEb,ef.inc.TMb))==0;refsub_TM_TE(zou)=ef.amplitude(ef.dif.TEb,ef.inc.TMb);end;
    R0sub_TE_TE(zou)=abs(refsub_TE_TE(zou))^2;R0sub_TE_TM(zou)=abs(refsub_TE_TM(zou))^2;R0sub_TM_TM(zou)=abs(refsub_TM_TM(zou))^2;R0sub_TM_TE(zou)=abs(refsub_TM_TE(zou))^2;

    if isempty(ef.amplitude(ef.dif.TEb,ef.inc.TEh))==0;trans_TE_TE(zou)=ef.amplitude(ef.dif.TEb,ef.inc.TEh);end;
    if isempty(ef.amplitude(ef.dif.TMb,ef.inc.TEh))==0;trans_TE_TM(zou)=ef.amplitude(ef.dif.TMb,ef.inc.TEh);end;
    if isempty(ef.amplitude(ef.dif.TMb,ef.inc.TMh))==0;trans_TM_TM(zou)=ef.amplitude(ef.dif.TMb,ef.inc.TMh);end;
    if isempty(ef.amplitude(ef.dif.TEb,ef.inc.TMh))==0;trans_TM_TE(zou)=ef.amplitude(ef.dif.TEb,ef.inc.TMh);end;
    T0_TE_TE(zou)=abs(trans_TE_TE(zou))^2;T0_TE_TM(zou)=abs(trans_TE_TM(zou))^2;T0_TM_TM(zou)=abs(trans_TM_TM(zou))^2;T0_TM_TE(zou)=abs(trans_TM_TE(zou))^2;

    if trace_champ==1
     
    [x,wx]=retgauss(-2*periodicity_x/2,2*periodicity_x/2,15,12,[-diameter_x/2,diameter_x/2]);
    [y,wy]=retgauss(-periodicity_y/2,periodicity_y/2,15,12,[-diameter_y/2,diameter_y/2]);    
        
    
    tab_norm=[0,1,1];
    if size(ef.inc.teta,2)==1
        inc=1;
    elseif size(ef.inc.teta,2)==2
        inc=[0,0];
        if pol==0;inc(ef.inc.TE)=1;elseif pol==2;inc(ef.inc.TM)=1;end;
    end
    if sens==1
        sb=rettronc(sb,[],0,-1);
        sb_norm=retb(init,ah,-0.1,0,[],[]);
        [ei,zi]=retchamp(init,{ah},sh,sb_norm,inc,{x,y},tab_norm,[],1:6,1,1,1:6);
        flux=retpoynting(ei,[0,0,-1],wx,wy,[]);
        inc=1/sqrt(flux).*inc;
        [einc,zinc]=retchamp(init,{ah},sh,sb_norm,inc,{0,0},tab_norm,[],1:6,1,1,1:6);
        Einc=squeeze(einc(:,:,:,1:3));
        Hinc=squeeze(einc(:,:,:,4:6));
    elseif sens==2
        sh=rettronc(sh,[],0,1);
        sh_norm=retb(init,ab,0.1,0,[],[]);
        [ei,zi]=retchamp(init,{ab},sh_norm,sb,inc,{x,y},tab_norm,[],1:6,1,1,1:6);
        flux=retpoynting(ei,[0,0,-1],wx,wy,[]);
        inc=1/sqrt(flux).*inc;
        [einc,zinc]=retchamp(init,{ab},sh_norm,sb,inc,{0,0},tab_norm,[],1:6,1,1,1:6);
        Einc=squeeze(einc(:,:,:,1:3));
        Hinc=squeeze(einc(:,:,:,4:6));
    end
        
    tab0=[[h_air,1,50];[hg,2,50];[hm,3,50]];
    struct0={ah,ag,am};
    for az=1:nrep
        tab0=[tab0;[h1,4,50];[h2,5,50]];
        struct0={struct0{:},a1,a2};
    end
    struct0={struct0{:},ab};
    tab0=[tab0;[h_sub,6,50]];
        
    if isempty(x0)==1&&isempty(z0)==1
        [e0,z,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x,y0},tab0,[],[1:6]+7.25*i,1,1,1:6);
    elseif isempty(y0)==1&&isempty(z0)==1
        [e0,z,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x0,y},tab0,[],[1:6]+7.25*i,1,1,1:6);
    elseif isempty(x0)==1&&isempty(y0)==1
        tab0(:,3)=0;
        HH=cumsum(tab0(:,1));
        numz=1;
        while (HH(end)-z0)>HH(numz);numz=numz+1;end;
        tab1=[tab0(1:numz-1,:);[tab0(numz,1)-(z0-sum(tab0(numz+1:end,1))),numz,0];[0,numz,1];[z0-sum(tab0(numz+1:end,1)),numz,0];tab0(numz+1:end,:)];
        [e0,z,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x,y},tab1,[],[1:6]+7.25*i,1,1,1:6);
    end
        
    for ii=1:3                                  
        o0(:,:,:,ii+3)=o0(:,:,:,ii+3)./o0(:,:,:,ii);
        o0(:,:,:,ii)=1;
    end
    indice=squeeze(sqrt(o0(:,:,:,4)));
    Ex=squeeze(e0(:,:,:,1));Ey=squeeze(e0(:,:,:,2));Ez=squeeze(e0(:,:,:,3));
    Hx=squeeze(e0(:,:,:,4));Hy=squeeze(e0(:,:,:,5));Hz=squeeze(e0(:,:,:,6));
    end
    
    %%%% Sauvegarde des data
    text='test';
    save(text,'wavelength','ref_TE_TE','ref_TE_TM','ref_TM_TM','ref_TM_TE','R0_TE_TE','R0_TE_TM','R0_TM_TM','R0_TM_TE','nh','ng','ngm','nm','n1','n2','nsub','periodicity_x','periodicity_y','diameter_x','diameter_y','hg','hm','h1','h2','Mx','My','pol',...
        'op_granet','x0','y0','z0','h_air','h_sub','Einc','Hinc','Ex','Ey','Ez','Hx','Hy','Hz','x','y','z','indice','trace_champ')

    retio
    toc
end

  Intensity = (abs(Ex).^2 + abs(Ez).^2) /(Einc' * Einc);

 %   Intensity=real(Ez);
%% figure props
width = 4;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize


%%%% Exemple de figure pour tracer le champ pour une calcul fait avec trace_champ=1, x0=[], y0=0, z0=[]
%%%% Hy en fonction de x et z
%%%% Les contours de l'objet sont en blancs (indice contient la carte des indices) 
if trace_champ==1
      save(sprintf('nf%s',int2str(sens)), 'x','z','Intensity' );
    figure
    pcolor(x,z,Intensity);shading flat;colorbar;colormap(hot);hold on;
    contour(x,z,indice,'w','linewidth',1.5)
else 
    
      save('ff', 'wavelength','R0_TM_TM','R0sub_TM_TM' );
    h=figure
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); % Set size
    set(gca, 'FontSize', fsz, 'LineWidth', alw); % Set properties
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

    plot(wavelength*1000, R0_TM_TM, wavelength*1000, R0sub_TM_TM,'LineWidth',lw,'MarkerSize',msz)
    xlim([min(wavelength) max(wavelength)]*1e3)
    ylim([0 1])
    xlabel('wavelength /nm')
    ylabel('R')
    title('Basic Tamm')
   % saveas(h,'apl','fig')
   % saveTightFigure(h,'apl.pdf')
end

