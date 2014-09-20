% Function for spectral estimation of material refr. index dispersion


function [n,varargout]=RefractiveIndexDispersion(w,unit,material,varaging)

% Sellmeier equation {http://en.wikipedia.org/wiki/Sellmeier_equation}

%INPUT
%w: (Px1) P different wavelength [unit]
%unit: for wavelength
%material: material type
%varargin

%OUTPUT
%n: (Px1) refractive index at the P different wavelength [adim]
%varargout;  {1} abbe number

%NOTE: Coeff.s have to be consistent with the wavelength ( typically in um)

switch unit    %Express the wavelength in  micrometer [um]
    case {'mm'}
        w=w.*1e3;
    case {'m'}
        w=w.*1e6;
    otherwise
        warning('Check if wavelength is expressed in [mm]/[m]')
        material='none';
end

switch material
    
    case {'BASF6';'BaSF6'} %http://refractiveindex.info/legacy/?group=HIKARI&material=E-BASF6
        %CUSTOMIZED
        B=[2.71377616,-9.41072521*1e-3,2.21914885*1e-2,7.35180213*1e-4,-2.11522864*1e-5,3.73099495*1e-6];
        n2=B(1)+B(2).*w.^2+B(3).*w.^(-2)+B(4).*w.^(-4)+B(5).*w.^(-6)+B(6).*w.^(-8);
        [n]=sqrt(n2);  
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=41.96; %abbe number
        varargout{1}=Vd;
        
    case {'BASF7';'BaSF7'} %http://refractiveindex.info/legacy/?group=HIKARI&material=E-BASF7
        %CUSTOMIZED
        B=[2.83515039,-1.78508143*1e-3,1.53978598*1e-2,3.71905515*1e-3,-4.88482182*1e-4,3.29017911*1e-5];
        n2=B(1)+B(2).*w.^2+B(3).*w.^(-2)+B(4).*w.^(-4)+B(5).*w.^(-6)+B(6).*w.^(-8);
        [n]=sqrt(n2);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=41.17; %abbe number
        varargout{1}=Vd;
    
    case {'BK7'}  %http://refractiveindex.info/legacy/?group=GLASSES&material=BK7
         n_ord=3;
        B=[1.03961212,0.231792344,1.01046945];
        C=[0.00600069867,0.0200179144,103.560653];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=64.17; %abbe number
        varargout{1}=Vd;
        
    case {'BK10'}  %http://refractiveindex.info/legacy/?group=LZOS&material=BK10
         n_ord=3;
        B=[0.888308131,0.328964475,0.984610769];
        C=[0.00516900822,0.0161190045,99.7575331];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=56.05; %abbe number
        varargout{1}=Vd;
        
        
    case{'F4','S-F4'}  %http://refractiveindex.info/legacy/?group=CDGM&material=F4
         %CUSTOMIZED
        B=[2.5564435,-9.5602169*1e-3,2.195072*1e-2,1.0014378*1e-4,-3.9758963*1e-5,5.4421828*1e-6];
        n2=B(1)+B(2).*w.^2+B(3).*w.^(-2)+B(4).*w.^(-4)+B(5).*w.^(-6)+B(6).*w.^(-8);
        [n]=sqrt(n2);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=36.35; %abbe number
        varargout{1}=Vd;
        dnl=-0.0865; % [um-1]refractive index dispersione
        varargout{2}=dnl;
        
    case{'FK5'}  %http://refractiveindex.info/legacy/?group=SCHOTT&material=N-FK5
         n_ord=3;
        B=[0.844309338,0.344147824,0.910790213];
        C=[0.00475111955,0.0149814849,97.86002930];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=70.40; %abbe number
        varargout{1}=Vd;
        
    case{'K3'}  %http://refractiveindex.info/legacy/?group=HIKARI&material=E-K3
        %CUSTOMIZED
        B=[2.27078330,-7.72009861*1e-3,1.26286148*1e-2,-4.60005382*1e-5,3.30622876*1e-5,-1.37462973*1e-6];
        n2=B(1)+B(2).*w.^2+B(3).*w.^(-2)+B(4).*w.^(-4)+B(5).*w.^(-6)+B(6).*w.^(-8);
        [n]=sqrt(n2);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=58.93; %abbe number
        varargout{1}=Vd;
    
    case {'LAF2';'LaF2'}  %http://refractiveindex.info/legacy/?group=SCHOTT&material=N-LAF2
         n_ord=3;
        B=[1.80984227,0.15729555,1.0930037];
        C=[0.0101711622,0.0442431765,100.687748];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=44.85; %abbe number
        varargout{1}=Vd;
    
    case {'LaK8';'LAK8'}  %http://refractiveindex.info/legacy/?group=SCHOTT&material=N-LAK8
         n_ord=3;
        B=[1.33183167,0.546623206,1.19084015];
        C=[0.00620023871,0.0216465439,82.5827736];
        [n]=sellmeier_func(w,n_ord,B,C);      
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=53.83; %abbe number
        varargout{1}=Vd;
        
    case {'LaK9';'LAK9'}  %http://refractiveindex.info/legacy/?group=SCHOTT&material=N-LAK9
         n_ord=3;
        B=[1.46231905,0.344399589,1.15508372];
        C=[0.00724270156,0.0243353131,85.4686868];
        [n]=sellmeier_func(w,n_ord,B,C); 
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=54.71; %abbe number
        varargout{1}=Vd;
        
        
        
     case {'SF2'}  %http://refractiveindex.info/?shelf=glass&book=SCHOTT-SF&page=SF2
         n_ord=3;
        B=[1.40301821,0.231767504,0.939056586];
        C=[0.0105795466,0.0493226978,112.405955];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=33.85; %abbe number
        varargout{1}=Vd;
    
     case {'SF5'}  %http://refractiveindex.info/legacy/?group=GLASSES&material=SF5
         n_ord=3;
        B=[1.52481889,0.187085527,1.42729015];
        C=[0.011254756,0.0588995392,129.141675];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=32.25; %abbe number
        varargout{1}=Vd;
        
        
     case {'SF19';'S-F19';'F19'}  %http://refractiveindex.info/?shelf=glass&book=SCHOTT-inquiry&page=N-SF19
         n_ord=3;
        B=[1.52005444,0.17573947,1.43623424];
        C=[0.01096144,0.0593248486,126.795151];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=33.12; %abbe number
        varargout{1}=Vd;        
        dnl=-0.10152; % [um-1]refractive index dispersione
        varargout{2}=dnl;
        
    case{'SF56'}  %http://refractiveindex.info/?shelf=glass&book=SCHOTT-inquiry&page=N-SF56
         n_ord=3;
        B=[1.73562085,0.317487012,1.95398203];
        C=[0.0129624742,0.0612884288,161.559441];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=26.09; %abbe number
        varargout{1}=Vd;
        
     case {'SF57'}  %http://refractiveindex.info/legacy/?group=SCHOTT&material=N-SF57
         n_ord=3;
        B=[1.87543831,0.37375749,2.30001797];
        C=[0.0141749518,0.0640509927,177.389795];
        [n]=sellmeier_func(w,n_ord,B,C);
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=23.78; %abbe number
        varargout{1}=Vd;
        
    case {'SK4'}  %http://refractiveindex.info/?shelf=glass&book=SCHOTT-SK&page=N-SK4
         n_ord=3;
        B=[1.32993741,0.228542996,0.988465211];
        C=[0.00716874107,0.0246455892,100.886364];
        [n]=sellmeier_func(w,n_ord,B,C);   
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=58.63; %abbe number
        varargout{1}=Vd;
    
    
    case {'SK11'} %http://refractiveindex.info/?shelf=glass&book=SCHOTT-SK&page=N-SK11
        n_ord=3;
        B=[1.17963631,0.229817295,0.935789652];
        C=[0.00680282081,0.0219737205,101.513232];
        [n]=sellmeier_func(w,n_ord,B,C);   
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=60.80; %abbe number
        varargout{1}=Vd;
        dnl=-0.0480; % [um-1]refractive index dispersione
        varargout{2}=dnl;
        
    case {'SK16'} %http://refractiveindex.info/legacy/?group=SCHOTT&material=N-SK16
        n_ord=3;
        B=[1.34317774,0.241144399,0.994317969];
        C=[0.00704687339,0.0229005,92.7508526];
        [n]=sellmeier_func(w,n_ord,B,C);   
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=60.20; %abbe number
        varargout{1}=Vd;
        dnl=-0.0532; %[um-1]refractive index dispersione
        varargout{2}=dnl;
        

        
    case {'air'}
        [n]=airindex(w,'um');
        % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=89.3; %abbe number
        varargout{1}=Vd;
        
    case {'air0'}
    [n]=repmat([1],1,size(w,1));
    % Abbe number (Vd = (nd-1)/(nF-nC);d: 0.58756 µm (Yellow helium line)
        % F: 0.48613 µm (Blue hydrogen line)C: 0.65627 µm (Red hydrogen line))
        Vd=89.3; %abbe number
        varargout{1}=Vd;
        
    otherwise
        error(['Refractive index of ',material,'is not AVAILABLE :-( !'])
        n=[];

        
end

n=n'; %column vector