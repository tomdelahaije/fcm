% This matlab script is used to generate the sharpening factors described
% in the paper:
%
%    A. Tristan-Vega, S. Aja-Fernandez, and C.-F. Westin
%    "Deblurring of Probabilistic ODFs in Quantitative Diffusion MRI"
%    In Proceedings of the IEEE Intl. Sym. on Biomed. Im., ISBI 2012
%    Barcelona (Spain), pp. 932-935 
%
% They are used to paliate the FRT-blurring in Q-Balls like estimators.
% Please, see also the file HARDI-ITK/GenerateDeconvolutionFactors.h

bDmin = 0.1;
bDmax = 100;
NS    = 40;
bD    = exp(  linspace( log(bDmin), log(bDmax), NS )  );

L     = 8;

NP    = 10000;
theta = pi*(0:NP-1)/NP;
LEG   = zeros( L/2+1, NP );
cmp   = 1;
for l=0:2:L
    dumb         = legendre( l, cos(theta) ) * sqrt( 4*pi./(2*l+1) ) * (2*pi)*(-1)^(l/2)*cmp;
    cmp          = cmp * (l+1)/(l+2);
    LEG(l/2+1,:) = (  (2*pi)*sqrt( (2*l+1)/4/pi )*(theta(2)-theta(1))  ) * dumb(1,:).*sin(theta); 
end

RESULTS = zeros( L/2+1, NS );

p1 = ( abs(cos(theta)) >= 1e-6 );
p2 = ( abs(cos(theta)) <  1e-6 );
CC = cos(theta).*cos(theta);
SS = 1-CC;

for n=1:NS
    
    F = 2*bD(n).*SS.*CC.*exp( -bD(n).*CC );
    F = F - ( 1 + CC ).*( 1 - exp(-bD(n).*CC) );
    
    F(p1) = F(p1)./( 8*pi*pi.*CC(p1) );
    F(p2) = bD(n)/( 8*pi*pi );
    
    F     = F + 1/8/pi/pi;
    
    F = ones( L/2+1, 1 ) * F;
    F = F .* LEG;
    
    
    RESULTS(:,n) = sum( F, 2 );
    %RESULTS(1,n) = 1;
    
end


COLORS = rand( size(RESULTS,1), 3 );

hf = figure(1);
hold on;
grid on;

hl  = zeros( 1, size(RESULTS,1) );
hlt = cell(  1, size(RESULTS,1) );
for k=1:size(RESULTS,1)
    hl(k)  = plot( bD, real(log(RESULTS(k,:)))/log(10), 'Color', COLORS(k,:), 'LineStyle', '-', 'LineWidth', 2 );
    hlt{k} = ['l = ',num2str(2*(k-1))];
end
hl = legend( hl, hlt{:}, 'Location', 'SouthEast' );
xlabel('b{\cdot}D','FontSize',20,'FontName','Times');
ylabel('log_{10}(f_{l})','FontSize',20,'FontName','Times');
set(hl,'FontSize',16,'FontName','Times');
hg = get(hf,'CurrentAxes');
set(hg,'FontSize',14);

hf = figure(2);
hold on;
grid on;

hl  = zeros( 1, size(RESULTS,1) );
hlt = cell(  1, size(RESULTS,1) );
for k=1:size(RESULTS,1)
    hl(k)  = plot( bD, RESULTS(k,:), 'Color', COLORS(k,:), 'LineStyle', '-', 'LineWidth', 2 );
    hlt{k} = ['l = ',num2str(2*(k-1))];
end
hl = legend( hl, hlt{:}, 'Location', 'NorthWest' );
xlabel('b{\cdot}D','FontSize',20,'FontName','Times');
ylabel('f_{l}','FontSize',20,'FontName','Times');
set(hl,'FontSize',16,'FontName','Times');
hg = get(hf,'CurrentAxes');
set(hg,'FontSize',14);

save(['sharpen_factors_LMAX',num2str(L),'.mat'],'RESULTS','bD');

for k=0:NS-1
    for l=1:L/2
        fprintf(1,'factors[%d][%d] = %1.8f; ',k,l,1./RESULTS(l+1,k+1));
    end
    fprintf(1,'\n');
end
fprintf(1,'\n');
for k=0:NS-1
    fprintf(1,'bD[%d]=%1.8f; ',k,bD(k+1));
    if( rem(k,5)==4 )
        fprintf(1,'\n');
    end
end
