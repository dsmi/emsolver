% Inductance of a straight wire, Rosa formula
function L = rosaind( l, r )
   h = l;
   s = sqrt( r.*r + h.*h ); 
   L = mu0 ./ ( 2.0 * pi ) * ( h.*log( ( h + s ) ./ r ) + ( r - s ) );
end    
