function SXPWrite(freq, s, FileName, unit_adj, comment)
% SXPWrite(freq, S, FileName, unit_adj, COMMENT)
%
% writes multiport parameter data S to an .sxp file data 
% using the MDIF format (a.k.a. HPEEsof format); for a detailed
% description of data format see SXPParse.m
%
% freq is multiplied by unit_adj and written as MHz,
% so make sure you pass the right number (default 1e-6, Hz -> MHz)
% some very optional params still have to be edited by hand
% COMMENT is a string to be added in the file header for future reference
%
% See also SXPParse.
%
% written by tudor dima, tudima@zahoo.com, change the z into y

% ver 1.41: 2009.09.05  - uSXPstrfit as subfunction, some cleanup

if nargin < 5
   comment = 'MDIF file - unknown source';
end;

if nargin < 4
   unit_adj = 1e-6;
end;

N = max(size(freq));
order = size(s,1);
if order~=1 && size(s,2) ~= order
   disp('data does not seem to be valid');
end;

s_digits = 6; r_digits = 9; % f_digits = 6;
fprintf(1,'\n%s', ['writing parameter data to file ' FileName '...']);

% --- start writing the data ---
fid = fopen(FileName, 'wt');

fprintf(fid, '%s\n', ['! ' comment]);
fprintf(fid, '%s\n', '# MHz S RI R 50');

if order > 2
   for i = 1:N
      candidate = num2str(freq(i)*unit_adj);
      word = uSXPstrfit(candidate, r_digits);
      for j = 1:order
         candidate = num2str(real(s(1,j,i)), s_digits);
         word = [word uSXPstrfit(candidate, r_digits)];
         candidate = num2str(imag(s(1,j,i)), s_digits);
         word = [word uSXPstrfit(candidate, r_digits)];
      end;
      fprintf(fid, '%s\n', word);
      
      for k = 2:order
         word = uSXPstrfit(' ', r_digits);
         for j = 1:order
            candidate = num2str(real(s(k,j,i)), s_digits);
            word = [word uSXPstrfit(candidate, r_digits)];
            candidate = num2str(imag(s(k,j,i)), s_digits);
            word = [word uSXPstrfit(candidate, r_digits)];     
         end;
         fprintf(fid, '%s\n', word);
      end;
   end;
elseif order == 2 
   for i = 1:N
      word = [uSXPstrfit(num2str(freq(i)*unit_adj), r_digits) ...
          uSXPstrfit(num2str(real(s(1,1,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(imag(s(1,1,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(real(s(2,1,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(imag(s(2,1,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(real(s(1,2,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(imag(s(1,2,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(real(s(2,2,i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(imag(s(2,2,i)), s_digits), r_digits)];
      fprintf(fid, '%s\n', word);
   end;
elseif order == 1
   for i = 1:N
      word = [uSXPstrfit(num2str(freq(i)*unit_adj), r_digits) ...
          uSXPstrfit(num2str(real(s(i)), s_digits), r_digits) ...
          uSXPstrfit(num2str(imag(s(i)), s_digits), r_digits)];
      fprintf(fid, '%s\n', word);
   end;
end;

fclose(fid);
fprintf(1,'\n%s', '...done');

end

function word = uSXPstrfit(candidate, nr_of_chars)

% fits a string to a set length
% (useful for writing matrices in ascii)
%
% older uSXPstrfit (11 may 1999)

lC = size(candidate,2);

if lC >= nr_of_chars
    word = [candidate(1:nr_of_chars-1) ' '];
else
    word(1:nr_of_chars) = ' ';
    word(1:lC) = candidate;
end

end