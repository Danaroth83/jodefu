function varargout = isString(s,varargin)
% ISSTRING Require a char row vector, or ''
%   T = ISSTRING(S) returns true if S is a 1-by-M character vector, where M
%   can be zero, or if S is the empty string '', and returns false otherwise.
%   
%   T = ISSTRING(S,ALLOWCELL), when ALLOWCELL is true, returns true if S is
%   either a string or a scalar cell array containing a single string, and
%   returns false otherwise.
%   
%   When T is true, [T,ASCHAR] = ISSTRING(S,ALLOWCELL) also returns either S,
%   or S{1} if S is a scalar cell containing a single string.  When T is
%   false, ISSTRING returns '' in ASCELLSTR.
%   
%   See also STRINGS.


%   Copyright 2011-2015 The MathWorks, Inc.

[varargout{1:nargout}] = sli_isString(s,varargin{:});