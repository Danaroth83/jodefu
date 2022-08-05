function checkSupportedNumeric(name,x,okInteger,okSparse,okComplex,okGpu)
%   checkSupportedNumeric(NAME,VAL) checks that the input named NAME has a
%   value VAL of a numeric type supported by code in the
%   Statistics and Machine Learning Toolbox
%
%   checkSupportedNumeric(NAME,VAL,TRUE) marks integer values as okay.
%
%   checkSupportedNumeric(NAME,VAL,OKINT,TRUE) marks sparse values as okay.
%
%   checkSupportedNumeric(NAME,VAL,OKINT,OKSPRS,TRUE) marks complex values as okay.

%   Copyright 2014-2018 The MathWorks, Inc.

m = [];

if ~isnumeric(x)
      m = message('stats:internal:utils:NumericRequiredNamed',name);
elseif isobject(x) && ~isa(x, 'gpuArray')
    m = message('stats:internal:utils:NoObjectsNamed',name);
elseif (nargin<6 || ~okGpu) && isa(x, 'gpuArray')
    m = message('stats:internal:utils:NoObjectsNamed',name);
elseif (nargin<5 || ~okComplex) && ~isreal(x)
    m = message('stats:internal:utils:NoComplexNamed',name);
elseif (nargin<4 || ~okSparse) && issparse(x)
    m = message('stats:internal:utils:NoSparseNamed',name);
elseif (nargin<3 || ~okInteger) && ~isfloat(x)
    m = message('stats:internal:utils:FloatRequiredNamed',name);
end
if ~isempty(m)
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
end


