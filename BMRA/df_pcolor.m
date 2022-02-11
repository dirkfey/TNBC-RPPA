function varargout = df_pcolor(varargin)

%[nx, ny] = size(data);
if nargin==1
    data = varargin{1};
elseif nargin==3
    X = varargin{1};
    Y = varargin{2};
    data = varargin{3}
else
    error('Wrong number of inputs')
end

data = [data data(:,end)];
data = [data; data(end,:)];

if nargin==1
    h = pcolor(data);
else
    h = pcolor(X,Y,data);
end

if nargout>0
    varargout{1} = h;
end

end