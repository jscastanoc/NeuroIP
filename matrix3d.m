classdef matrix3d < double
    
    properties
        d3 = 1;
    end
    methods
        function obj = matrix3d(data,d3)
            obj = obj@double(data);
            obj.d3 = d3;
        end
        function obj = subsref(A,S)
            SS = S;
            SS.subs ={S.subs{1},S.subs{2},1};
            obj = builtin('subsref',double(A),SS);
        end
        function varargout = size(obj,varargin)
            if numel(varargin)
                siz = builtin('size',double(obj),varargin{:});
                if varargin{1} == 3
                    siz = obj.d3;
                end
            else
                siz = builtin('size',double(obj));
                siz(3) = obj.d3;
            end
            if (nargout == 1)
                varargout{1} = siz;
            else
                for i = 1:nargout
                    varargout{i} = siz(i);
                end
            end
        end
    end
end