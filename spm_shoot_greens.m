function varargout = spm_shoot_greens(varargin)
% Build and apply FFT of Green's function (to map from momentum to velocity)
% FORMAT v = spm_shoot_greens(m,K,prm,bnd)
% m    - Momentum field n1*n2*n3*3 (single prec. float)
% K    - Fourier transform representation of Green's function
%        - either size n1*n2*n3 or n1*n2*n3*3*3
% prm  - Differential operator parameters (3 voxel sizes, 5 hyper-parameters)
%        - only needed when K is of size n1*n2*n3, in which case, voxel sizes
%          are necessary for dealing with each component individually
% bnd  - Boundary type:
%           0 = circulant [default]
%           1 = neumann
%           2 = dirichlet
%           3 = sliding
% v    - velocity field
%
% FORMAT [K,ld] = spm_shoot_greens('kernel',dm,prm,bnd)
% dm  - dimensions n1*n2*n3
% prm - Differential operator parameters (3 voxel sizes, 5 hyper-parameters)
% bnd - Boundary type [0]
% K   - Fourier transform representation of Green's function
%        - either size n1*n2*n3 or n1*n2*n3*3*3
% ld(1)  - Log determinant of operator
% ld(2)  - Number of degrees of freedom
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2012)

% Reference for symmetric convolution (i.e. which DCT/DST order to use):
% https://en.wikipedia.org/wiki/Symmetric_convolution

% John Ashburner
% $Id: spm_shoot_greens.m 7054 2017-04-04 12:09:32Z john $

if isa(varargin{1},'char') && strcmp(varargin{1},'kernel'),
    d   = varargin{2};
    prm = varargin{3};
    if nargin < 4
        bnd = 0;
    else
        bnd = varargin{4};
    end
    spm_diffeo('boundary',bnd);
    
    switch bnd
        case 0   % Circulant
            % The differential operator is symmetric, so the Fourier 
            % transform should be real
            dt = @(a, dim) real(fft(a, [], dim));
        case {1,2,3}   % Neumann/Dirichlet/Sliding
            dt = @(a, dim) dctr(a, dim, 'Type', 2);
        otherwise
            error('Boundary type %d does not exist. Should be in 0..3', bnd);
    end

    F = spm_diffeo('kernel',d,prm);
    if size(F,4) == 1
        % Diagonal and off-diagonal conditions are the same
        F = dt(dt(dt(F,1),2),3);
        sm = numel(F);
        if nargout >=2
            ld = log(F);
            if prm(4)==0, ld(1,1,1) = 0; end
            ld = -sum(ld(:));
        end
        if prm(4)==0
            F(1,1,1) = 0;
            sm       = sm - 1;
        end
        if nargout >=2
           ld = 3*ld + sm*sum(2*log(prm(1:3)));
        end
    else
        for i=1:size(F,4),
            F(:,:,:,i,1) = dt(dt(dt(F(:,:,:,i,1),1),2),3);
            F(:,:,:,i,2) = dt(dt(dt(F(:,:,:,i,2),1),2),3);
            F(:,:,:,i,3) = dt(dt(dt(F(:,:,:,i,3),1),2),3);
        end
        ld = 0;
        sm = 0;
        for k=1:size(F,3),
            % Compare the following with inverting a 3x3 matrix...
            A   = F(:,:,k,:,:);
            dt  = A(:,:,:,1,1).*(A(:,:,:,2,2).*A(:,:,:,3,3) - A(:,:,:,2,3).*A(:,:,:,3,2)) +...
                  A(:,:,:,1,2).*(A(:,:,:,2,3).*A(:,:,:,3,1) - A(:,:,:,2,1).*A(:,:,:,3,3)) +...
                  A(:,:,:,1,3).*(A(:,:,:,2,1).*A(:,:,:,3,2) - A(:,:,:,2,2).*A(:,:,:,3,1));
            msk     = dt<=0;
            if prm(4)==0 && k==1, msk(1,1,1) = true; end
            dt      = 1./dt;
            dt(msk) = 0;
            if nargout>=2
                sm      = sm + sum(sum(~msk));
                ld      = ld - sum(log(dt(~msk)));
            end
            F(:,:,k,1,1) = (A(:,:,:,2,2).*A(:,:,:,3,3) - A(:,:,:,2,3).*A(:,:,:,3,2)).*dt;
            F(:,:,k,2,1) = (A(:,:,:,2,3).*A(:,:,:,3,1) - A(:,:,:,2,1).*A(:,:,:,3,3)).*dt;
            F(:,:,k,3,1) = (A(:,:,:,2,1).*A(:,:,:,3,2) - A(:,:,:,2,2).*A(:,:,:,3,1)).*dt;

            F(:,:,k,1,2) = (A(:,:,:,1,3).*A(:,:,:,3,2) - A(:,:,:,1,2).*A(:,:,:,3,3)).*dt;
            F(:,:,k,2,2) = (A(:,:,:,1,1).*A(:,:,:,3,3) - A(:,:,:,1,3).*A(:,:,:,3,1)).*dt;
            F(:,:,k,3,2) = (A(:,:,:,1,2).*A(:,:,:,3,1) - A(:,:,:,1,1).*A(:,:,:,3,2)).*dt;

            F(:,:,k,1,3) = (A(:,:,:,1,2).*A(:,:,:,2,3) - A(:,:,:,1,3).*A(:,:,:,2,2)).*dt;
            F(:,:,k,2,3) = (A(:,:,:,1,3).*A(:,:,:,2,1) - A(:,:,:,1,1).*A(:,:,:,2,3)).*dt;
            F(:,:,k,3,3) = (A(:,:,:,1,1).*A(:,:,:,2,2) - A(:,:,:,1,2).*A(:,:,:,2,1)).*dt;
        end
    end
    varargout{1} = F;
    if nargout>=2
        varargout{2} = [ld, sm];
    end

else
    % Convolve with the Green's function via Fourier methods
    m = varargin{1};
    F = varargin{2};
    if nargin < 4
        bnd = 0;
    else
        bnd = varargin{4};
    end
    switch bnd
        case 0   % Circulant
            dtd = @(a, dim) fft(a, [], dim);
            dto = @(a, dim) fft(a, [], dim);
            itd = @(a, dim) ifft(a, [], dim, 'symmetric');
            ito = @(a, dim) ifft(a, [], dim, 'symmetric');
        case 1   % Neumann
            dtd = @(a, dim) dctr(a, dim, 'Type', 2);
            dto = @(a, dim) dctr(a, dim, 'Type', 2);
            itd = @(a, dim) idctr(a, dim, 'Type', 1);
            ito = @(a, dim) idctr(a, dim, 'Type', 1);
        case 2   % Dirichlet
            dtd = @(a, dim) dstr(a, dim, 'Type', 2);
            dto = @(a, dim) dstr(a, dim, 'Type', 2);
            itd = @(a, dim) idstr(a, dim, 'Type', 1);
            ito = @(a, dim) idstr(a, dim, 'Type', 1);
        case 3   % Sliding
            dtd = @(a, dim) dstr(a, dim, 'Type', 2);
            dto = @(a, dim) dctr(a, dim, 'Type', 2);
            itd = @(a, dim) idstr(a, dim, 'Type', 1);
            ito = @(a, dim) idctr(a, dim, 'Type', 1);
        otherwise
            error('Boundary type %d does not exist. Should be in 0..3', bnd);
    end
    
    v = zeros(size(m),'single');
    if size(F,4) == 1,
        % Simple case where convolution is done one field at a time
        prm = varargin{3};
        v(:,:,:,1) = ito(ito(itd(F.*dto(dto(dtd(m(:,:,:,1),1),2),3)*prm(1)^2,1),2),3);
        v(:,:,:,2) = ito(itd(ito(F.*dto(dtd(dto(m(:,:,:,2),1),2),3)*prm(2)^2,1),2),3);
        v(:,:,:,3) = itd(ito(ito(F.*dtd(dto(dto(m(:,:,:,3),1),2),3)*prm(3)^2,1),2),3);
    else
        % More complicated case for dealing with linear elasticity, where
        % convolution is not done one field at a time
        m(:,:,:,1) = dto(dto(dtd(m(:,:,:,1),1),2),3);
        m(:,:,:,2) = dto(dtd(dto(m(:,:,:,2),1),2),3);
        m(:,:,:,3) = dtd(dto(dto(m(:,:,:,3),1),2),3);
        for k=1:size(m,3),
            a = m(:,:,k,:);
            m(:,:,k,:) = 0;
            for j=1:3,
                for i=1:3,
                    m(:,:,k,j) = m(:,:,k,j) + F(:,:,k,j,i).*a(:,:,:,i);
                end
            end
        end
        v(:,:,:,1) = ito(ito(itd(m(:,:,:,1),1),2),3);
        v(:,:,:,2) = ito(itd(ito(m(:,:,:,2),1),2),3);
        v(:,:,:,3) = itd(ito(ito(m(:,:,:,3),1),2),3);
    end
    varargout{1} = v;
end

%__________________________________________________________________________________

% Robust DCT/DST that work with slices

function a = dctr(a, dim, varargin)
    if size(a, dim) > 1
        a = dct(a, [], dim, varargin{:});
    end


function a = dstr(a, dim, varargin)
    if size(a, dim) > 1
        a = dst(a, [], dim, varargin{:});
    end


function a = idctr(a, dim, varargin)
    if size(a, dim) > 1
        a = idct(a, [], dim, varargin{:});
    end


function a = idstr(a, dim, varargin)
    if size(a, dim) > 1
        a = idst(a, [], dim, varargin{:});
    end
