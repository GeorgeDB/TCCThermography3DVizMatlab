function [self,u1]=refinePDEMmesh(self,u,it,mode)
%REFINEPDEMMESH Refine a triangular mesh.
%
%       [PDEM]=REFINEPDEMMESH(PDEM) returns a refined version
%       of the triangular mesh specified within pde model PDEM.
%
%       To refine the mesh, the triangular mesh data P, E, and T are
%       extracted from PDEM.
%
%       [PDEM,U1]=REFINEPDEMMESH(PDEM,U) refines the mesh and also
%       extends the function U to the new mesh by linear interpolation.
%       The number of rows in U should correspond to the number of
%       nodes in the PDEM mesh, and U1 will have as many rows as there are
%       points in the new refined mesh. 
%       Each column of U is interpolated separately.
%
%       An extra input argument is interpreted as a list of
%       subdomains to refine, if it is a row vector, or a list of
%       triangles to refine, if it is a column vector.
%
%       The default refinement method is regular refinement, where
%       all of the specified triangles are divided into four triangles
%       of the same shape. Longest edge refinement, where
%       the longest edge of each specified triangle is bisected, can
%       be demanded by giving 'longest' as a final parameter. Using
%       'regular' as the final parameter results in regular refinement.
%       Some triangles outside of the specified set may also be refined,
%       in order to preserve the triangulation and its quality.
%
%       See also INITMESH, PDEGEOM
%
%       Note: This function does not support quadratic triangular elements
%
%       Original function refinemesh() by
%       A. Nordmark 4-26-94, AN 6-21-94.
%       Copyright 1994-2017 The MathWorks, Inc.
%
%       Functionality ported to Matlab's new PDE toolbox/pde model
%       environment by
%       I. Jentz 10-2-2020

%%%%%%%%%%%%% Additions by I. Jentz, 10-2-20
% old output was [p1,e1,t1,u1]
g = self.Geometry;
[p,e,t] = self.Mesh.meshToPet;
%%%%%%%%%%%%% End of additions by I. Jentz, 10-2-20

% if nargin > 0
%     g = convertStringsToChars(g);
% end

if nargin > 1
    u = convertStringsToChars(u);
end

if nargin > 2
    it = convertStringsToChars(it);
end

if nargin > 3
    mode = convertStringsToChars(mode);
end

np=size(p,2);
ne=size(e,2);
nt=size(t,2);
% only linear elements are currently supported
%quadraticTriNotSupported(t);
if size(t, 1) == 7
        error(message('pde:pdeSanityChecker:QuadraticTriNotSupported'));
end

if nargout==2
  intp=1;
else
  intp=0;
end

if nargin-intp==1
  it=(1:nt)';                           % All triangles
  mode='regular';
end

if (~intp) && nargin==2
  it=u;
  if ischar(it)
    mode=it;
    it=(1:nt)';                         % All triangles
  else
    mode='regular';
  end
end

if (~intp) && nargin==3
  mode=it;
  it=u;
end

if intp && nargin==3
  if ischar(it)
    mode=it;
    it=(1:nt)';                         % All triangles
  else
    mode='regular';
  end
end

if strcmp(mode,'regular')==0 && strcmp(mode,'longest')==0
  error(message('pde:refinemesh:InvalidRefineMode'));
end

if size(it,1)>1                        % Triangles
  it=it';
else                                    % Subdomains
  it=pdesdt(t,it);
end

% Cannot use matrix indices that exceeds the size of a signed int
[comp,maxsize]=computer;
indexproblem=np^2>maxsize;

% Find longest side of each triangle
ls=3*ones(1,nt);
d1=(p(1,t(1,:))-p(1,t(2,:))).^2+(p(2,t(1,:))-p(2,t(2,:))).^2;
d=(p(1,t(2,:))-p(1,t(3,:))).^2+(p(2,t(2,:))-p(2,t(3,:))).^2;
ii=find(d>d1);
ls(ii)=1*ones(size(ii));
d1=max(d,d1);
d=(p(1,t(3,:))-p(1,t(1,:))).^2+(p(2,t(3,:))-p(2,t(1,:))).^2;
ii=find(d>d1);
ls(ii)=2*ones(size(ii));
% Permute so longest side is 3
ii=find(ls==1);
d=t(1,ii);
t(1,ii)=t(2,ii);
t(2,ii)=t(3,ii);
t(3,ii)=d;
ii=find(ls==2);
d=t(1,ii);
t(1,ii)=t(3,ii);
t(3,ii)=t(2,ii);
t(2,ii)=d;

itt1=ones(1,nt);
itt1(it)=zeros(size(it));
it1=find(itt1);                         % Triangles not yet to be refined
it=find(itt1==0);                       % Triangles whos longest side is to be bisected

% Make a connectivity matrix, with edges to be refined.
% -1 means no point is yet allocated
ip1=t(1,it);
ip2=t(2,it);
if strcmp(mode,'regular')
  ip3=t(3,it);
end
A=sparse(ip1,ip2,-1,np,np);
if strcmp(mode,'regular')
  A=A+sparse(ip2,ip3,-1,np,np);
  A=A+sparse(ip3,ip1,-1,np,np);
end
A=-((A+A.')<0);
newpoints=1;

% loop until no additional hanging nodes are introduced
while newpoints
  newpoints=0;
  n=length(it1);
  ip1=t(1,it1);
  ip2=t(2,it1);
  ip3=t(3,it1);
  m1=zeros(1,n);
  m2=m1;
  m3=m1;
  for i=1:n
    m3(i)=A(ip1(i),ip2(i));
    m1(i)=A(ip2(i),ip3(i));
    m2(i)=A(ip3(i),ip1(i));
  end
  ii=find(m3);
  if length(ii)>0
    itt1(it1(ii))=zeros(size(ii));
  end
  ii=find((m1 | m2) & (~m3));
  if length(ii)>0
    A=A+sparse(ip1(ii),ip2(ii),-1,np,np);
    A=-((A+A.')<0);
    newpoints=1;
    itt1(it1(ii))=zeros(size(ii));
  end
  it1=find(itt1);                       % Triangles not yet fully refined
  it=find(itt1==0);                     % Triangles fully refined
end

% Find edges to be refined
if ~indexproblem
  ie=full(A(e(1,:)+(e(2,:)-1)*np))==-1;
else
  ie=l_extract(A,e(1,:),e(2,:))==-1;
end

ie1=find(ie==0);                        % Edges not to be refined
ie=find(ie);                            % Edges to be refined

%%%%%%%%%%%%% Additions by I. Jentz, 10-2-20
% Get the edge "midpoint" coordinates
% this is editied to pass in the pde.AnalyticGeometry object g.geom 
[x,y]=pdeigeom(g.geom,e(5,ie),(e(3,ie)+e(4,ie))/2);
%%%%%%%%%%%%% End of additions by I. Jentz, 10-2-20

% Create new points
p1=[p [x;y]];
if intp
  u1=[u;(u(e(1,ie),:)+u(e(2,ie),:))/2];
end
ip=(np+1):(np+length(ie));
np1=np+length(ie);
% Create new edges
e1=[e(:,ie1) ...
        [e(1,ie);ip;e(3,ie);(e(3,ie)+e(4,ie))/2;e(5:7,ie)] ...
        [ip;e(2,ie);(e(3,ie)+e(4,ie))/2;e(4,ie);e(5:7,ie)]];
% Fill in the new points
if ~indexproblem
  A(e(1,ie)+np*(e(2,ie)-1))=ip;
  A(e(2,ie)+np*(e(1,ie)-1))=ip;
else
  A=l_assign(A,[e(1,ie) e(2,ie)],[e(2,ie) e(1,ie)],[ip ip]);
end

% Generate points on interior edges
[i1,i2]=find(A==-1 & A.'==-1);
i=find(i2>i1);
i1=i1(i);
i2=i2(i);
p1=[p1 ((p(1:2,i1)+p(1:2,i2))/2)];
if intp
  u1=[u1;(u(i1,:)+u(i2,:))/2];
end
ip=(np1+1):(np1+length(i));
np1=np1+length(i);
% Fill in the new points
if ~indexproblem
  A(i1+np*(i2-1))=ip;
  A(i2+np*(i1-1))=ip;
else
  A=l_assign(A,[i1 i2],[i2 i1],[ip ip]);
end

% Lastly form the triangles
ip1=t(1,it);
ip2=t(2,it);
ip3=t(3,it);
if ~indexproblem
  mp1=full(A(ip2+np*(ip3-1)));
  mp2=full(A(ip3+np*(ip1-1)));
  mp3=full(A(ip1+np*(ip2-1)));
else
  mp1=l_extract(A,ip2,ip3);
  mp2=l_extract(A,ip3,ip1);
  mp3=l_extract(A,ip1,ip2);
end

% Find out which sides are refined
bm=1*(mp1>0)+2*(mp2>0);
% The number of new triangles
nt1=length(it1)+length(it)+sum(mp1>0)+sum(mp2>0)+sum(mp3>0);
t1=zeros(4,nt1);
t1(:,1:length(it1))=t(:,it1);           % The unrefined triangles
nnt1=length(it1);
if isempty(bm)
	i = bm;
else
	i=find(bm==3);                          % All sides are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp2(i);mp1(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[mp1(i);mp2(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
if isempty(bm)
	i = bm;
else
	i=find(bm==2);                          % Sides 2 and 3 are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp2(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
if isempty(bm)
	i = bm;
else
	i=find(bm==1);                          % Sides 3 and 1 are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp3(i);mp1(i);t(4,it(i))];
nnt1=nnt1+length(i);
if isempty(bm)
	i = bm;
else
	i=find(bm==0);                          % Side 3 is refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);

%%%%%%%%%%%%% Additions by I. Jentz, 10-2-20
%%%%% see pde\+pde\@EquationModel\generateMeshFacettedAnalytic.m
Hmax = self.Mesh.MaxElementSize;
Hmin = self.Mesh.MinElementSize;
Hgrad = self.Mesh.MeshGradation;
geomOrder = self.Mesh.GeometricOrder;
verts = 1:self.Geometry.NumVertices;
vas = self.Mesh.findNodes('region','Vertex',verts);
% F = self.facetAnalyticGeometry(g, Hmax);
% [nodes, tri, cas, fas, eas, vas, Hmax, Hmin, Hgrad, esense] = genMeshFromFacetRep(F, Hmax, Hmin, Hgrad, geomOrder);
[nodes, ndparams, tri, cas, fas, eas, esense] = assocInputsfromPET(p1,e1,t1,t);

assoc = pde.FEMeshAssociation(tri, cas, fas, eas, vas, esense); 
[nodes, ndparams] = projectEdgeNodesToCurves(nodes, g, assoc);
msh = pde.FEMesh(nodes,tri,Hmax,Hmin, Hgrad, geomOrder,assoc, ndparams);

self.Mesh.delete();
self.Mesh = msh;
addlistener(msh,'ObjectBeingDestroyed',@pde.PDEModel.delistMesh); 


% Want to see what methods are available in the obfuscated pde.FEMesh
% class? just run this stuff.

% meta = ?pde.FEMesh;
% 
% for i = 1:size(meta.MethodList,1)
%     method = meta.MethodList(i).Name;
%     try
%         val = msh.(method);
%         disp([method,' returns ',num2str(size(val,1)),'x',num2str(size(val,2)),' ',class(val)]);
%     catch
%         disp(['error with ',method])
%     end
% end

%%%%%%%%%%%%% End of additions by I. Jentz, 10-2-20

function k=l_extract(A,i,j)

if numel(i)~=numel(j)
  error(message('pde:refinemesh:ijNumel'))
end

k=zeros(size(i));

for l=1:numel(i)
  k(l)=A(i(l),j(l));
end

end

function A=l_assign(A,i,j,k)

if numel(i)~=numel(j) || numel(i)~=numel(k) 
  error(message('pde:refinemesh:ijkNumel'))
end

for l=1:numel(i)
  A(i(l),j(l))=k(l);
end

end

%%%%%%%%%%%%% Additions by I. Jentz, 10-2-20
%%%%Function appear below have been added

% projectEdgeNodesToCurves - project edge nodes to curves and compute
%                            nodal parametric positions.
function [nodes, ndparams] = projectEdgeNodesToCurves(nodes, gm, ma)
    numEdges = gm.NumEdges();
    ndparams = cell(1,numEdges);
    geom = gm.geom;
    isfcn = false;
    if(isa(geom, 'function_handle'))
        isfcn = true;
    end
    NumRealArcFacets = 12;
    for edgeId = 1:numEdges
        en = ma.getNodes('Edge',edgeId);
        xy = nodes(:,en)';
        darcL = sqrt(sum((diff(xy,1,1).^2), 2)); %#ok<UDIM> % Facet, lengths
        s = [0; cumsum(darcL)]; % Arc length at each node in polycurve
        s0 = 0; % Start parameter of polycurve
        s1 = s(end); % Total length of edge == end parameter of polycurve
        D=pdeigeom(geom,edgeId);
        sp = D(1);
        ep = D(2);        
        P = linspace(sp, ep, NumRealArcFacets);
        P(1) = sp;
        P(end) = ep;
        [X,Y]=pdeigeom(geom,edgeId,P);
        XY = [X; Y];
        pp = pdearcl(P,XY,s,s0,s1);
        % start and end parameters may be off by eps
        % Adjust to avoid NaNs
        pp(1) = sp;
        pp(end) = ep;
        ndparams{edgeId} = [pp(1:end-1); pp(2:end)];
        if isfcn || (geom(1,edgeId) ~= 2)
            [Xn,Yn]=pdeigeom(geom,edgeId,pp);        
            nodes(1,en) = Xn;
            nodes(2,en) = Yn;                
        end
    end
end


    function [nodes, ndparams, tri, cas, fas, eas, esense] = assocInputsfromPET(p1,e1,t1,t)
        
        % nodes are the simplest,
        nodes = p1;
        
        % tri is t1 re-organized so that the regions are in increasing
        % order, i.e. t1(4,:) is in ascending order
        ncell = max(t1(4,:));
        ntb = size(t,2);
        nt1 = size(t1,2);
        its = 1:ntb;        ts = t1(:,its);
        ite = (ntb+1):nt1;  te = t1(:,ite);
        itri = its;         tsort = ts;  
        itri = [itri(:,(tsort(4,:)<=1)),ite(:,(te(4,:)==1)),itri(:,(tsort(4,:)>1))];
        tsort = [tsort(:,(tsort(4,:)<=1)),te(:,(te(4,:)==1)),tsort(:,(tsort(4,:)>1))];
        if ncell > 1
        for ic = 2:ncell
            itri = [itri(:,(tsort(4,:)<=ic)),ite(:,(te(4,:)==ic)),itri(:,(tsort(4,:)>ic))];
            tsort = [tsort(:,(tsort(4,:)<=ic)),te(:,(te(4,:)==ic)),tsort(:,(tsort(4,:)>ic))];
        end
        end
        tri = tsort(1:3,:); % we only want to output points, can remove region info now

        
        % cas is cell association, which is the first and last triangle in
        % a region
        cas = zeros(2,ncell);
        for ic = 1:ncell
           cas(1,ic) = find(tsort(4,:)==ic,1,'first');
           cas(2,ic) = find(tsort(4,:)==ic,1,'last');
        end
        
        % fas is face association, which is the same as cas for 2D problems
        fas = cas;
               
        
        % esense is the direction relative to e1, true is in the same
        % direction
        esense = true(ne,1); % the way we calculated eas makes the direction the same
        % eas is edge association, which is a cell of edges with each cell
        % containing the triangles of the edge and the non-edge node of
        % each traignle
        ne = max(e1(5,:));
        eas = cell(ne,1);
        % ndparams is edge length parameters, which is a cell of edges with
        % each cell containing the proportion of length, from 0 to 1, of
        % the edge segment
        ndparams = cell(1,ne);
        % loop over all the edges
        for ie = 1:ne   
            
            % need to sort the edges so that they line up
            iedg = (e1(5,:)==ie);
            etmp = e1(:,iedg); 
            test = etmp(1,2:end)-etmp(2,1:(end-1));
            if any(test ~= 0)
                js = find(~ismember(etmp(1,:),etmp(2,:)));
            if js ~= 1
                ecnxt = etmp(:,1);
                etmp(:,1) = etmp(:,js);
                etmp(:,js) = ecnxt;           
            end
            for j = 1:size(etmp,2)-1
              js = find(etmp(1,:) == etmp(2,j));
              ecnxt = etmp(:,j+1);
              try
              etmp(:,j+1) = etmp(:,js);
              catch
                  keyboard
              end
              etmp(:,js) = ecnxt;
            end
            end
            
            if etmp(6,1) == 0
               etmp = flip(etmp,2); 
               etmp([3,4],:) = 1-flip(etmp([3,4],:),1);
               esense(ie) = false;
            end
            
            % set the eas{} of the edge
            ic = max(etmp(6:7,1)); % pick a region that neighbors the edge
            itc = (tsort(4,:) == ic); % triangle within this region
            eas{ie} = zeros(2,size(etmp,2));
            for j = 1:size(etmp,2)
                tt = find( (sum((tri == etmp(1,j)) + (tri == etmp(2,j))) == 2) & itc ); % find triangle that is on the edge
                eas{ie}(1,j) = tt;
                nen = tri(~( (tri(:,tt) == etmp(1,j)) | (tri(:,tt) == etmp(2,j)) ),tt); % triangle node that is not on the edge
                eas{ie}(2,j) = find(tri(:,tt)==nen);
            end
            
            % set the npparams{} of the edge
            ndparams{ie} = etmp(3:4,:);
        end 
        
    end
%%%%%%%%%%%%% End of additions by I. Jentz, 10-2-20

end

