function [N,A] = patchnormals(FV)
%Vertex normals of a triangulated mesh, area weighted, left-hand-rule 
% N = patchnormals(FV) -struct with fields, faces Nx3 and vertices Mx3 
%N: vertex normals as Mx3
%A: vertex areas as Mx1
%face corners index 
A = FV.faces(:,1); 
B = FV.faces(:,2); 
C = FV.faces(:,3);
%face normals 
n = cross(FV.vertices(A,:)-FV.vertices(B,:),FV.vertices(C,:)-FV.vertices(A,:)); %area weighted
A = dot(n,n,2);
N = n./A;
%vertice normals 
% N = zeros(size(FV.faces,1),3); %init vertix normals 
% for i = 1:size(FV.faces,1) %step through faces (a vertex can be reference any number of times) 
% N(A(i),:) = N(i,:)+n(i,:); %sum face normals 
% N(B(i),:) = N(i,:)+n(i,:); 
% N(C(i),:) = N(i,:)+n(i,:); 
% end[