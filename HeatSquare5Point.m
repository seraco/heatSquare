classdef HeatSquare5Point < HeatSquareBase
    % Heat conduction for 5 point discretisation
    % Author: Sebastia Ramon
    
    methods
        
        function A = assemble_A(obj)
            A = sparse(zeros(obj.Nx * obj.Ny, obj.Nx * obj.Ny));
            a = 1 / obj.Delta_x / obj.Delta_x;
            b = 1 / obj.Delta_y / obj.Delta_y;
            c = -2 * (a + b);
            for i=1:obj.Ny
                for j=1:obj.Nx
                    inode = obj.Nx * (i-1) + j;
                    if i == 1
                        A(inode, inode) = 1;
                    elseif i == obj.Ny
                        A(inode, inode) = 1;
                    elseif j == 1
                        A(inode, inode) = 1;
                    elseif j == obj.Nx
                        A(inode, inode) = 1;
                    else
                        A(inode, inode - obj.Nx) = b;
                        A(inode, inode - 1) = a;
                        A(inode, inode) = c;
                        A(inode, inode + 1) = a;
                        A(inode, inode + obj.Nx) = b;
                    end
                end
            end
        end
        
        function b = assemble_b(obj)
            b = zeros(obj.Nx * obj.Ny, 1);
            for i=1:obj.Ny
                for j=1:obj.Nx
                    inode = obj.Nx * (i-1) + j;
                    if i == 1
                        x = (j - 1) * obj.Delta_x;
                        b(inode) = 1 + x;
                    elseif i == obj.Ny
                        b(inode) = 1;
                    elseif j == 1
                        b(inode) = 1;
                    elseif j == obj.Nx
                        y = (i - 1) * obj.Delta_y;
                        b(inode) = cos(6 * 3/2 * pi * y) + 1;
                    else
                        b(inode) = 0;
                    end
                end
            end
        end
        
    end
    
end

