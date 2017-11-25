classdef HeatSquare9Point < HeatSquareBase
    % Heat conduction for 9 point discretisation
    % Author: Sebastia Ramon
    
    methods
        
        function A = assemble_A(obj)
            if obj.Nx ~= obj.Ny
                error('Nx is not equal to Ny');
            end
            N = obj.Nx;
            Delta_x = obj.Delta_x;
            A = sparse(zeros(N * N, N * N));
            a = 2 / Delta_x / Delta_x;
            b = 1 / Delta_x / Delta_x;
            c = -12 / Delta_x / Delta_x;
            for i=1:N
                for j=1:N
                    inode = N * (i-1) + j;
                    if i == 1
                        A(inode, inode) = 1;
                    elseif i == N
                        A(inode, inode) = 1;
                    elseif j == 1
                        A(inode, inode) = 1;
                    elseif j == N
                        A(inode, inode) = 1;
                    else
                        A(inode, inode - N - 1) = b;
                        A(inode, inode - N) = a;
                        A(inode, inode - N + 1) = b;
                        A(inode, inode - 1) = a;
                        A(inode, inode) = c;
                        A(inode, inode + 1) = a;
                        A(inode, inode + N - 1) = b;
                        A(inode, inode + N) = a;
                        A(inode, inode + N + 1) = b;
                    end
                end
            end
        end
        
        function b = assemble_b(obj)
            if obj.Nx ~= obj.Ny
                error('Nx is not equal to Ny');
            end
            N = obj.Nx;
            Delta_x = obj.Delta_x;
            b = zeros(N * N, 1);
            for i=1:N
                for j=1:N
                    inode = N * (i-1) + j;
                    if i == 1
                        x = (j - 1) * Delta_x;
                        b(inode) = 1 + x;
                    elseif i == N
                        b(inode) = 1;
                    elseif j == 1
                        b(inode) = 1;
                    elseif j == N
                        y = (i - 1) * Delta_x;
                        b(inode) = cos(6 * 3/2 * pi * y) + 1;
                    else
                        b(inode) = 0;
                    end
                end
            end
        end
        
    end
    
end

