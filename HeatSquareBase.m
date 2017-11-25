classdef HeatSquareBase
    % Heat conduction base class
    % Author: Sebastia Ramon
    
    properties
        L = 1
        Nx = 80
        Ny = 80
    end
    properties (Dependent)
        Delta_x
        Delta_y
    end
    
    methods
        
        function Delta_x = get.Delta_x(obj)
            Delta_x = obj.L / (obj.Nx - 1);
        end
        
        function Delta_y = get.Delta_y(obj)
            Delta_y = obj.L / (obj.Ny - 1);
        end
        
        function T = solve(obj, plot)
            A = obj.assemble_A();
            b = obj.assemble_b();
            xVec = obj.Delta_x * (0:(obj.Nx - 1));
            yVec = obj.Delta_y * (0:(obj.Ny - 1));
            T = A\b;
            T = reshape(T, [obj.Ny, obj.Nx]);
            T = permute(T, [2, 1]);
            if exist('plot','var')
                if plot ~= 0
                    figure;
                    surf(xVec, yVec, T);
                    xlabel('x');
                    ylabel('y');
                    legend('T');
                    figure;
                    contour(xVec, yVec, T, 50);
                    xlabel('x');
                    ylabel('y');
                    legend('T');
                end
            end
        end
        
        function e = findAvgError(obj, size, sizeExact, exact)
            obj.Nx = size + 1;
            obj.Ny = size + 1;
            inexact = obj.solve();
            multiplier = sizeExact/size;
            e = 0;
            count = 0;
            for i=1:size+1
                iexact = (i-1)*multiplier + 1;
                for j=1:size+1
                    jexact = (j-1)*multiplier + 1;
                    error = abs(inexact(i,j)-exact(iexact,jexact));
                    e = e + error;
                    count = count + 1;
                end
            end
            e = e / count;
        end
        
        function e = findMaxError(obj, size, sizeExact, exact)
            obj.Nx = size + 1;
            obj.Ny = size + 1;
            inexact = obj.solve();
            multiplier = sizeExact/size;
            e = 0;
            for i=1:size+1
                iexact = (i-1)*multiplier + 1;
                for j=1:size+1
                    jexact = (j-1)*multiplier + 1;
                    error = abs(inexact(i,j)-exact(iexact,jexact));
                    if error > e
                        e = error;
                    end
                end
            end
        end
        
        function [errorVec, dispVec] = error(obj)
            % tic
            sizeExact = 128;
            sizeInexact = [4,8,16,32,64];
            dispVec = obj.L ./ sizeInexact;
            nErr = size(sizeInexact,2);
            obj.Nx = sizeExact + 1;
            obj.Ny = sizeExact + 1;
            exact = obj.solve();
            errorVec = zeros(1, nErr);
            for i = 1:nErr
                errorVec(i) = obj.findAvgError(sizeInexact(i), sizeExact, exact);
            end
            % toc
            loglog(dispVec, errorVec);
            polyfit(log(dispVec), log(errorVec), 1)
        end
        
    end
    
end

