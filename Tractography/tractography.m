function tractography(ellipsoideMatrix,I0,it,step)

if nargin==3
    step = 0.2;
end

slice(I0,[],[],30);
hold on;
%slice(I0,[],[],55);
        
for nbLoop=1:it
    nbLoop
    %seed = [39 44 35];
    %seed = [randi([30 50],1,1) randi([34 54],1,1) randi([30 40],1,1)]
    seed = [randi([20 60],1,1) randi([20 70],1,1) randi([20 50],1,1)]
    
    if(I0(seed) < 150)
        I0(seed)
        continue;
    end
    beginSeed = seed;
    
    previousNormalizeVector = [0 0 0];

    while 1
        
        % Interpolation du nouveau point
        for i=1:3
            for j=1:3
                interpVec(i,j) = interp3(squeeze(ellipsoideMatrix(i,j,:,:,:)), seed(2),seed(1),seed(3));
            end
        end
        
        if isnan(interpVec) == 1
            break;
        end
        % Calcul des valeurs et vecteurs propres
        [V,Diag] = eig(interpVec);
        eigenValue = [Diag(1,1) Diag(2,2) Diag(3,3)];
        eigenVector = transpose(V(:,3));

        % Calcul de la position de la nouvelle position de la seed
        normalizeVector = eigenVector/norm(eigenVector);
        
        if dot(previousNormalizeVector,normalizeVector) < 0
            normalizeVector = -normalizeVector;
        end
        
        newSeed = seed + step*normalizeVector;
        
        % Affichage de la fibre

        plot3([seed(2) newSeed(2)], [seed(1) newSeed(1)], [seed(3) newSeed(3)], 'Color', abs(normalizeVector));
        drawnow;

        % Calcul de l'angle entre le vector propre courant et le précedent
        crossProduct = cross(previousNormalizeVector,normalizeVector);
        dotProduct = dot(previousNormalizeVector,normalizeVector);
        ThetaInDegrees = atan2d( norm(crossProduct),dotProduct );

        % Si l'angle est supérieur à 45 degré on recalcul l'angle avec un vector
        % propre courant dans l'autre sens
        if ThetaInDegrees > 45
%             normalizeVector = normalizeVector*-1;
%             crossProduct = cross(previousNormalizeVector,normalizeVector);
%             dotProduct = dot(previousNormalizeVector,normalizeVector);
%             ThetaInDegrees = atan2d( norm(crossProduct),dotProduct );
% 
%             % Si l'angle est toujours supérieur à 45 degré on stop la tractography
%             % de la fibre à cette extrémité
%             if ThetaInDegrees > 45
               break;
            %end
        end

        % Calcul du Facteur d'Anisotropie (FA)
        MD = (eigenValue(1)+eigenValue(2)+eigenValue(3)) / 3;

        num   = sqrt((eigenValue(1)-MD).^2 + (eigenValue(2)-MD).^2 + (eigenValue(3)-MD).^2);
        denum = sqrt(eigenValue(1).^2 + eigenValue(2).^2 + eigenValue(3).^2);

        FA = sqrt(3/2) * num/denum;

        if FA < 0.2
           break; 
        end
        
        % Mise à jour des données
        previousNormalizeVector = normalizeVector;
        seed = newSeed;
        
    end
    
    
    SecondSense = true;
    first = true;
    seed = beginSeed;
    previousNormalizeVector = [0 0 0];
    
    while 1
        
        % Interpolation du nouveau point
        for i=1:3
            for j=1:3
                interpVec(i,j) = interp3(squeeze(ellipsoideMatrix(i,j,:,:,:)), seed(2),seed(1),seed(3));
            end
        end
        
        % Calcul des valeurs et vecteurs propres
        if isnan(interpVec) == 1
            break;
        end
        [V,Diag] = eig(interpVec);
        eigenValue = [Diag(1,1) Diag(2,2) Diag(3,3)];
        eigenVector = transpose(V(:,3));

        % Calcul de la position de la nouvelle position de la seed
        normalizeVector = eigenVector/norm(eigenVector);
        
        if dot(previousNormalizeVector,normalizeVector) < 0 && first == false
            normalizeVector = -normalizeVector;
        end
        
        if first == true
            normalizeVector = -normalizeVector;
            first = false;
        end
        

        
        newSeed = seed + step*normalizeVector;
        
        % Affichage de la fibre

        plot3([seed(2) newSeed(2)], [seed(1) newSeed(1)], [seed(3) newSeed(3)], 'Color', abs(normalizeVector));
        drawnow;

        % Calcul de l'angle entre le vector propre courant et le précedent
        crossProduct = cross(previousNormalizeVector,normalizeVector);
        dotProduct = dot(previousNormalizeVector,normalizeVector);
        ThetaInDegrees = atan2d( norm(crossProduct),dotProduct );

        % Si l'angle est supérieur à 45 degré on recalcul l'angle avec un vector
        % propre courant dans l'autre sens
        if ThetaInDegrees > 45
%             normalizeVector = normalizeVector*-1;
%             crossProduct = cross(previousNormalizeVector,normalizeVector);
%             dotProduct = dot(previousNormalizeVector,normalizeVector);
%             ThetaInDegrees = atan2d( norm(crossProduct),dotProduct );
% 
%             % Si l'angle est toujours supérieur à 45 degré on stop la tractography
%             % de la fibre à cette extrémité
%             if ThetaInDegrees > 45
               break;
            %end
        end

        % Calcul du Facteur d'Anisotropie (FA)
        MD = (eigenValue(1)+eigenValue(2)+eigenValue(3)) / 3;

        num   = sqrt((eigenValue(1)-MD).^2 + (eigenValue(2)-MD).^2 + (eigenValue(3)-MD).^2);
        denum = sqrt(eigenValue(1).^2 + eigenValue(2).^2 + eigenValue(3).^2);

        FA = sqrt(3/2) * num/denum;

        if FA < 0.2
           break; 
        end
        
        % Mise à jour des données
        previousNormalizeVector = normalizeVector;
        seed = newSeed;
        
    end
    %drawnow;
end