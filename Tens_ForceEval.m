clc
clear all
close all
directory ='D:\BiAx\20%\21_2_24\Set3\Processing\';
directory_save = directory;
images = dir([directory,'*.png']);
files = dir([directory,'*solved.mat']);
for ii = 1:numel(files)
    load([directory,files(ii).name]);
    neighbond = [];
    signstore = [];
    fx = 0;
    fy = 0;
    particle = pres;
    bond = [particle.bond];
    par = particle(bond);
    erase = ([par.z]' == 0);

    b = cat(1,par.id);
    for k = 1:numel(particle)
        if bond(k) == false
            continue
        end
        beta = [particle(k).betas];
        f = zeros(size(beta));
        alphanew = zeros(size(beta));
        neigh = particle(k).neighbours;
        for i = 1:length(f)
            if (neigh(i) < 0)
                continue

            end
            ind = find(particle(neigh(i)).neighbours ==k);
            f(i) = particle(neigh(i)).forces(ind);
            alphanew(i) = particle(neigh(i)).alphas(ind);
            f_cell = num2cell(f);
            [particle(:).forces] = deal(f_cell{:});
        end
        particle(k).forces = f;
        particle(k).alphas = alphanew;
    end
    par2 = particle(bond);


    for i = 1:numel(particle)
        m = isempty(particle(i).parBonded);
        if bond(i) == false || m == true

            continue
        end
        c = isempty(find([particle(i).neighbours] == particle(i).parBonded));
        if c == true
            continue
        end
        sum1=0;
        sum2=0;
        fn=0;
        fx = 0;
        fy = 0;
        fnew = 0;
        ind = find(particle(i).neighbours == particle(i).parBonded);
        for j = 1:particle(i).z
            if particle(i).neighbours(j) == particle(i).parBonded
                storebeta = particle(i).betas(j);
                continue
            end
            b = particle(i).betas(j);
            a = particle(i).alphas(j);
            angle = inclination(a,b);
            angle = particle(i).alphas(j) - particle(i).betas(j);
            fx = particle(i).forces(j)*cos(angle) +fx;
            fy = particle(i).forces(j)*sin(angle) +fy;
            sum1 = sum1 - particle(i).forces(j)*cos(alpha(j)+beta(i)-beta(j));
            sum2 = sum2 - particle(i).forces(j)*sin(alpha(j)+beta(i)-beta(j));
        end
        fn = sqrt(sum1^2 + sum2^2);
        incl2=atan2(-fy,-fx);
        fnew = sqrt(fx^2 + fy^2);
        inclin = atan2(-fy,-fx);
        index = find(particle(i).neighbours == particle(i).parBonded);
        particle(i).forces(index) =  fnew;
        particle(i).alphas(index) = attackangle(storebeta,inclin);

    end
    for k = 1: numel(particle)
        if bond(k) == false
            continue
        end

        forces = [particle(k).forces]';
        alpha = [particle(k).alphas]';
        beta = particle(k).betas;
        fsigma = particle.fsigma;
        radius_in_meter = particle(k).rm;
        p = [forces;alpha];
        optim_method = 1;
        scaling = 1;
        particle_template = imresize(particle(k).forceImage, scaling);
        particle_template_size = size(particle_template,1);
        imgFit = joForceImg(p, beta, fsigma, radius_in_meter,...
            particle_template_size*(1/scaling), optim_method);
        particle(k).synthImg = imgFit;
    end
    pres = particle;
    save([directory_save, files(ii).name(1:end-4),'.mat'],'pres');
end






