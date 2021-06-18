tic
addpath('./scGEAToolbox-master');

mtxf='./matrix.mtx'; %read expression matrix
genf='./gene.tsv'; %read gene name
[X,genelist]=sc_readmtxfile(mtxf,genf,[],2);

[Xnorm]=sc_norm(X,'type','deseq');
[T,~]=sc_splinefit(Xnorm,genelist,true);
writetable(T,"scGEAToolbox_l_default.csv") %save result
toc