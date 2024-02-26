library("Rsubread")
library(edgeR)
library(Rsubread)



mouse_psi.df<-data.frame(mouse_psi)
rownames(mouse_psi.df)<-mouse_psi$Gene

mouse_psi.df<-subset(mouse_psi.df, select=-Gene)

expr.design <- read.table('experiment_design.txt', header=T, sep='\t')
rownames(expr.design) <- expr.design$SampleID

#order the design in the same ordering as the counts object
expr.design <- expr.design[colnames(mouse_psi.df),]

expr.design

samples <- as.character(expr.design$SampleID)
group <- factor(expr.design$group)

library(edgeR)

dim(mouse_psi.df)

design <- model.matrix(~0 + group)
colnames(design) <- c(levels(group))
design


library(limma)
group.colours <- c('grey30','red4','red')[group];
par(mar=c(5.1, 5, 4.1, 7), xpd=TRUE)

dpsi <- DGEList(mouse_psi.df)


boxplot(dge$counts, 
        col=group.colours,
        main="PSI Variability",
        xlab="",
        ylab="PSI values",
        las=2,cex.axis=0.8)
####################################################
contrast.psi<- makeContrasts(WT - C, levels = design) #G298S - C, G298S - WT


fit <- lmFit(mouse_psi.df, design)
fit2 <- contrasts.fit(fit, contrast.psi)
fit2 <- eBayes(fit2, trend=TRUE)
fit2.decide<-decideTests(fit2)
summary(fit2.decide)

library(EnhancedVolcano)
top.table <- topTable(fit2, sort.by = "P", adjust.method="BH",n = Inf)
top15<-c(head(rownames(top.table), 15))
EnhancedVolcano(top.table,
                lab = rownames(top.table),
                x = 'logFC',
                y = 'P.Value',
                selectLab = c(top15),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 7,
                #labFace = 'bold',
                #boxedLabels = TRUE,
                legendPosition = 'bottom',
                legendLabSize = 5,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                col=c('grey', 'grey', 'grey', 'red4'),
                colAlpha = 1,
                widthConnectors = 0.5)

######################################
contrast.psi<- makeContrasts(G298S - C, levels = design) #G298S - C, G298S - WT


fit <- lmFit(mouse_psi.df, design)
fit2 <- contrasts.fit(fit, contrast.psi)
fit2 <- eBayes(fit2, trend=TRUE)

fit2.decide<-decideTests(fit2)
summary(fit2.decide)

library(EnhancedVolcano)
top.table <- topTable(fit2, sort.by = "P", adjust.method="BH",n = Inf)
top15<-c(head(rownames(top.table), 15))
EnhancedVolcano(top.table,
                lab = rownames(top.table),
                x = 'logFC',
                y = 'P.Value',
                selectLab = c(top15),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 7,
                #labFace = 'bold',
                #boxedLabels = TRUE,
                legendPosition = 'bottom',
                legendLabSize = 5,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                col=c('grey', 'grey', 'grey', 'red'),
                colAlpha = 1,
                widthConnectors = 0.5)
##############################################
contrast.psi<- makeContrasts(G298S - WT, levels = design) #G298S - C, G298S - WT


fit <- lmFit(mouse_psi.df, design)
fit2 <- contrasts.fit(fit, contrast.psi)
fit2 <- eBayes(fit2, trend=TRUE)

fit2.decide<-decideTests(fit2)
summary(fit2.decide)

library(EnhancedVolcano)
top.table <- topTable(fit2, sort.by = "P", adjust.method="BH",n = Inf)
top15<-c(head(rownames(top.table), 15))
EnhancedVolcano(top.table,
                lab = rownames(top.table),
                x = 'logFC',
                y = 'P.Value',
                selectLab = c(top15),
                pCutoff = 10e-2,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 7,
                #labFace = 'bold',
                #boxedLabels = TRUE,
                legendPosition = 'bottom',
                legendLabSize = 5,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                col=c('grey', 'grey', 'grey', 'blue'),
                colAlpha = 1,
                widthConnectors = 0.5)
