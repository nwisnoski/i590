setwd('~/i590/ind_proj/pairwise/')
chr1.fixed = read.table('chr1_all_win_res_fixed_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr2.fixed = read.table('chr2_all_win_res_fixed_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr3.fixed = read.table('chr3_all_win_res_fixed_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr4.fixed = read.table('chr4_all_win_res_fixed_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr22.fixed = read.table('chr22_all_win_res_fixed_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
fixed = rbind(chr1.fixed, chr2.fixed, chr3.fixed, chr4.fixed, chr22.fixed)

chr1.out = read.table('chr1_all_win_res_outside_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr2.out = read.table('chr2_all_win_res_outside_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr3.out = read.table('chr3_all_win_res_outside_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr4.out = read.table('chr4_all_win_res_outside_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
chr22.out = read.table('chr22_all_win_res_outside_pairwise_overlap.txt', sep=' ', col.names=c('win', 'rho'))
out = rbind(chr1.out, chr2.out, chr3.out, chr4.out, chr22.out)

chr1.pol = read.table('chr1_all_win_res_polymorphic_pairwise_rho_overlap.txt')
chr2.pol = read.table('chr2_all_win_res_polymorphic_pairwise_rho_overlap.txt')
chr3.pol = read.table('chr3_all_win_res_polymorphic_pairwise_rho_overlap.txt')
chr4.pol = read.table('chr4_all_win_res_polymorphic_pairwise_rho_overlap.txt')
chr22.pol = read.table('chr22_all_win_res_polymorphic_pairwise_rho_overlap.txt')
pol = rbind(chr1.pol, chr2.pol, chr3.pol, chr4.pol, chr22.pol)
names(pol) = c('id', 'chr', 'pos', 'ref', 'alt', 'strand', 'type', 'start', 'end', 'refb', 'ref_freq', 'altb', 'alt_freq', 'daf', 'rho')

corr.1000.f.1 = read.table('chr1_all_win_res_fixed_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.f.2 = read.table('chr2_all_win_res_fixed_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.f.3 = read.table('chr3_all_win_res_fixed_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.f.4 = read.table('chr4_all_win_res_fixed_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.f.22 = read.table('chr22_all_win_res_fixed_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.f = rbind(corr.1000.f.1, corr.1000.f.2, corr.1000.f.3, corr.1000.f.4, corr.1000.f.22)

corr.2000.f.1 = read.table('chr1_all_win_res_fixed_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.f.2 = read.table('chr2_all_win_res_fixed_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.f.3 = read.table('chr3_all_win_res_fixed_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.f.4 = read.table('chr4_all_win_res_fixed_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.f.22 = read.table('chr22_all_win_res_fixed_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.f = rbind(corr.2000.f.1, corr.2000.f.2, corr.2000.f.3, corr.2000.f.4, corr.2000.f.22)

corr.4000.f.1 = read.table('chr1_all_win_res_fixed_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.f.2 = read.table('chr2_all_win_res_fixed_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.f.3 = read.table('chr3_all_win_res_fixed_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.f.4 = read.table('chr4_all_win_res_fixed_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.f.22 = read.table('chr22_all_win_res_fixed_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.f = rbind(corr.4000.f.1, corr.4000.f.2, corr.4000.f.3, corr.4000.f.4, corr.4000.f.22)

corr.8000.f.1 = read.table('chr1_all_win_res_fixed_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.f.2 = read.table('chr2_all_win_res_fixed_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.f.3 = read.table('chr3_all_win_res_fixed_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.f.4 = read.table('chr4_all_win_res_fixed_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.f.22 = read.table('chr22_all_win_res_fixed_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.f = rbind(corr.8000.f.1, corr.8000.f.2, corr.8000.f.3, corr.8000.f.4, corr.8000.f.22)

corr.1000.o.1 = read.table('chr1_all_win_res_outside_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.o.2 = read.table('chr2_all_win_res_outside_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.o.3 = read.table('chr3_all_win_res_outside_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.o.4 = read.table('chr4_all_win_res_outside_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.o.22 = read.table('chr22_all_win_res_outside_pairwise_overlap_corr_1000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.1000.o = rbind(corr.1000.o.1, corr.1000.o.2, corr.1000.o.3, corr.1000.o.4, corr.1000.o.22)

corr.2000.o.1 = read.table('chr1_all_win_res_outside_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.o.2 = read.table('chr2_all_win_res_outside_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.o.3 = read.table('chr3_all_win_res_outside_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.o.4 = read.table('chr4_all_win_res_outside_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.o.22 = read.table('chr22_all_win_res_outside_pairwise_overlap_corr_2000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.2000.o = rbind(corr.2000.o.1, corr.2000.o.2, corr.2000.o.3, corr.2000.o.4, corr.2000.o.22)

corr.4000.o.1 = read.table('chr1_all_win_res_outside_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.o.2 = read.table('chr2_all_win_res_outside_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.o.3 = read.table('chr3_all_win_res_outside_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.o.4 = read.table('chr4_all_win_res_outside_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.o.22 = read.table('chr22_all_win_res_outside_pairwise_overlap_corr_4000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.4000.o = rbind(corr.4000.o.1, corr.4000.o.2, corr.4000.o.3, corr.4000.o.4, corr.4000.o.22)

corr.8000.o.1 = read.table('chr1_all_win_res_outside_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.o.2 = read.table('chr2_all_win_res_outside_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.o.3 = read.table('chr3_all_win_res_outside_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.o.4 = read.table('chr4_all_win_res_outside_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.o.22 = read.table('chr22_all_win_res_outside_pairwise_overlap_corr_8000.txt', sep='\t', col.names=c('win', 'rho', 'nearby.wins'))
corr.8000.o = rbind(corr.8000.o.1, corr.8000.o.2, corr.8000.o.3, corr.8000.o.4, corr.8000.o.22)

pdf('~/i590/ind_proj/plots/1_2_3_4_22_no_c_d_f_phased_overlap.pdf', height=5, width=7)
par(mar=c(8, 4, 2, 0.3))
boxplot(list(out$rho,
             corr.1000.o$rho[corr.1000.o$nearby.wins>=1],
             corr.1000.o$rho[corr.1000.o$nearby.wins>=2],
             corr.1000.o$rho[corr.1000.o$nearby.wins>=3],
             corr.1000.o$rho[corr.1000.o$nearby.wins>=4],
             pol$rho[pol$type=='d'],
             pol$rho[pol$type=='r'], fixed$rho,
             corr.1000.f$rho[corr.1000.f$nearby.wins>=1],
             corr.1000.f$rho[corr.1000.f$nearby.wins>=2],
             corr.1000.f$rho[corr.1000.f$nearby.wins>=3],
             corr.1000.f$rho[corr.1000.f$nearby.wins>=4]),
        outline=F,
        names=c(paste('No motif',  paste0('(', length(out$rho), ')')),
          paste('No motif 1', paste0('(', length(corr.1000.o$rho[corr.1000.o$nearby.wins>=1]), ')')),
          paste('No motif 2', paste0('(', length(corr.1000.o$rho[corr.1000.o$nearby.wins>=2]), ')')),
          paste('No motif 3', paste0('(', length(corr.1000.o$rho[corr.1000.o$nearby.wins>=3]), ')')),
          paste('No motif 4', paste0('(', length(corr.1000.o$rho[corr.1000.o$nearby.wins>=4]), ')')),
          paste('Create', paste0('(', length(pol$rho[pol$type=='r']), ')')),
          paste('Disrupt', paste0('(', length(pol$rho[pol$type=='d']), ')')),
          paste('Fixed', paste0('(', length(fixed$rho), ')')),
          paste('Fixed 1', paste0('(', length(corr.1000.f$rho[corr.1000.f$nearby.wins>=1]), ')')),
          paste('Fixed 2', paste0('(', length(corr.1000.f$rho[corr.1000.f$nearby.wins>=2]), ')')),
          paste('Fixed 3', paste0('(', length(corr.1000.f$rho[corr.1000.f$nearby.wins>=3]), ')')),
          paste('Fixed 4', paste0('(', length(corr.1000.f$rho[corr.1000.f$nearby.wins>=4]), ')'))),
        las=3, ylab='Rho', main='(No motif & Fixed 1-4: within 1 kb)')
dev.off()

pol$bin[pol$daf <= 0.14] = 1
pol$bin[pol$daf > 0.14 & pol$daf <= 0.23] = 2
pol$bin[pol$daf > 0.23 & pol$daf <= 0.32] = 3
pol$bin[pol$daf > 0.32 & pol$daf <= 0.41] = 4
pol$bin[pol$daf > 0.41 & pol$daf <= 0.50] = 5
pol$bin[pol$daf > 0.50 & pol$daf <= 0.59] = 6
pol$bin[pol$daf > 0.59 & pol$daf <= 0.68] = 7
pol$bin[pol$daf > 0.68 & pol$daf <= 0.77] = 8
pol$bin[pol$daf > 0.77 & pol$daf <= 0.86] = 9
pol$bin[pol$daf > 0.86 & pol$daf <= 0.95] = 10

pol.c = pol[pol$type=='r',]
pol.c$type=rep('c', length(pol.c$type))
pol.d = pol[pol$type=='d',]
pol.d$type=rep('d', length(pol.d$type))

nc.1 = length(pol.c[pol.c$bin==1,2])
nc.2 = length(pol.c[pol.c$bin==2,2])
nc.3 = length(pol.c[pol.c$bin==3,2])
nc.4 = length(pol.c[pol.c$bin==4,2])
nc.5 = length(pol.c[pol.c$bin==5,2])
nc.6 = length(pol.c[pol.c$bin==6,2])
nc.7 = length(pol.c[pol.c$bin==7,2])
nc.8 = length(pol.c[pol.c$bin==8,2])
nc.9 = length(pol.c[pol.c$bin==9,2])
nc.10 = length(pol.c[pol.c$bin==10,2])

nd.1 = length(pol.d[pol.d$bin==1,2])
nd.2 = length(pol.d[pol.d$bin==2,2])
nd.3 = length(pol.d[pol.d$bin==3,2])
nd.4 = length(pol.d[pol.d$bin==4,2])
nd.5 = length(pol.d[pol.d$bin==5,2])
nd.6 = length(pol.d[pol.d$bin==6,2])
nd.7 = length(pol.d[pol.d$bin==7,2])
nd.8 = length(pol.d[pol.d$bin==8,2])
nd.9 = length(pol.d[pol.d$bin==9,2])
nd.10 = length(pol.d[pol.d$bin==10,2])

pdf('~/i590/ind_proj/plots/1_2_3_4_22_c_pairwise_phased_overlap.pdf')
boxplot(list(pol.c$rho[pol.c$bin == 1], pol.c$rho[pol.c$bin == 2], pol.c$rho[pol.c$bin == 3], pol.c$rho[pol.c$bin == 4], pol.c$rho[pol.c$bin == 5], pol.c$rho[pol.c$bin == 6], pol.c$rho[pol.c$bin == 7], pol.c$rho[pol.c$bin == 8], pol.c$rho[pol.c$bin == 9], pol.c$rho[pol.c$bin == 10]), names=c(paste('0.05 - 0.14', paste0('(', nc.1, ')'), sep='\n'), paste('0.14 - 0.23', paste0('(', nc.2, ')'), sep='\n'), paste('0.23 - 0.32', paste0('(', nc.3, ')'), sep='\n'), paste('0.32 - 0.41', paste0('(', nc.4, ')'), sep='\n'), paste('0.41 - 0.5', paste0('(', nc.5, ')'), sep='\n'), paste('0.5 - 0.59', paste0('(', nc.6, ')'), sep='\n'), paste('0.59 - 0.68', paste0('(', nc.7, ')'), sep='\n'), paste('0.68 - 0.77', paste0('(', nc.8, ')'), sep='\n'), paste('0.77 - 0.86', paste0('(', nc.9, ')'), sep='\n'), paste('0.86 - 0.95', paste0('(', nc.10, ')'), sep='\n')), ylab='Mean SNP rho', las=2, main='SNPs that create motifs', outline=F)
dev.off()

pdf('~/i590/ind_proj/plots/1_2_3_4_22_d_pairwise_phased_overlap.pdf')
boxplot(list(pol.d$rho[pol.d$bin == 1], pol.d$rho[pol.d$bin == 2], pol.d$rho[pol.d$bin == 3], pol.d$rho[pol.d$bin == 4], pol.d$rho[pol.d$bin == 5], pol.d$rho[pol.d$bin == 6], pol.d$rho[pol.d$bin == 7], pol.d$rho[pol.d$bin == 8], pol.d$rho[pol.d$bin == 9], pol.d$rho[pol.d$bin == 10]), names=c(paste('0.05 - 0.14', paste0('(', nd.1, ')'), sep='\n'), paste('0.14 - 0.23', paste0('(', nd.2, ')'), sep='\n'), paste('0.23 - 0.32', paste0('(', nd.3, ')'), sep='\n'), paste('0.32 - 0.41', paste0('(', nd.4, ')'), sep='\n'), paste('0.41 - 0.5', paste0('(', nd.5, ')'), sep='\n'), paste('0.5 - 0.59', paste0('(', nd.6, ')'), sep='\n'), paste('0.59 - 0.68', paste0('(', nd.7, ')'), sep='\n'), paste('0.68 - 0.77', paste0('(', nd.8, ')'), sep='\n'), paste('0.77 - 0.86', paste0('(', nd.9, ')'), sep='\n'), paste('0.86 - 0.95', paste0('(', nd.10, ')'), sep='\n')), ylab='Mean SNP rho', las=2, main='SNPs that disrupt motifs', outline=F)
dev.off()

fixed.df = data.frame(fixed$rho, rep('f', length(fixed$rho)))
names(fixed.df) = c('rho', 'type')
out.df = data.frame(out$rho, rep('o', length(out$rho)))
names(out.df) = c('rho', 'type')
pol.c.df = data.frame(pol.c$rho, pol.c$type)
names(pol.c.df) = c('rho', 'type')
pol.d.df = data.frame(pol.d$rho, pol.d$type)
names(pol.d.df) = c('rho', 'type')ta
all = rbind(fixed.df, out.df, pol.c.df, pol.d.df)

kruskal.test(rho~type, data=all)
kruskalmc(all$rho, all$type)

kruskal.test(rho~bin, data=pol.c)
kruskalmc(pol.c$bin, pol.c$rho)
kruskal.test(rho~bin, data=pol.d)

pdf('~/i590/ind_proj/plots/no_motif_hist.pdf')
hist(out$rho, nclass=50, main='No motifs (YRI, 86 individuals)', xlab='Rho for 20 SNP-windows')
dev.off()

pdf('~/i590/ind_proj/plots/fixed_motif_hist.pdf')
hist(fixed$rho, nclass=50, main='Fixed motifs (YRI, 86 individuals)', xlab='Rho for 20 SNP-windows')
dev.off()

pdf('~/i590/ind_proj/plots/create_motif_hist.pdf')
hist(pol.c$rho, nclass=50, main='SNPs that create motifs (YRI, 86 individuals)', xlab='Rho for 21 SNP-windows')
dev.off()

pdf('~/i590/ind_proj/plots/disrupt_motif_hist.pdf')
hist(pol.d$rho, nclass=50, main='SNPs that disrupt motifs (YRI, 86 individuals)', xlab='Rho for 21 SNP-windows')
dev.off()

hist(corr.2000.f$nearby.wins, nclass=max(corr.2000.f$nearby.wins), main='Fixed: number of motif windows within 2kb')
hist(corr.4000.f$nearby.wins, nclass=max(corr.4000.f$nearby.wins), main='Fixed: number of motif windows within 4kb')
hist(corr.8000.f$nearby.wins, nclass=max(corr.8000.f$nearby.wins), main='Fixed: number of motif windows within 8kb')
plot(corr.2000.f$rho ~ corr.2000.f$nearby.wins)
plot(corr.4000.f$rho ~ corr.4000.f$nearby.wins)
plot(corr.8000.f$rho ~ corr.8000.f$nearby.wins)

linear.model.1000.o = lm(corr.1000.o$rho ~ corr.1000.o$nearby.wins)
summary(linear.model.1000.o)

linear.model.2000.o = lm(corr.2000.o$rho ~ corr.2000.o$nearby.wins)
summary(linear.model.2000.o)

linear.model.4000.o = lm(corr.4000.o$rho ~ corr.4000.o$nearby.wins)
summary(linear.model.4000.o)

linear.model.8000.o = lm(corr.8000.o$rho ~ corr.8000.o$nearby.wins)
summary(linear.model.8000.o)

linear.model.1000.f = lm(corr.1000.f$rho ~ corr.1000.f$nearby.wins)
summary(linear.model.1000.f)

linear.model.2000.f = lm(corr.2000.f$rho ~ corr.2000.f$nearby.wins)
summary(linear.model.2000.f)

linear.model.4000.f = lm(corr.4000.f$rho ~ corr.4000.f$nearby.wins)
summary(linear.model.4000.f)

linear.model.8000.f = lm(corr.8000.f$rho ~ corr.8000.f$nearby.wins)
summary(linear.model.8000.f)

pdf('~/i590/ind_proj/plots/correlation_rho_win_proximity_outside.pdf', width=10, height=10)
par(mfrow=c(2,2))
plot(corr.1000.o$nearby.wins, corr.1000.o$rho, xlab = 'Number of windows within 1 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.1000.o)[2]))
abline(linear.model.2000.o)
plot(corr.2000.o$nearby.wins, corr.2000.o$rho, xlab = 'Number of windows within 2 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.2000.o)[2]))
abline(linear.model.2000.o)
plot(corr.4000.o$nearby.wins, corr.4000.o$rho, xlab = 'Number of windows within 4 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.4000.o)[2]))
abline(linear.model.4000.o)
plot(corr.8000.o$nearby.wins, corr.8000.o$rho, xlab = 'Number of windows within 8 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.8000.o)[2]))
abline(linear.model.8000.o)
dev.off()

pdf('~/i590/ind_proj/plots/correlation_rho_win_proximity_fixed.pdf', width=10, height=10)
par(mfrow=c(2,2))
plot(corr.1000.f$nearby.wins, corr.1000.f$rho, xlab = 'Number of windows within 1 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.1000.f)[2]))
abline(linear.model.1000.f)
plot(corr.2000.f$nearby.wins, corr.2000.f$rho, xlab = 'Number of windows within 2 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.2000.f)[2]))
abline(linear.model.2000.f)
plot(corr.4000.f$nearby.wins, corr.4000.f$rho, xlab = 'Number of windows within 4 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.4000.f)[2]))
abline(linear.model.4000.f)
plot(corr.8000.f$nearby.wins, corr.8000.f$rho, xlab = 'Number of windows within 8 kb', ylab = 'Rho', main=paste('coefficient: ', coef(linear.model.8000.f)[2]))
abline(linear.model.8000.f)
dev.off()
