# Part C: population-referenced Kinhom and ZIP income map.
# Run from assignment4/ or assignment4/code/. All distances are US survey feet.
run_hw4_part_c <- function(root = c('.', '..')) {
  suppressPackageStartupMessages({library(sf); library(spatstat.geom)
    library(spatstat.explore); library(spatstat.random); library(ggplot2)})
  root <- root[file.exists(file.path(root, 'data', 'Philadelphia.shp'))][1]
  if (is.na(root)) stop('Run from assignment4/ or assignment4/code/.')
  root <- normalizePath(root)
  inp <- file.path(root,'data'); out <- file.path(root,'output')
  proc <- file.path(inp,'processed_kinhom')
  dir.create(proc,recursive=TRUE,showWarnings=FALSE)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  if (!file.exists(file.path(inp,'processed_prepare','markets_ppp.Rds'))) {
    source(file.path(root,'code','HW4_prepare_data.R'),local=TRUE)
    prepare_hw4_data(root)
  }
  X <- readRDS(file.path(inp,'processed_prepare','markets_ppp.Rds'))
  W <- readRDS(file.path(inp,'processed_prepare','city_owin.Rds'))
  z <- st_read(file.path(inp,'processed_prepare','zip_sf.gpkg'),quiet=TRUE)
  city <- st_read(file.path(inp,'Philadelphia.shp'),quiet=TRUE)
  markets <- st_read(file.path(inp,'Philadelphia_Farmers_Markets201302.shp'),quiet=TRUE)
  stopifnot(npoints(X)==62L,st_crs(z)$epsg==2272L)
  # ZIP density uses the full supplied ZIP area; do not invent a population
  # for zero-population ZIP 19112 or remove it from Philadelphia's window.
  z$area_ft2 <- as.numeric(st_area(z))
  z$population_density <- z$Pop2000/z$area_ft2
  mask <- as.mask(W,eps=250)
  grid <- expand.grid(x_ft=mask$xcol,y_ft=mask$yrow)
  inside <- inside.owin(grid$x_ft,grid$y_ft,W)
  gsf <- st_as_sf(grid[inside,],coords=c('x_ft','y_ft'),crs=2272)
  hits <- st_intersects(gsf,z)
  iz <- vapply(hits,function(v) if(length(v)) v[1] else NA_integer_,integer(1))
  # Only fill tiny city/ZIP boundary mismatches, recording the cell count.
  missing <- is.na(iz)
  if(any(missing)) iz[missing] <- st_nearest_feature(gsf[missing,],z)
  raw <- rep(NA_real_,nrow(grid)); raw[inside] <- z$population_density[iz]
  pop_im <- im(matrix(raw,nrow=length(mask$yrow),byrow=TRUE),
               xcol=mask$xcol,yrow=mask$yrow)
  # normalise removes artificial boundary attenuation. A rectangular image
  # supports boundary interpolation; the ORIGINAL polygon still defines W.
  smoothed <- blur(pop_im,sigma=5280,normalise=TRUE,bleed=TRUE)
  lambda <- smoothed * (npoints(X)/integral.im(smoothed,domain=W))
  stopifnot(all(is.finite(lambda[X])),all(lambda[X]>0),
            abs(integral.im(lambda,domain=W)-62)<1e-7)
  saveRDS(lambda,file.path(proc,'lambda_im.Rds'))
  lambda_values <- lambda[as.ppp(as.matrix(grid[inside,]),W)]
  write.csv(data.frame(grid[inside,],population_density_ft2=raw[inside],
    lambda_markets_ft2=lambda_values,lambda_markets_sqmi=lambda_values*5280^2),
    file.path(proc,'lambda_grid.csv'),row.names=FALSE)

  # Use the same dense radius grid as B, report ten 2500-ft upper band edges.
  # Population-proportional Poisson reference, expected count 62; count varies.
  # Fixed externally specified intensity: no additional Kinhom renormalization.
  r <- seq(0,25000,100); bands <- seq(0,25000,2500); idx <- match(bands,r)
  estimate <- function(Y) {
    lv <- lambda[Y]
    stopifnot(all(is.finite(lv)),all(lv>0))
    as.data.frame(Kinhom(Y,lambda=lv,r=r,correction='isotropic',
                        renormalise=FALSE))$iso
  }
  obs <- estimate(X)
  nsim <- 99L; set.seed(671)
  sims <- matrix(NA_real_,length(r),nsim); sim_n <- integer(nsim)
  for(i in seq_len(nsim)) {
    # Function-based thinning uses the exact polygon, avoiding pixel-edge
    # leakage from an image window in rpoispp's faster raster algorithm.
    Y <- rpoispp(as.function(lambda),lmax=max(lambda$v,na.rm=TRUE),win=W)
    stopifnot(all(inside.owin(Y$x,Y$y,W)))
    sim_n[i] <- npoints(Y); sims[,i] <- estimate(Y)
  }
  lo <- apply(sims,1,min); hi <- apply(sims,1,max)
  all_results <- data.frame(distance_ft=r,Kinhom_obs=obs,K_theory=pi*r^2,
    envelope_lo=lo,envelope_hi=hi,L_inhom=sqrt(pmax(obs,0)/pi)-r,
    L_lo=sqrt(pmax(lo,0)/pi)-r,L_hi=sqrt(pmax(hi,0)/pi)-r,
    p_cluster=(1+rowSums(sims>=obs))/(nsim+1),
    p_disperse=(1+rowSums(sims<=obs))/(nsim+1),
    envelope_position=ifelse(obs>hi,'above',ifelse(obs<lo,'below','inside')))
  stopifnot(all(is.finite(as.matrix(all_results[,2:10]))))
  write.csv(all_results[idx,],file.path(out,'table_04_kinhom_values.csv'),row.names=FALSE)
  write.csv(all_results,file.path(proc,'kinhom_results.csv'),row.names=FALSE)
  saveRDS(list(r=r,K=sims,counts=sim_n,seed=671,nsim=nsim),file.path(proc,'simulations.Rds'))

  # ZIP counts use geometry, not the market's postal-address string.
  memberships <- st_intersects(markets,z)
  stopifnot(all(lengths(memberships)==1))
  z$market_count <- tabulate(unlist(memberships),nbins=nrow(z))
  z$income_plot <- ifelse(z$MedIncome>0,z$MedIncome,NA_real_)
  zs <- st_drop_geometry(z)[,c('ZipCode','Pop2000','MedIncome','area_ft2','market_count')]
  write.csv(zs,file.path(proc,'zip_market_summary.csv'),row.names=FALSE)
  policy <- data.frame(field=c('ACCEPT_SNA','EBT_MACHIN','EBT_MACHIN','EBT_MACHIN'),
     category=c('Yes','Yes','Yes*','No'),count=c(sum(markets$ACCEPT_SNA=='Yes'),
       sum(markets$EBT_MACHIN=='Yes'),sum(markets$EBT_MACHIN=='Yes*'),sum(markets$EBT_MACHIN=='No')))
  write.csv(policy,file.path(proc,'payment_attributes.csv'),row.names=FALSE)
  map_theme <- theme_void(base_size=12)+theme(plot.title=element_text(face='bold',size=15),
    plot.subtitle=element_text(size=11),plot.caption=element_text(size=9,hjust=0),
    legend.position='right',plot.margin=margin(15,15,15,15))
  g <- data.frame(grid[inside,],value=lambda_values*5280^2)
  intensity_map <- ggplot()+geom_raster(data=g,aes(x_ft,y_ft,fill=value))+
    geom_sf(data=city,fill=NA,color='grey40',linewidth=.3)+
    geom_sf(data=markets,shape=21,fill='white',color='#923922',size=1.6,stroke=.6)+
    scale_fill_viridis_c(name='Expected markets\nper square mile',option='C')+
    coord_sf(crs=st_crs(city),datum=NA)+labs(title='Population-based reference intensity',
      subtitle='2000 ZIP population | 250-ft grid | 1-mile Gaussian smoothing',
      caption='Circles: 62 market records (2013); one coincident pair. Integral of intensity = 62.\nSurface is a modeling assumption, not observed demand. North is up.')+map_theme
  ggsave(file.path(out,'fig_07_lambda_intensity_map.png'),intensity_map,width=7.3,height=7.4,dpi=180,bg='white')
  labels <- z[z$ZipCode %in% c(19111,19114,19115,19116,19149,19152,19154,19121,19133,19103,19112),]
  label_points <- suppressWarnings(st_point_on_surface(labels))
  income_map <- ggplot()+geom_sf(data=z,aes(fill=income_plot),color='white',linewidth=.35)+
    scale_fill_viridis_c(name='Median household\nincome (2000 USD)',option='D',
      labels=function(x) paste0('$',format(x,big.mark=',',trim=TRUE)),na.value='grey80')+
    geom_sf(data=markets,shape=21,color='#84351F',fill='white',size=2,stroke=.7)+
    geom_sf_text(data=label_points,aes(label=ZipCode),size=2.6,color='black',check_overlap=TRUE)+
    coord_sf(datum=NA)+labs(title='Farmers markets and neighborhood income',
      subtitle='2013 market locations over supplied ZIP-level median household income',
      caption='White circles: market records; the two Reading Terminal records overlap.\nGray ZIP 19112: recorded income is 0, displayed as unavailable. Selected ZIPs labeled; north is up.')+map_theme
  ggsave(file.path(out,'fig_09_income_map_markets.png'),income_map,width=7.3,height=7.4,dpi=180,bg='white')
  png(file.path(out,'fig_08_kinhom_curves.png'),width=1600,height=780,res=150)
  par(mfrow=c(1,2),mar=c(4.4,4.8,3.5,1),mgp=c(2.8,1,0))
  for(center in c(FALSE,TRUE)) {
    v <- if(center) all_results$L_inhom else obs/1e9
    lower <- if(center) all_results$L_lo else lo/1e9
    upper <- if(center) all_results$L_hi else hi/1e9
    theory <- if(center) rep(0,length(r)) else pi*r^2/1e9
    plot(r,v,type='n',ylim=range(v,lower,upper,theory),xlab='Distance (feet)',
      ylab=if(center) 'Centered L inhom (feet)' else 'K inhom (billion square feet)',
      main=if(center) 'Population-adjusted centered L' else 'Population-adjusted K')
    polygon(c(r,rev(r)),c(lower,rev(upper)),col='#CFDDE4',border=NA)
    lines(r,theory,lty=3,col='grey25',lwd=1.5)
    lines(r,v,col='#AD4025',lwd=2); points(bands,v[idx],pch=16,col='#AD4025',cex=.55)
    legend('topleft',c('Observed','99 Poisson simulations (min/max)','Poisson reference'),
      col=c('#AD4025','#CFDDE4','grey25'),lty=c(1,1,3),lwd=c(2,6,1.5),bty='n',cex=.72)
  }
  dev.off()
  summary <- c('HW4 Part C: population-based inhomogeneous K',
    '62 observed records; original Philadelphia polygon; EPSG:2272 (US survey feet).',
    '250-ft raster; sigma 5280 ft; boundary-normalized Gaussian smoothing.',
    sprintf('Grid cells without a ZIP match, filled from nearest ZIP: %d',sum(missing)),
    sprintf('Integrated lambda: %.8f; point-intensity range: %.4g to %.4g per ft2',
      integral.im(lambda,domain=W),min(lambda[X]),max(lambda[X])),
    sprintf('99 plug-in Poisson simulations, seed 671; count range %d-%d.',min(sim_n),max(sim_n)),
    'Isotropic correction; renormalise=FALSE; min/max pointwise envelope, not simultaneous.',
    paste('Above envelope at positive band edges (ft):',paste(bands[bands>0 & all_results$envelope_position[idx]=='above'],collapse=', ')),
    paste('Below envelope at positive band edges (ft):',paste(bands[bands>0 & all_results$envelope_position[idx]=='below'],collapse=', ')),
    'Income zero in 19112 displayed as unavailable; zero population retained before smoothing.',
    'No causal inference, bandwidth validation, or direct measure of travel access is claimed.')
  writeLines(summary,file.path(out,'analysis_handoff_summary_kinhom.txt'))
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info_kinhom.txt'))
  cat(paste(summary,collapse='\n'),'\n')
  invisible(list(table=all_results[idx,],zip=zs,lambda=lambda))
}
if(sys.nframe()==0L) run_hw4_part_c()
