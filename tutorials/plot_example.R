require(plotly)
require(htmlwidgets)
require(reticulate)

np <- import ('numpy')
pd <- import("pandas")


unbiass <- pd$read_pickle("github_repos/dingo/tutorials/sars_samples/unbiased_sars_model_samples.pckl")
biomass <- pd$read_pickle("github_repos/dingo/tutorials/sars_samples/max_human_biomass_sars_model_samples.pckl")
vbof <- pd$read_pickle("github_repos/dingo/tutorials/sars_samples/max_vbof_sars_model_samples.pckl")


# In python index we add +1 
vbof_index    <- 3394
biomass_index <- 3393
dgk_1_index   <- 1201


# Find reaction indexes of interest; for example when the differenc in the mean flux is greater than an order of magnitude
for (i in 1: 3394)
{
  x = mean(biomass[i,])
  y = mean(vbof[i,])
  if (x > 0 && y >0) 
  {
    min_val = min(x,y)
    dif = abs(x-y)
    if (dif > 10*min_val)
    {
      print(i)
    }    
  }
}


## plot style ##
m = 10
axx <- list(
  title = '% GK1',
  ticketmode = 'array',
  ticktext = c("0", "0.5", "1"),
  tickvals = c(0, m/2, m-1)
)
axy <- list(
  title = '% VBOF',
  ticketmode = 'array',
  ticktext = c("0", "0.5", "1"),
  tickvals = c(0,m/2,m-1)
)
axz <- list(
  title = 'probability',
  ticketmode = 'array',
  ticktext = c(" "),
  tickvals = c(0)
)


# UNBIASED CASE
#--------------

unb_gk1_d  = unbiass[dgk_1_index, ]
unb_vbof_d = unbiass[vbof_index, ]
unb_human_biom_d = unbiass[biomass_index, ]

# VBOF - human biomass
axx[["title"]] <- "% VBOF"; 
axy[["title"]] <- "% human_biomass";
cop_vbof_human_biom = compute_copula(unb_vbof_d, unb_human_biom_d, m)
title_char = "Unbiased case, VBOF ~ human biomass"
fig = plotly::plot_ly(z = ~cop_vbof_human_biom) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title=title_char, 
                     font=list(size=16), 
                     scene = list(xaxis=axx, yaxis=axy, zaxis=axz),
                     autosize = F, 
                     height = 900,
                     margin = list(l=50, r=50, b=100, t=100, pad=4)
                     )
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/tutorials/sars_copulas/",title_char,".html"))

# VBOF- GK1
axx[["title"]] <- "% VBOF"; 
axy[["title"]] <- "% GK1";
cop_gk1_vbof = compute_copula( unb_vbof_d, unb_gk1_d, m )
title_char = "Unbiased case, VBOF ~ GK1"
fig = plotly::plot_ly(z = ~cop_gk1_vbof) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title=title_char, 
                     font=list(size=16), 
                     scene = list(xaxis=axx, yaxis=axy, zaxis=axz),
                     autosize = F, 
                     height = 900,
                     margin = list(l=50, r=50, b=100, t=100, pad=4)
)
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/tutorials/sars_copulas/",title_char,".html"))


# human biomass - GK1 
axx[["title"]] <- "% human_biomass"; 
axy[["title"]] <- "% GK1";
cop_gk1_vbof = compute_copula( unb_human_biom_d, unb_gk1_d, m )
title_char = "Unbiased case, human biomass ~ GK1"
fig = plotly::plot_ly(z = ~cop_gk1_vbof) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title=title_char, 
                     font=list(size=16), 
                     scene = list(xaxis=axx, yaxis=axy, zaxis=axz),
                     autosize = F, 
                     height = 900,
                     margin = list(l=50, r=50, b=100, t=100, pad=4)
)
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/tutorials/sars_copulas/",title_char,".html"))


# AFTER MAXIMIZING FOR VBOF
#--------------------------

max_vbof_gk1_d  = vbof[dgk_1_index, ]
max_vbof_human_biom_d = vbof[biomass_index, ]

# human biomass - GK1 
axx[["title"]] <- "% human_biomass"; 
axy[["title"]] <- "% GK1";
cop_gk1_vbof = compute_copula( max_vbof_human_biom_d, max_vbof_gk1_d, m )
title_char = "After maximizing VBOF, human biomass ~ GK1"
fig = plotly::plot_ly(z = ~cop_gk1_vbof) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title=title_char, 
                     font=list(size=16), 
                     scene = list(xaxis=axx, yaxis=axy, zaxis=axz),
                     autosize = F, 
                     height = 900,
                     margin = list(l=50, r=50, b=100, t=100, pad=4)
)
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/tutorials/sars_copulas/",title_char,".html"))


# AFTER MAXIMIZING FOR HUMAN BIOMASS 
#-----------------------------------

max_hb_gk1_d = biomass[dgk_1_index, ]
max_hb_vbof_d = biomass[vbof_index, ]

# VBOF- GK1 
axx[["title"]] <- "% VBOF"; 
axy[["title"]] <- "% GK1";
cop_gk1_vbof = compute_copula( max_hb_vbof_d, max_hb_gk1_d, m )
title_char = "After maximizing human biomass, VBOF ~ GK1"
fig = plotly::plot_ly(z = ~cop_gk1_vbof) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title=title_char, 
                     font=list(size=16), 
                     scene = list(xaxis=axx, yaxis=axy, zaxis=axz),
                     autosize = F, 
                     height = 900,
                     margin = list(l=50, r=50, b=100, t=100, pad=4)
)
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/tutorials/sars_copulas/",title_char,".html"))

#-----------------------------------------------------------------------

# nets = unbiass[biomass_index,]
# mar = compute_marginal(unb_vbof_d, unb_human_biom_d, unb_gk1_d, m)
# 
# h1 = hist(mar,     
#           main="",     
#           xlab="Flux (mmol/gDW/h)",     
#           border="black",     
#           col="pink",     
#           xlim=c(min(mar), max(mar)),     
#           las=1,     
#           breaks=10,
#           prob = TRUE) 

