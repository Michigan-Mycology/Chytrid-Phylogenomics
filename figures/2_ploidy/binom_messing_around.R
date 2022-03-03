min = 0
max = 1000000
coverage = 11
step = 1

x = seq(min, max, step)

#y = dbinom(x, coverage, 0.5)

### Normal estimation
dist_mean = coverage*0.5
dist_sd = sqrt(coverage*0.5*0.5)
y = dnorm(x = x, mean = dist_mean, sd = dist_sd)

df = as_tibble(cbind(x,y)) %>%
  mutate(x = (x*(max/coverage))/max)

#raw_df = as_tibble(cbind(x,y))

stdev = (sqrt(coverage*0.5*(1-0.5)))
scaled_stdev = (stdev*(max/coverage))/max

pbiom()
snp_fit_range_unscaled = c( (coverage/2)-stdev, (coverage/2)+stdev )
snp_fit_range = c(0.5-scaled_stdev, 0.5+scaled_stdev)
snp_fit_range

#p_outside_lower_bound = pbinom(, 100, 0.5)
#p_outside_lower_bound

ggplot(df, aes(x=x, y=y)) +
  geom_point() +
  geom_vline(xintercept=snp_fit_range[1]) +
  geom_vline(xintercept=snp_fit_range[2]) +
  xlim(0,1) +
  theme_bw()

prob = 0
for (i in 0:1) {
  this = choose(100,i)*(0.5^i)*(0.5^(100-i))
  print(this)
  prob = prob + this
}
prob
