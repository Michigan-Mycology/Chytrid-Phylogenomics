library(tidyverse)

args = commandArgs(trailingOnly = F)
files = list.files(args[1])

for (f in files) {

  f_out = gsub("[.]tsv",".pdf", f)
  strain = gsub("[_]snps[_]in[_]windows[.]tsv", "", f)

  counts = read_delim(f, delim="\t", col_names = F) %>%
    rename(contig = X1, window_start = X2, window_end = X3, num_snps = X4, contig_start_pos = X5) %>%
    mutate(start_pos_cumm = contig_start_pos + window_start) %>%
    mutate(end_pos_cumm = contig_start_pos + window_end) %>%
    mutate(mean_pos_cumm = (end_pos_cumm - ((end_pos_cumm-start_pos_cumm)/2) )) %>%
    select(-window_start, -window_end, -contig_start_pos)

  original_order = unique(counts$contig)
  counts = counts %>%
    mutate(contig=factor(contig, levels=original_order))

  contig_ranges = counts %>%
    group_by(contig) %>%
    summarise(start = min(start_pos_cumm), end = max(end_pos_cumm), max_snps = max(num_snps)) %>%
    mutate(y_max = max(max_snps))
  stripes = rep(c(0,1), ceiling(dim(contig_ranges)[1]/2))
  contig_ranges = contig_ranges %>%
    mutate(to_fill = as.factor(stripes[1:dim(contig_ranges)[1]])) %>%
    mutate(middle = end - ((end-start)/2))

  plt = ggplot(counts, aes(x=mean_pos_cumm, y=num_snps)) +
    geom_rect(data = contig_ranges, aes(xmin = start, xmax = end, ymin = 0, ymax = y_max+10, fill = to_fill), alpha = 0.5, inherit.aes = F) +
    scale_fill_manual(values=c("white", "grey")) +
    geom_line(alpha=0.5) +
    geom_text(data = contig_ranges, aes(label=contig, x=middle, y=y_max+5), inherit.aes = F, size=3, angle=60) +
    geom_point(size=0.5) +
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    guides(fill=F)+
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
    ggtitle(strain)

  ggsave(filename = f_out, plot = plt, device="pdf", width=11, height=8.5, units="in")
}
