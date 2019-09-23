library(tidyverse)
df <- tibble(
  Month = rep(c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun' ,'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), 3),
  Year = c(rep(2018, 12), rep(2019, 12), rep(2020, 12)),
  Events = NA)
df


df <- tribble(~Event, ~Start_date, ~End_date,
        "Tissue collection & cell sorting", "2017-01-01", "2018-12-31",
        "Array prep", "2019-01-01", "2019-04-22",
        "Array run", "2019-04-22", "2019-05-03",
        "Data QC", "2019-05-03", "2019-06-01",
        "Analysis", "2019-06-01", "2020-12-30")
df %>%
  mutate_at(c('Start_date', 'End_date'), as.Date)%>%
  ggplot() +
  geom_segment(aes(x = Start_date, xend = End_date, y= 0, yend = 0, color = Event), size = 10) +
  scale_y_continuous(limits = c(-0.0005,0.0005)) +
  scale_x_date()
