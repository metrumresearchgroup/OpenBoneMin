library(tidyverse)
library(mrgsolve)
 
#' converts PO4 mM to phosphorous mg/dL
convPhos <- 3.095975232198142

#' for mass conversion: 1 mmol phosphorus = 31 mg phosphorus
convCtriol <- 34.61538461538461 / 90 # http://www.endmemo.com/medical/unitconvert/Vitamin_D.php

#' Regimens
reg <- tribble(
  ~n, ~dose, ~freq, ~trt, 
  1,  16, 24, "16 mg QD", 
  2,  32, 24, "32 mg QD", 
  3,   8, 12, "8 mg BID", 
  4,  16, 12, "16 mg BID",
  5,  16, 24, "16 mg QD (5/2 on/off)", 
  6,  32, 24, "32 mg QD (5/2 on/off)", 
  7,   8, 12, "8 mg BID (5/2 on/off)", 
  8,  16, 12, "16 mg BID (5/2 on/off)",
  9,   8, 12, "8 mg BID (4/3 on/off)",
 10,  16, 12, "16 mg BID (4/3 on/off)",
) %>% mutate(tdd = dose * 24 / freq, trt = fct_inorder(trt))
reg <- mutate(reg, TDD = factor(paste(tdd, " mg per day")))

#' once and twice daily doses
qd <- expand.ev(amt=c(16,32), ii=24, addl = 27) %>% mutate(n = c(1,2))
bid <- mutate(qd, amt = c(8,16), ii = 12, addl = 55, n = c(3,4))

#' once daily doses, 5 on / 2 off
qd_52_16 <- ev_days(ev(amt = 16, n = 5), ii = 7*24, addl = 3, n = 9, days = "m,t,w,th,f")
qd_52_32 <- mutate(qd_52_16, amt = 32, n = 6)

#' twice daily doses, 5 on / 2 off
am <- ev_days(ev(amt = 8, n = 7), ii = 7*24, addl = 3, days = "m,t,w,th,f")
pm <- mutate(am, time = time + 12)
bid_52_8 <- bind_rows(am,pm)
bid_52_16 <- mutate(bid_52_8, amt = 16, n = 8)

#' twice daily doses, 4 on / 3 off
am <- ev_days(ev(amt = 8, n = 9), ii = 7*24, addl = 3, days = "m,t,w,th")
pm <- mutate(am, time = time + 12)
bid_43_8 <- bind_rows(am,pm) 
bid_43_16 <- mutate(bid_43_8, amt = 16, n = 10)

#' combined data 
data <- bind_rows(qd,bid,bid_52_8,bid_52_16,bid_43_8,bid_43_16,qd_52_16,qd_52_32) 
data <- arrange(data,n,time) %>% mutate(rate = -2, ID = n)

#' Simulate
mod <- mread(model="asp5878_SysPcolModel", end = 28*24, delta = 0.1)
out <- mrgsim_d(mod, data, carry_out = "n", output = "df")

sims <- left_join(out, reg, by = "n") %>% mutate(Phos = ECCPhos/14*convPhos)

#' Plot
ggplot(data = sims) +
  geom_line(aes(color=trt,x=time,y=Phos)) +
  geom_hline(yintercept = 7.0, lty = 2) + 
  geom_hline(yintercept = 5.5, lty = 3) +
  geom_hline(yintercept = 6.0, lty = 4) +
  labs(x="Time from first dose (hours)", y="Phosphate (mg/dL)") +
  facet_wrap(~TDD+trt,nrow=2) + theme(legend.position ='none' )

ggsave("fgf-phos.png",width = 8, height = 4, dpi = 300)

