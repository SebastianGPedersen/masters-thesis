a_m_func <- function(m) {sqrt(2*log(m))}
b_m_func <- function(m) {a_m_func(m) - 1/2*log(pi*log(m))/a_m_func(m)}

qgumbel <- function(p) {-log(-log(p))} #from: https://stats.stackexchange.com/questions/71197/usable-estimators-for-parameters-in-gumbel-distribution

q95 <- function(m) {qgumbel(p=0.95)/a_m_func(m)+b_m_func(m)}

#q95(60*6.5) = 3.89. Dvs. sætter bandwidth ratio til 10 hvis hvert minut
#q95(12*60*6.5) = 4.43. Dvs. sætter bandwidth ratio til 15 hvis hvert 5. sekund
